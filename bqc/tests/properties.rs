// Copyright (c) 2026 The Hiller Lab at the Senckenberg Gesellschaft für Naturforschung
// Distributed under the terms of the GNU General Public License, Version 3.0.

//! Property tests for the invariants the design treats as contractual.

#![allow(clippy::float_cmp)]

mod common;

use bqc::adapter::{Adapter, AdapterParams, AdapterStage, find_three_prime};
use bqc::filter::{FilterStage, complexity};
use bqc::io::Schema;
use bqc::process::Workflow;
use bqc::read::{ReadView, Span};
use bqc::segment::{SegmentScratch, SegmentStage, Terminal, segment};
use bqc::trim::{MateTrim, PolyParams, QualityCut, TrimStage};
use common::{Record, read_cbq, write_cbq};
use proptest::prelude::*;

fn sequence() -> impl Strategy<Value = Vec<u8>> {
    prop::collection::vec(
        prop::sample::select(vec![b'A', b'C', b'G', b'T', b'N']),
        0..80,
    )
}

/// A workflow exercising every stage at once.
fn everything() -> Workflow {
    let adapter = AdapterStage::new(
        vec![Adapter::new("a", b"AGATCGGAAGAGCACACGT").unwrap()],
        vec![Adapter::new("b", b"AGATCGGAAGAGCGTCGTG").unwrap()],
        AdapterParams {
            min_overlap: 6,
            max_error_rate: 0.2,
            max_errors: None,
            allow_indels: false,
        },
        None,
    )
    .unwrap();
    let mate = MateTrim {
        front: 2,
        tail: 1,
        quality_front: Some(QualityCut {
            minimum_phred: 20,
            window: 4,
        }),
        quality_tail: Some(QualityCut {
            minimum_phred: 20,
            window: 4,
        }),
        terminal_n: true,
        poly_g: Some(PolyParams {
            min_length: 5,
            max_mismatch_rate: 0.2,
        }),
        poly_x: Some(PolyParams {
            min_length: 5,
            max_mismatch_rate: 0.2,
        }),
        max_length: Some(30),
        ..MateTrim::default()
    };
    let filter = FilterStage {
        min_length: Some(5),
        max_length: Some(60),
        max_n: Some(3),
        max_n_fraction: Some(0.5),
        qualified_quality: 15,
        max_unqualified_bases: Some(10),
        max_unqualified_fraction: Some(0.5),
        min_mean_quality: Some(15),
        min_complexity: Some(0.2),
    };
    Workflow::new(
        Some(adapter),
        None,
        Some(trim_stage(mate)),
        Some(filter),
        None,
    )
    .unwrap()
}

fn trim_stage(mate: MateTrim) -> TrimStage {
    TrimStage { r1: mate, r2: mate }
}

proptest! {
    // Integration tests cannot resolve a persistence path for regressions.
    #![proptest_config(ProptestConfig { failure_persistence: None, ..ProptestConfig::default() })]

    /// Spans stay inside the original read and never grow it.
    #[test]
    fn spans_are_always_valid_sub_ranges(
        sequence in sequence(),
    ) {
        let quality = vec![b'5'; sequence.len()];
        let read = ReadView::new(&sequence, Some(&quality), 0, "R1").unwrap();
        let result = everything().process(0, read, Some(read)).unwrap();
        for (_, mate) in result.mates() {
            let span = mate.retained;
            prop_assert!(span.start <= span.end);
            prop_assert!(span.end <= sequence.len());
            prop_assert!(span.len() <= sequence.len());
            prop_assert_eq!(mate.original_length, sequence.len());
            // Sequence and quality are sliced with identical coordinates.
            prop_assert_eq!(
                span.sequence(read).len(),
                span.quality(read).unwrap().len()
            );
            // The retained bases are a contiguous window of the original.
            prop_assert_eq!(span.sequence(read), &sequence[span.start..span.end]);
        }
    }

    /// Trimming and adapter removal never add or reorder bases.
    #[test]
    fn transformations_only_remove_bases(
        sequence in sequence(),
        seed in 0u8..64,
    ) {
        let quality: Vec<u8> = (0..sequence.len())
            .map(|i| 33 + ((i as u8).wrapping_add(seed) % 41))
            .collect();
        let read = ReadView::new(&sequence, Some(&quality), 0, "R1").unwrap();
        let result = everything().process(7, read, None).unwrap();
        let retained = result.r1.retained.sequence(read);
        prop_assert!(retained.len() <= sequence.len());
        prop_assert!(
            sequence.windows(retained.len().max(1)).any(|window| window == retained)
                || retained.is_empty()
        );
        // Stage accounting is consistent.
        prop_assert!(result.r1.adapter_trimmed_length <= result.r1.original_length);
        prop_assert!(result.r1.final_length() <= result.r1.adapter_trimmed_length);
        prop_assert_eq!(
            result.r1.adapter_trimmed_length - result.r1.final_length(),
            result.r1.trimming.total() as usize
        );
    }

    /// An adapter hit always points at a coordinate that can be trimmed to.
    #[test]
    fn adapter_hits_are_within_the_read(
        sequence in sequence(),
        adapter in prop::collection::vec(prop::sample::select(vec![b'A', b'C', b'G', b'T']), 6..20),
        min_overlap in 1usize..8,
        rate in 0.0f64..0.5,
    ) {
        let adapters = vec![Adapter::new("a", &adapter).unwrap()];
        let params = AdapterParams { min_overlap, max_error_rate: rate, max_errors: None, allow_indels: false };
        if let Some(hit) = find_three_prime(&adapters, params, &sequence) {
            prop_assert!(hit.start < sequence.len());
            prop_assert!(hit.overlap >= min_overlap);
            prop_assert!(hit.start + hit.overlap <= sequence.len());
            prop_assert!(hit.errors <= params.error_limit(hit.overlap));
            // The reported overlap is the largest possible at that coordinate.
            prop_assert_eq!(hit.overlap, adapter.len().min(sequence.len() - hit.start));
        }
    }

    /// Indel-aware hits obey the same coordinate and budget invariants. The
    /// overlap counts adapter bases, so it need not fit inside the read.
    #[test]
    fn indel_hits_obey_the_same_invariants(
        sequence in sequence(),
        adapter in prop::collection::vec(prop::sample::select(vec![b'A', b'C', b'G', b'T']), 6..20),
        min_overlap in 1usize..8,
        rate in 0.0f64..0.5,
    ) {
        let adapters = vec![Adapter::new("a", &adapter).unwrap()];
        let params = AdapterParams { min_overlap, max_error_rate: rate, max_errors: None, allow_indels: true };
        if let Some(hit) = find_three_prime(&adapters, params, &sequence) {
            prop_assert!(hit.start < sequence.len());
            prop_assert!(hit.overlap >= min_overlap);
            prop_assert!(hit.overlap <= adapter.len());
            prop_assert!(hit.errors <= params.error_limit(hit.overlap));
        }
    }

    /// Fragments partition the read: ordered, disjoint, inside the read, and
    /// never covering a base an accepted delimiter consumed.
    #[test]
    fn fragments_partition_the_read_around_the_delimiters(
        sequence in sequence(),
        adapter in prop::collection::vec(prop::sample::select(vec![b'A', b'C', b'G', b'T']), 6..12),
        min_overlap in 1usize..7,
        rate in 0.0f64..0.4,
        min_segment_length in 0usize..5,
        max_segments in 1usize..6,
        discard in any::<bool>(),
    ) {
        let params = AdapterParams { min_overlap, max_error_rate: rate, max_errors: None, allow_indels: false };
        let stage = SegmentStage::new(
            vec![Adapter::new("a", &adapter).unwrap()],
            params,
            if discard { Terminal::Discard } else { Terminal::Keep },
            min_segment_length,
            max_segments,
        );
        // A delimiter shorter than the minimum overlap can never match, and the
        // stage rejects that configuration rather than silently doing nothing.
        let Ok(stage) = stage else {
            prop_assert!(adapter.len() < min_overlap);
            return Ok(());
        };
        let mut scratch = SegmentScratch::default();
        let outcome = segment(&stage, &sequence, &mut scratch);

        prop_assert!(scratch.fragments.len() <= max_segments);
        prop_assert_eq!(scratch.fragments.len(), outcome.emitted());
        let mut previous_end = 0usize;
        for (index, fragment) in scratch.fragments.iter().enumerate() {
            prop_assert_eq!(fragment.index, index, "indices are dense and ordered");
            prop_assert!(fragment.span.start >= previous_end, "fragments are disjoint");
            prop_assert!(fragment.span.start < fragment.span.end, "no empty fragment");
            prop_assert!(fragment.span.end <= sequence.len(), "inside the read");
            prop_assert!(fragment.span.len() >= min_segment_length);
            if discard {
                prop_assert!(fragment.internal(), "terminals were discarded");
            }
            previous_end = fragment.span.end;
        }
        // Every accepted delimiter's own bases lie outside every fragment.
        for hit in scratch.delimiters() {
            let delimiter = hit.start..hit.start + hit.consumed;
            for fragment in &scratch.fragments {
                prop_assert!(
                    fragment.span.end <= delimiter.start || fragment.span.start >= delimiter.end,
                    "fragment {:?} overlaps delimiter {:?}", fragment.span, delimiter
                );
            }
        }
        // Every candidate was either accepted as a delimiter or suppressed.
        prop_assert_eq!(outcome.candidates, outcome.boundaries + outcome.suppressed);
    }

    /// Reported fragment counts reconcile with the records and sidecar rows
    /// actually written, and every fragment maps to exactly one source record.
    #[test]
    fn segment_counts_reconcile_with_the_written_output(
        lengths in prop::collection::vec(1usize..50, 1..25),
        threads in 1usize..4,
        min_segment_length in 1usize..6,
    ) {
        const DELIMITER: &str = "AGGTCAGTCTACAT";
        let schema = Schema { paired: false, quality: true, headers: true, flags: true };
        let dir = tempfile::tempdir().unwrap();
        let input = dir.path().join("in.cbq");
        let records: Vec<Record> = lengths
            .iter()
            .enumerate()
            .map(|(index, &length)| {
                let filler: Vec<u8> = (0..length).map(|i| b"ACGTN"[(i + index) % 5]).collect();
                // Every other read carries a delimiter, so both the split and the
                // whole-read paths are exercised.
                let read = if index % 2 == 0 {
                    [filler.as_slice(), DELIMITER.as_bytes(), &filler].concat()
                } else {
                    filler
                };
                let quality: Vec<u8> = (0..read.len()).map(|i| 33 + ((i * 7 + index) % 41) as u8).collect();
                Record::new(schema, index, &read, &read).with_quality(&quality, &quality)
            })
            .collect();
        write_cbq(&input, schema, &records, 512);

        let output = dir.path().join("fragments.cbq");
        let sidecar = dir.path().join("segments.tsv");
        let report = dir.path().join("report.json");
        common::bqc_ok(&[
            "segment",
            input.to_str().unwrap(),
            "-o", output.to_str().unwrap(),
            "--adapter-r1", DELIMITER,
            "--segments", sidecar.to_str().unwrap(),
            "--report", report.to_str().unwrap(),
            "--min-segment-length", &min_segment_length.to_string(),
            "-T", &threads.to_string(),
        ]);

        let (out_schema, fragments) = read_cbq(&output);
        prop_assert_eq!(out_schema, schema);
        let rows: Vec<String> = std::fs::read_to_string(&sidecar).unwrap()
            .lines().skip(1).map(str::to_string).collect();
        let json: serde_json::Value =
            serde_json::from_str(&std::fs::read_to_string(&report).unwrap()).unwrap();
        let segment = &json["segment"];

        // With no filter stage every emitted fragment is written.
        prop_assert_eq!(segment["fragments_emitted"].as_u64().unwrap(), fragments.len() as u64);
        prop_assert_eq!(rows.len(), fragments.len(), "one sidecar row per fragment");
        prop_assert_eq!(segment["source_records"].as_u64().unwrap(), records.len() as u64);
        prop_assert_eq!(
            segment["internal_fragments"].as_u64().unwrap()
                + segment["terminal_fragments"].as_u64().unwrap(),
            fragments.len() as u64
        );

        let mut seen = vec![0usize; records.len()];
        let mut previous: Option<(u64, u64)> = None;
        for (row, fragment) in rows.iter().zip(&fragments) {
            let columns: Vec<&str> = row.split('\t').collect();
            let source: u64 = columns[0].parse().unwrap();
            let index: u64 = columns[1].parse().unwrap();
            let (start, end): (usize, usize) =
                (columns[3].parse().unwrap(), columns[4].parse().unwrap());
            // Exactly one source record, and the ordering contract.
            prop_assert!((source as usize) < records.len());
            seen[source as usize] += 1;
            if let Some(last) = previous {
                prop_assert!((source, index) > last, "output must be strictly ordered");
            }
            previous = Some((source, index));
            // Sequence and quality are sliced with identical coordinates.
            let record = &records[source as usize];
            prop_assert_eq!(&fragment.s_seq, &record.s_seq[start..end].to_vec());
            prop_assert_eq!(
                fragment.s_qual.as_ref().unwrap(),
                &record.s_qual.as_ref().unwrap()[start..end].to_vec()
            );
            prop_assert!(end - start >= min_segment_length);
        }
        // A read whose every fragment was dropped contributes none; no read may
        // contribute more fragments than the report counted for it.
        prop_assert_eq!(seen.iter().sum::<usize>(), fragments.len());
        prop_assert!(seen.iter().all(|&count| count <= segment["max_fragments_per_source"].as_u64().unwrap() as usize));
    }

    /// Complexity is a fraction, and constant reads are minimally complex.
    #[test]
    fn complexity_is_bounded(sequence in sequence()) {
        let value = complexity(&sequence);
        prop_assert!((0.0..=1.0).contains(&value));
        if sequence.len() >= 2 && sequence.iter().all(|&b| b == sequence[0]) {
            prop_assert_eq!(value, 0.0);
        }
    }

    /// Filtering never depends on the order in which reasons are collected.
    #[test]
    fn filter_reasons_are_deterministic(
        sequence in sequence(),
        quality_seed in 0u8..40,
    ) {
        let quality = vec![33 + quality_seed; sequence.len()];
        let filter = FilterStage {
            min_length: Some(10),
            max_n: Some(1),
            qualified_quality: 20,
            max_unqualified_fraction: Some(0.3),
            min_mean_quality: Some(25),
            min_complexity: Some(0.4),
            ..FilterStage::default()
        };
        let first = filter.evaluate(&sequence, Some(&quality));
        let second = filter.evaluate(&sequence, Some(&quality));
        prop_assert_eq!(first, second);
        prop_assert_eq!(first.label(), second.label());
    }
}

proptest! {
    // Each case writes and reads real CBQ files, so fewer cases are run.
    #![proptest_config(ProptestConfig {
        cases: 24,
        failure_persistence: None,
        ..ProptestConfig::default()
    })]

    /// Accepted and rejected records reconcile with the input, in input order,
    /// with all metadata preserved, for any thread count.
    #[test]
    fn end_to_end_accounting_reconciles(
        lengths in prop::collection::vec(1usize..60, 1..40),
        threads in 1usize..5,
    ) {
        let schema = Schema { paired: true, quality: true, headers: true, flags: true };
        let dir = tempfile::tempdir().unwrap();
        let input = dir.path().join("in.cbq");
        let records: Vec<Record> = lengths
            .iter()
            .enumerate()
            .map(|(index, &length)| {
                let r1: Vec<u8> = (0..length).map(|i| b"ACGTN"[(i + index) % 5]).collect();
                let r2: Vec<u8> = (0..length).map(|i| b"ACGTN"[(i + index + 2) % 5]).collect();
                let q1: Vec<u8> = (0..length).map(|i| 33 + ((i * 7 + index) % 41) as u8).collect();
                let q2 = q1.clone();
                Record::new(schema, index, &r1, &r2).with_quality(&q1, &q2)
            })
            .collect();
        write_cbq(&input, schema, &records, 512);

        let accepted = dir.path().join("accepted.cbq");
        let failed = dir.path().join("failed.cbq");
        common::bqc_ok(&[
            "workflow",
            input.to_str().unwrap(),
            "-o",
            accepted.to_str().unwrap(),
            "--failed",
            failed.to_str().unwrap(),
            "--front",
            "1",
            "--quality-tail",
            "20",
            "--min-length",
            "10",
            "-T",
            &threads.to_string(),
        ]);

        let (accepted_schema, accepted_records) = read_cbq(&accepted);
        let (_, failed_records) = read_cbq(&failed);
        // Every input record is accounted for exactly once.
        prop_assert_eq!(accepted_records.len() + failed_records.len(), records.len());
        prop_assert_eq!(accepted_schema, schema);
        // Rejected records are verbatim copies, kept in input order.
        for record in &failed_records {
            prop_assert!(records.contains(record), "failed records must be verbatim copies");
        }
        let rejected_headers: Vec<_> =
            failed_records.iter().map(|r| r.s_header.clone()).collect();
        let expected_order: Vec<_> = records
            .iter()
            .filter(|record| failed_records.contains(record))
            .map(|record| record.s_header.clone())
            .collect();
        prop_assert_eq!(&rejected_headers, &expected_order);
        // Headers, flags and pairing survive on the accepted side.
        for record in &accepted_records {
            let original = records
                .iter()
                .find(|candidate| candidate.s_header == record.s_header)
                .expect("accepted record must derive from an input record");
            prop_assert_eq!(record.flag, original.flag);
            prop_assert!(record.x_seq.is_some());
            prop_assert_eq!(
                record.s_seq.len(),
                record.s_qual.as_ref().unwrap().len()
            );
        }
    }
}

/// A regression guard for span arithmetic on the empty read.
#[test]
fn empty_reads_produce_empty_spans() {
    let read = ReadView::new(b"", Some(b""), 0, "R1").unwrap();
    let result = everything().process(0, read, Some(read)).unwrap();
    assert_eq!(result.r1.retained, Span { start: 0, end: 0 });
    assert!(
        !result.accepted(),
        "--min-length rejects a zero-length read"
    );
}

proptest! {
    #![proptest_config(ProptestConfig {
        cases: 200,
        failure_persistence: None,
        ..ProptestConfig::default()
    })]

    /// Correction only ever touches eligible mismatch positions, and everything
    /// it writes agrees with the donor mate.
    #[test]
    fn correction_only_changes_eligible_positions(
        insert in prop::collection::vec(prop::sample::select(vec![b'A', b'C', b'G', b'T']), 60..90),
        breakages in prop::collection::vec((0usize..40, 0u8..45, 0u8..45), 0..6),
    ) {
        use bqc::correct::{CorrectionScratch, CorrectionStage, plan};
        use bqc::overlap::{OverlapParams, complement, find_overlap};
        use bqc::process::Mate;

        let read_length = insert.len() * 2 / 3;
        let offset = insert.len() - read_length;
        let mut r1 = insert[..read_length].to_vec();
        let r2: Vec<u8> = insert[offset..].iter().rev().map(|&b| complement(b)).collect();
        let mut q1 = vec![b'!' + 40; r1.len()];
        let mut q2 = vec![b'!' + 40; r2.len()];
        // Plant disagreements at arbitrary positions with arbitrary qualities.
        for (position, score1, score2) in &breakages {
            let position = position % read_length;
            r1[position] = match r1[position] { b'A' => b'C', _ => b'A' };
            q1[position] = b'!' + score1;
            let mirrored = position % r2.len();
            q2[mirrored] = b'!' + score2;
        }

        let params = OverlapParams { min_overlap: 20, max_error_rate: 0.5 };
        let Some(overlap) = find_overlap(&r1, &r2, params) else { return Ok(()); };
        let stage = CorrectionStage::default();
        let original1 = ReadView::unchecked(&r1, Some(&q1));
        let original2 = ReadView::unchecked(&r2, Some(&q2));

        let mut scratch = CorrectionScratch::default();
        let summary = plan(stage, original1, original2, overlap, &mut scratch.edits);
        let edits = scratch.edits.clone();
        let (corrected1, corrected2) = scratch.apply_pair(original1, original2);

        // Lengths never change and quality tracks sequence.
        prop_assert_eq!(corrected1.sequence.len(), r1.len());
        prop_assert_eq!(corrected2.sequence.len(), r2.len());
        prop_assert_eq!(corrected1.quality.unwrap().len(), corrected1.sequence.len());
        prop_assert_eq!(corrected2.quality.unwrap().len(), corrected2.sequence.len());

        // Exactly the planned positions changed, and nothing else did.
        for (mate, before, after) in [
            (Mate::R1, (&r1, &q1), (corrected1.sequence, corrected1.quality.unwrap())),
            (Mate::R2, (&r2, &q2), (corrected2.sequence, corrected2.quality.unwrap())),
        ] {
            let planned: Vec<usize> = edits
                .iter()
                .filter(|edit| edit.target == mate)
                .map(|edit| edit.target_position)
                .collect();
            for position in 0..before.0.len() {
                let changed = after.0[position] != before.0[position]
                    || after.1[position] != before.1[position];
                prop_assert_eq!(
                    changed,
                    planned.contains(&position),
                    "{:?} position {} changed without a plan",
                    mate,
                    position
                );
            }
        }

        // Every correction agrees with the donor after reverse complementation,
        // carries the donor's exact quality byte, and is in range.
        for edit in &edits {
            prop_assert!(edit.target_position < before_len(edit.target, r1.len(), r2.len()));
            let (target_sequence, target_quality, donor_sequence, donor_quality) =
                match edit.target {
                    Mate::R1 => (
                        corrected1.sequence,
                        corrected1.quality.unwrap(),
                        corrected2.sequence,
                        corrected2.quality.unwrap(),
                    ),
                    Mate::R2 => (
                        corrected2.sequence,
                        corrected2.quality.unwrap(),
                        corrected1.sequence,
                        corrected1.quality.unwrap(),
                    ),
                };
            prop_assert!(edit.donor_position < donor_sequence.len());
            prop_assert_eq!(
                target_sequence[edit.target_position],
                complement(donor_sequence[edit.donor_position])
            );
            prop_assert_eq!(
                target_quality[edit.target_position],
                donor_quality[edit.donor_position],
                "the donor's raw quality byte is copied"
            );
            // Donor evidence is never ambiguous.
            prop_assert!(matches!(
                donor_sequence[edit.donor_position],
                b'A' | b'C' | b'G' | b'T'
            ));
        }

        // The summary reconciles with the edit list.
        prop_assert_eq!(summary.corrected() as usize, edits.len());
        prop_assert_eq!(
            summary.corrected_r1 as usize,
            edits.iter().filter(|e| e.target == Mate::R1).count()
        );
        prop_assert_eq!(
            summary.corrected_r2 as usize,
            edits.iter().filter(|e| e.target == Mate::R2).count()
        );
        prop_assert!(summary.corrected() + summary.unresolved + summary.skipped_noncanonical
            <= summary.mismatches);
    }
}

/// Length of the mate an edit targets.
fn before_len(mate: bqc::process::Mate, r1: usize, r2: usize) -> usize {
    match mate {
        bqc::process::Mate::R1 => r1,
        bqc::process::Mate::R2 => r2,
    }
}

// ------------------------------------------------------------------- sniffing

use bqc::sniff::sample::{SamplePlan, block_work};
use bqc::sniff::{Confidence, Decision, related};

proptest! {
    /// A sample plan always selects distinct, ascending, in-range indices.
    ///
    /// This is what makes the result reproducible without a seed: the plan is
    /// arithmetic, so it cannot depend on scheduling or on a generator.
    #[test]
    fn a_sample_plan_is_a_valid_selection(
        start in 0u64..10_000,
        length in 0u64..100_000,
        requested in 0u64..5_000,
    ) {
        let plan = SamplePlan::new(start, start + length, requested);
        prop_assert_eq!(plan.selected, requested.min(length));

        let indices: Vec<u64> = plan.indices().collect();
        prop_assert_eq!(indices.len() as u64, plan.selected);
        for pair in indices.windows(2) {
            prop_assert!(pair[0] < pair[1], "indices must strictly ascend");
        }
        for &index in &indices {
            prop_assert!(index >= start, "index below the range");
            prop_assert!(index < start + length, "index past the range");
        }
    }

    /// Evaluating a plan twice gives the same answer.
    #[test]
    fn a_sample_plan_is_reproducible(
        start in 0u64..1_000,
        length in 1u64..50_000,
        requested in 1u64..2_000,
    ) {
        let plan = SamplePlan::new(start, start + length, requested);
        let once: Vec<u64> = plan.indices().collect();
        let again: Vec<u64> = plan.indices().collect();
        prop_assert_eq!(once, again);
    }

    /// Relatedness is reflexive and symmetric.
    ///
    /// It decides whether a second candidate is a rival library or another
    /// spelling of the first, so an asymmetric answer would make the decision
    /// depend on candidate order.
    #[test]
    fn relatedness_is_reflexive_and_symmetric(
        a in prop::collection::vec(prop::sample::select(vec![b'A', b'C', b'G', b'T']), 1..40),
        b in prop::collection::vec(prop::sample::select(vec![b'A', b'C', b'G', b'T']), 1..40),
    ) {
        prop_assert!(related(&a, &a), "a sequence must be related to itself");
        prop_assert_eq!(related(&a, &b), related(&b, &a));
    }
}

/// Every block a plan touches is listed exactly once, and every sampled record
/// is claimed by exactly one block.
#[test]
fn block_work_partitions_the_sample() {
    let fixture = tempfile::tempdir().expect("tempdir");
    let path = fixture.path().join("input.cbq");
    let schema = Schema {
        paired: false,
        quality: true,
        headers: true,
        flags: false,
    };
    let records: Vec<Record> = (0..5000)
        .map(|index| {
            let sequence = vec![b"ACGT"[index % 4]; 60];
            let quality = vec![b'I'; 60];
            Record::new(schema, index, &sequence, b"").with_quality(&quality, b"")
        })
        .collect();
    write_cbq(&path, schema, &records, 1 << 13);

    let input = bqc::io::CbqInput::open(&path).expect("open");
    for requested in [1u64, 7, 100, 999, 5000, 9999] {
        let plan = SamplePlan::for_input(&input, None, requested);
        let work = block_work(&plan, &input);

        // Each block appears at most once, and blocks ascend.
        for pair in work.windows(2) {
            assert!(pair[0].block < pair[1].block, "blocks repeat or descend");
        }
        // The work items cover the plan exactly, with no gaps or overlaps.
        let mut expected = 0u64;
        for item in &work {
            assert_eq!(item.first_ordinal, expected, "a gap in the ordinals");
            expected = item.end_ordinal;
        }
        assert_eq!(
            expected, plan.selected,
            "not every sampled record is covered"
        );

        // Every ordinal really lives in the block claiming it.
        let blocks = input.blocks();
        for item in &work {
            let block = blocks[item.block];
            for ordinal in item.first_ordinal..item.end_ordinal {
                let record = plan.index(ordinal);
                assert!(
                    record >= block.first_record && record < block.end_record(),
                    "record {record} is not in block {}",
                    item.block
                );
            }
        }
    }
}

/// Confidence and decision names round-trip through the report vocabulary.
#[test]
fn the_reported_vocabulary_is_closed() {
    for confidence in [Confidence::Low, Confidence::Medium, Confidence::High] {
        assert!(["low", "medium", "high"].contains(&confidence.name()));
    }
    for decision in [Decision::Confident, Decision::Mixed, Decision::Inconclusive] {
        assert!(["confident", "mixed", "inconclusive"].contains(&decision.name()));
        assert_eq!(decision.is_confident(), decision == Decision::Confident);
    }
}
