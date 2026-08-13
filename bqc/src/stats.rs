// Copyright (c) 2026 The Hiller Lab at the Senckenberg Gesellschaft für Naturforschung
// Distributed under the terms of the GNU General Public License, Version 3.0.

//! Operational statistics.
//!
//! Every worker accumulates a local [`Stats`] value; the ordered committer
//! merges them. Merging is pure addition, so the totals are independent of
//! thread count and chunk scheduling.

use std::collections::BTreeMap;

use crate::filter::FilterReason;
use crate::process::{Mate, PairResult};
use crate::read::Span;
use crate::trim::{TRIM_OPS, TrimOp, TrimOutcome};

/// Per-run counters.
#[derive(Debug, Clone, Default)]
pub struct Stats {
    /// Records (pairs, for paired files) read from the input.
    pub records_in: u64,
    /// Records written to the accepted output.
    pub records_out: u64,
    /// Records written to the failed output.
    pub records_rejected: u64,
    /// Records whose surviving mate was written to an orphan output, per mate.
    pub records_orphaned: [u64; 2],
    /// Bases written to the orphan outputs, per mate. Kept apart from
    /// `bases_out` so that the accepted-output record and base counts always
    /// describe the same set of records.
    pub bases_orphaned: [u64; 2],

    /// Bases read, per mate.
    pub bases_in: [u64; 2],
    /// Bases written to the accepted output, per mate.
    pub bases_out: [u64; 2],

    /// Reads whose 3' end was shortened by adapter removal, per mate.
    pub adapter_reads: [u64; 2],
    /// Bases removed by adapter removal, per mate.
    pub adapter_bases: [u64; 2],
    /// Hits per declared adapter, per mate.
    pub adapter_hits: [Vec<u64>; 2],
    /// Reads shortened by paired-overlap inference, per mate.
    pub overlap_reads: [u64; 2],
    /// Bases removed by paired-overlap inference, per mate.
    pub overlap_bases: [u64; 2],

    /// Reads shortened by each trimming operation.
    pub trim_reads: [u64; TRIM_OPS],
    /// Bases removed by each trimming operation.
    pub trim_bases: [u64; TRIM_OPS],

    /// Reads (mates) failing each individual reason.
    pub reason_counts: [u64; 7],
    /// Records failing each unique combination of reasons, keyed by bits.
    pub reason_combinations: BTreeMap<u32, u64>,
    /// Pairs where only R1 failed.
    pub r1_only_failed: u64,
    /// Pairs where only R2 failed.
    pub r2_only_failed: u64,
    pub both_failed: u64,

    /// Pairs the correction stage examined.
    pub correction_examined: u64,
    /// Pairs with an accepted overlap.
    pub correction_overlaps: u64,
    /// Pairs whose accepted overlap contained at least one mismatch.
    pub correction_mismatched_pairs: u64,
    /// Mismatches seen inside accepted overlaps.
    pub correction_mismatches: u64,
    /// Pairs where at least one base was corrected.
    pub corrected_pairs: u64,
    /// Pairs where both mates were corrected.
    pub corrected_both_mates: u64,
    /// Reads corrected, per mate.
    pub corrected_reads: [u64; 2],
    /// Bases corrected, per mate. R1 is corrected *from* R2 and vice versa, so
    /// `corrected_bases[0]` is also the R2-to-R1 direction count.
    pub corrected_bases: [u64; 2],
    /// Mismatches the quality rule declined to correct.
    pub correction_unresolved: u64,
    /// Mismatches skipped because the donor base was not canonical.
    pub correction_noncanonical: u64,
    /// Corrections per pair, for pairs with at least one correction.
    pub correction_histogram: BTreeMap<u32, u64>,
    /// Original-to-corrected substitutions, indexed `[original][corrected]` over
    /// `A, C, G, T, N`.
    pub correction_substitutions: [[u64; 5]; 5],
    /// Corrected pairs by where they were finally routed.
    pub corrected_by_disposition: BTreeMap<&'static str, u64>,

    /// Reads where both linked flanks matched, per mate.
    pub linked_both: [u64; 2],
    /// Reads where only the 5' flank matched, per mate.
    pub linked_five_only: [u64; 2],
    /// Reads where only the 3' flank matched, per mate.
    pub linked_three_only: [u64; 2],
    /// Reads whose flanks appeared in the wrong order, per mate.
    pub linked_out_of_order: [u64; 2],
    /// Reads whose insert was below the minimum, per mate.
    pub linked_short_insert: [u64; 2],
    /// Reads where neither flank matched, per mate.
    pub linked_neither: [u64; 2],
    /// Bases removed before the retained insert, per mate.
    pub linked_leading_bases: [u64; 2],
    /// Bases removed after the retained insert, per mate.
    pub linked_trailing_bases: [u64; 2],
    /// Matches per linked definition, per mate.
    pub linked_hits: [Vec<u64>; 2],

    /// Source records the segmentation stage examined.
    pub segment_sources: u64,
    /// Source records that produced at least one delimiter.
    pub segment_split_sources: u64,
    /// Delimiters accepted across all source records.
    pub segment_boundaries: u64,
    /// Candidates dropped for overlapping an accepted delimiter.
    pub segment_suppressed: u64,
    /// Fragments emitted with an adapter on both sides.
    pub segment_internal: u64,
    /// Fragments emitted that touch a read end.
    pub segment_terminal: u64,
    /// Fragments discarded for being empty.
    pub segment_empty: u64,
    /// Fragments discarded for being shorter than the minimum.
    pub segment_short: u64,
    /// Fragments discarded because terminals are discarded.
    pub segment_terminal_discarded: u64,
    /// Fragments discarded because the source read hit the safety limit.
    pub segment_over_limit: u64,
    /// Fragments emitted per source record.
    pub segment_histogram: BTreeMap<u32, u64>,
    /// Most fragments any one source record produced.
    pub segment_max_fragments: u64,
    /// Length summary over emitted fragments, before trimming: total, shortest,
    /// longest. Zero is the "no fragment yet" sentinel for the minimum, which is
    /// unambiguous because an empty fragment is never emitted.
    pub segment_length_total: u64,
    pub segment_length_min: u64,
    pub segment_length_max: u64,
    /// Records whose headers received a UMI tag.
    pub records_tagged: u64,
    /// Bases removed from the front of each mate by UMI extraction.
    pub umi_bases_removed: [u64; 2],
}

/// Adds per-key hit counts, growing `into` to cover keys only `from` has seen.
fn merge_counts(into: &mut Vec<u64>, from: &[u64]) {
    if into.len() < from.len() {
        into.resize(from.len(), 0);
    }
    for (total, count) in into.iter_mut().zip(from) {
        *total += count;
    }
}

impl Stats {
    /// Creates a statistics accumulator sized for the declared adapter counts.
    #[must_use]
    pub fn new(r1_adapters: usize, r2_adapters: usize) -> Self {
        Self {
            adapter_hits: [vec![0; r1_adapters], vec![0; r2_adapters]],
            ..Self::default()
        }
    }

    /// Sizes the per-linked-definition counters.
    #[must_use]
    pub fn with_linked(mut self, r1: usize, r2: usize) -> Self {
        self.linked_hits = [vec![0; r1], vec![0; r2]];
        self
    }

    /// Folds another accumulator's segmentation counters into this one.
    fn merge_segments(&mut self, other: &Self) {
        self.segment_sources += other.segment_sources;
        self.segment_split_sources += other.segment_split_sources;
        self.segment_boundaries += other.segment_boundaries;
        self.segment_suppressed += other.segment_suppressed;
        self.segment_internal += other.segment_internal;
        self.segment_terminal += other.segment_terminal;
        self.segment_empty += other.segment_empty;
        self.segment_short += other.segment_short;
        self.segment_terminal_discarded += other.segment_terminal_discarded;
        self.segment_over_limit += other.segment_over_limit;
        for (fragments, count) in &other.segment_histogram {
            *self.segment_histogram.entry(*fragments).or_insert(0) += count;
        }
        self.segment_max_fragments = self.segment_max_fragments.max(other.segment_max_fragments);
        self.segment_length_total += other.segment_length_total;
        if other.segment_length_min != 0
            && (self.segment_length_min == 0 || other.segment_length_min < self.segment_length_min)
        {
            self.segment_length_min = other.segment_length_min;
        }
        self.segment_length_max = self.segment_length_max.max(other.segment_length_max);
    }

    /// Records one segmented source record and the fragments it produced.
    ///
    /// The output counters are the same ones [`Stats::record`] maintains, because
    /// a fragment is an output record like any other; only the segmentation
    /// counters distinguish source records from fragments.
    pub fn record_segments(
        &mut self,
        original_length: usize,
        outcome: crate::segment::SegmentOutcome,
        delimiters: &[crate::adapter::AdapterHit],
        fragments: &[crate::process::FragmentResult],
    ) {
        self.records_in += 1;
        self.bases_in[0] += original_length as u64;
        self.segment_sources += 1;
        if outcome.boundaries > 0 {
            self.segment_split_sources += 1;
        }
        self.segment_boundaries += outcome.boundaries as u64;
        self.segment_suppressed += outcome.suppressed as u64;
        self.segment_internal += outcome.internal as u64;
        self.segment_terminal += outcome.terminal as u64;
        self.segment_empty += outcome.empty as u64;
        self.segment_short += outcome.too_short as u64;
        self.segment_terminal_discarded += outcome.terminal_discarded as u64;
        self.segment_over_limit += outcome.over_limit as u64;
        *self
            .segment_histogram
            .entry(outcome.emitted() as u32)
            .or_insert(0) += 1;
        self.segment_max_fragments = self.segment_max_fragments.max(outcome.emitted() as u64);

        // A delimiter is credited to the adapter that matched it, in the same
        // per-adapter table ordinary matching fills. Counted from the delimiters
        // rather than from the fragments they bound, because a delimiter at a read
        // end bounds a fragment that is discarded for being empty.
        for hit in delimiters {
            if let Some(count) = self.adapter_hits[0].get_mut(hit.adapter) {
                *count += 1;
            }
        }
        for fragment in fragments {
            // Fragment lengths are summarised as cut, before trimming: they
            // describe segmentation, not what a later stage did to the pieces.
            let length = fragment.fragment.span.len() as u64;
            self.segment_length_total += length;
            if self.segment_length_min == 0 || length < self.segment_length_min {
                self.segment_length_min = length;
            }
            self.segment_length_max = self.segment_length_max.max(length);
            self.count_trim(&fragment.trimming);
            self.count_reasons(fragment.reasons);
            if fragment.passed() {
                self.records_out += 1;
                self.bases_out[0] += fragment.retained.len() as u64;
            } else {
                self.records_rejected += 1;
                *self
                    .reason_combinations
                    .entry(fragment.reasons.bits())
                    .or_insert(0) += 1;
            }
        }
    }

    /// Records one processed record.
    pub fn record(&mut self, result: &PairResult) {
        self.records_in += 1;
        if result.umi_tagged {
            self.records_tagged += 1;
        }
        self.record_correction(result);
        for (mate, mate_result) in result.mates() {
            let index = mate.index();
            self.bases_in[index] += mate_result.original_length as u64;
            self.umi_bases_removed[index] += mate_result.umi_clip as u64;

            let adapter_removed = mate_result.adapter_bases_removed() as u64;
            if adapter_removed > 0 {
                self.adapter_reads[index] += 1;
                self.adapter_bases[index] += adapter_removed;
            }
            if let Some(hit) = mate_result.adapter_hit
                && let Some(count) = self.adapter_hits[index].get_mut(hit.adapter)
            {
                *count += 1;
            }
            if mate_result.overlap_removed > 0 {
                self.overlap_reads[index] += 1;
                self.overlap_bases[index] += mate_result.overlap_removed as u64;
            }

            self.count_trim(&mate_result.trimming);

            self.record_linked(mate, mate_result);

            self.count_reasons(mate_result.reasons);
        }

        match result.disposition {
            crate::process::PairDisposition::Accepted => {
                self.records_out += 1;
                for (mate, mate_result) in result.mates() {
                    self.bases_out[mate.index()] += mate_result.final_length() as u64;
                }
            }
            disposition => {
                let combined = result
                    .mates()
                    .fold(FilterReason::empty(), |acc, (_, m)| acc | m.reasons);
                *self.reason_combinations.entry(combined.bits()).or_insert(0) += 1;

                let r1_failed = !result.r1.passed();
                let r2_failed = result.r2.as_ref().is_some_and(|r| !r.passed());
                match (r1_failed, r2_failed) {
                    (true, true) => self.both_failed += 1,
                    (true, false) => self.r1_only_failed += 1,
                    (false, true) => self.r2_only_failed += 1,
                    (false, false) => {}
                }

                match disposition {
                    crate::process::PairDisposition::Rejected => self.records_rejected += 1,
                    crate::process::PairDisposition::OrphanR1 => {
                        self.records_orphaned[0] += 1;
                        self.bases_orphaned[0] += result.r1.final_length() as u64;
                    }
                    crate::process::PairDisposition::OrphanR2 => {
                        self.records_orphaned[1] += 1;
                        if let Some(r2) = result.r2.as_ref() {
                            self.bases_orphaned[1] += r2.final_length() as u64;
                        }
                    }
                    crate::process::PairDisposition::Accepted => unreachable!(),
                }
            }
        }
    }

    /// Accumulates the correction facts for one record.
    ///
    /// Counts describe corrections *applied to the read*, before trimming: a
    /// corrected base that a later stage trims away is still counted here, and
    /// the disposition breakdown says where those pairs ended up.
    fn record_correction(&mut self, result: &PairResult) {
        if result.r2.is_none() {
            return; // correction is paired-only
        }
        let summary = result.correction;
        self.correction_examined += 1;
        if result.overlap.is_some() {
            self.correction_overlaps += 1;
        }
        self.correction_mismatches += u64::from(summary.mismatches);
        if summary.mismatches > 0 {
            self.correction_mismatched_pairs += 1;
        }
        self.correction_unresolved += u64::from(summary.unresolved);
        self.correction_noncanonical += u64::from(summary.skipped_noncanonical);
        if summary.corrected() == 0 {
            return;
        }
        self.corrected_pairs += 1;
        if summary.both_mates() {
            self.corrected_both_mates += 1;
        }
        for mate in [Mate::R1, Mate::R2] {
            let corrected = summary.corrected_in(mate);
            if corrected > 0 {
                self.corrected_reads[mate.index()] += 1;
                self.corrected_bases[mate.index()] += u64::from(corrected);
            }
        }
        *self
            .correction_histogram
            .entry(summary.corrected())
            .or_insert(0) += 1;
        *self
            .corrected_by_disposition
            .entry(result.disposition.name())
            .or_insert(0) += 1;
    }

    /// Accumulates the linked-segmentation facts for one mate.
    fn record_linked(&mut self, mate: Mate, result: &crate::process::MateResult) {
        use crate::linked::{LinkedOutcome, Rejected};
        let Some(outcome) = result.linked else {
            return;
        };
        let index = mate.index();
        match outcome {
            LinkedOutcome::Matched(found) => {
                if found.both() {
                    self.linked_both[index] += 1;
                } else if found.five_prime.is_some() {
                    self.linked_five_only[index] += 1;
                } else {
                    self.linked_three_only[index] += 1;
                }
                let searched = Span::full(result.original_length);
                self.linked_leading_bases[index] += found.leading_removed(searched) as u64;
                self.linked_trailing_bases[index] += found.trailing_removed(searched) as u64;
                if let Some(count) = self.linked_hits[index].get_mut(found.definition) {
                    *count += 1;
                }
            }
            LinkedOutcome::Unmatched(reason) => match reason {
                Rejected::Neither => self.linked_neither[index] += 1,
                Rejected::FiveOnly => self.linked_five_only[index] += 1,
                Rejected::ThreeOnly => self.linked_three_only[index] += 1,
                Rejected::OutOfOrder => self.linked_out_of_order[index] += 1,
                Rejected::InsertTooShort => self.linked_short_insert[index] += 1,
            },
        }
    }

    /// Records one applied substitution, before the base is mutated.
    pub fn record_substitution(&mut self, original: u8, corrected: u8) {
        self.correction_substitutions[base_slot(original)][base_slot(corrected)] += 1;
    }

    /// Counts one mate's trim operations into the per-op totals.
    fn count_trim(&mut self, trimming: &TrimOutcome) {
        for op in TrimOp::ALL {
            let removed = u64::from(trimming.removed[op as usize]);
            if removed > 0 {
                self.trim_reads[op as usize] += 1;
                self.trim_bases[op as usize] += removed;
            }
        }
    }

    /// Counts one mate's filter reasons into the per-reason totals.
    fn count_reasons(&mut self, reasons: FilterReason) {
        for (bit, (reason, _)) in crate::filter::REASONS.iter().enumerate() {
            if reasons.contains(*reason) {
                self.reason_counts[bit] += 1;
            }
        }
    }

    /// Adds `other` into `self`.
    pub fn merge(&mut self, other: &Self) {
        self.merge_segments(other);
        self.records_in += other.records_in;
        self.records_out += other.records_out;
        self.records_rejected += other.records_rejected;
        for mate in 0..2 {
            self.records_orphaned[mate] += other.records_orphaned[mate];
            self.bases_in[mate] += other.bases_in[mate];
            self.bases_out[mate] += other.bases_out[mate];
            self.bases_orphaned[mate] += other.bases_orphaned[mate];
            self.adapter_reads[mate] += other.adapter_reads[mate];
            self.adapter_bases[mate] += other.adapter_bases[mate];
            self.overlap_reads[mate] += other.overlap_reads[mate];
            self.overlap_bases[mate] += other.overlap_bases[mate];
            merge_counts(&mut self.adapter_hits[mate], &other.adapter_hits[mate]);
        }
        for op in 0..TRIM_OPS {
            self.trim_reads[op] += other.trim_reads[op];
            self.trim_bases[op] += other.trim_bases[op];
        }
        for bit in 0..self.reason_counts.len() {
            self.reason_counts[bit] += other.reason_counts[bit];
        }
        for (bits, count) in &other.reason_combinations {
            *self.reason_combinations.entry(*bits).or_insert(0) += count;
        }
        self.correction_examined += other.correction_examined;
        self.correction_overlaps += other.correction_overlaps;
        self.correction_mismatched_pairs += other.correction_mismatched_pairs;
        self.correction_mismatches += other.correction_mismatches;
        self.corrected_pairs += other.corrected_pairs;
        self.corrected_both_mates += other.corrected_both_mates;
        self.correction_unresolved += other.correction_unresolved;
        self.correction_noncanonical += other.correction_noncanonical;
        for mate in 0..2 {
            self.corrected_reads[mate] += other.corrected_reads[mate];
            self.corrected_bases[mate] += other.corrected_bases[mate];
        }
        for (corrections, count) in &other.correction_histogram {
            *self.correction_histogram.entry(*corrections).or_insert(0) += count;
        }
        for (original, row) in other.correction_substitutions.iter().enumerate() {
            for (corrected, count) in row.iter().enumerate() {
                self.correction_substitutions[original][corrected] += count;
            }
        }
        for (disposition, count) in &other.corrected_by_disposition {
            *self
                .corrected_by_disposition
                .entry(disposition)
                .or_insert(0) += count;
        }
        for mate in 0..2 {
            self.linked_both[mate] += other.linked_both[mate];
            self.linked_five_only[mate] += other.linked_five_only[mate];
            self.linked_three_only[mate] += other.linked_three_only[mate];
            self.linked_out_of_order[mate] += other.linked_out_of_order[mate];
            self.linked_short_insert[mate] += other.linked_short_insert[mate];
            self.linked_neither[mate] += other.linked_neither[mate];
            self.linked_leading_bases[mate] += other.linked_leading_bases[mate];
            self.linked_trailing_bases[mate] += other.linked_trailing_bases[mate];
            merge_counts(&mut self.linked_hits[mate], &other.linked_hits[mate]);
        }
        self.r1_only_failed += other.r1_only_failed;
        self.r2_only_failed += other.r2_only_failed;
        self.both_failed += other.both_failed;
        self.records_tagged += other.records_tagged;
        self.umi_bases_removed[0] += other.umi_bases_removed[0];
        self.umi_bases_removed[1] += other.umi_bases_removed[1];
    }

    /// Total bases read across both mates.
    #[must_use]
    pub fn total_bases_in(&self) -> u64 {
        self.bases_in[0] + self.bases_in[1]
    }

    /// Hits recorded for one declared adapter.
    #[must_use]
    pub fn adapter_hit_count(&self, mate: Mate, index: usize) -> u64 {
        self.adapter_hits[mate.index()]
            .get(index)
            .copied()
            .unwrap_or(0)
    }
}

/// Substitution-matrix slot of a base: `A, C, G, T` then anything else.
fn base_slot(base: u8) -> usize {
    match base {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => 4,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::adapter::{Adapter, AdapterParams, AdapterStage};
    use crate::filter::FilterStage;
    use crate::process::Workflow;
    use crate::read::ReadView;
    use crate::trim::{MateTrim, TrimStage};

    const ADAPTER: &[u8] = b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";

    fn workflow() -> Workflow {
        let adapter = AdapterStage::new(
            vec![Adapter::new("r1", ADAPTER).unwrap()],
            vec![Adapter::new("r2", ADAPTER).unwrap()],
            AdapterParams::default(),
            None,
        )
        .unwrap();
        let trim = TrimStage {
            r1: MateTrim {
                tail: 2,
                ..MateTrim::default()
            },
            r2: MateTrim {
                tail: 2,
                ..MateTrim::default()
            },
        };
        let filter = FilterStage {
            min_length: Some(16),
            qualified_quality: 15,
            ..FilterStage::default()
        };
        Workflow::new(Some(adapter), None, Some(trim), Some(filter), None, None).unwrap()
    }

    #[test]
    fn statistics_track_every_stage() {
        let workflow = workflow();
        let mut stats = Stats::new(1, 1);

        let with_adapter = [b"ACGTACGTACGTACGTACGT".as_slice(), ADAPTER].concat();
        let r1 = ReadView::unchecked(&with_adapter, None);
        let clean = b"ACGTACGTACGTACGTACGTACGT";
        let r2 = ReadView::unchecked(clean, None);
        stats.record(&workflow.process(0, r1, Some(r2)).unwrap());

        assert_eq!(stats.records_in, 1);
        assert_eq!(
            stats.bases_in,
            [with_adapter.len() as u64, clean.len() as u64]
        );
        assert_eq!(stats.adapter_reads, [1, 0]);
        assert_eq!(stats.adapter_bases, [ADAPTER.len() as u64, 0]);
        assert_eq!(stats.adapter_hits[0], vec![1]);
        assert_eq!(stats.adapter_hits[1], vec![0]);
        assert_eq!(
            stats.trim_reads[TrimOp::Fixed as usize],
            2,
            "both mates cut"
        );
        assert_eq!(stats.trim_bases[TrimOp::Fixed as usize], 4);
        assert_eq!(stats.records_out, 1);
        assert_eq!(stats.bases_out, [18, 22]);
    }

    #[test]
    fn rejected_records_are_attributed_per_mate() {
        let filter = FilterStage {
            min_length: Some(8),
            qualified_quality: 15,
            ..FilterStage::default()
        };
        let workflow = Workflow::new(None, None, None, Some(filter), None, None).unwrap();
        let mut stats = Stats::new(0, 0);
        let long = ReadView::unchecked(b"ACGTACGTACGT", None);
        let short = ReadView::unchecked(b"ACGT", None);

        stats.record(&workflow.process(0, long, Some(short)).unwrap());
        stats.record(&workflow.process(1, short, Some(long)).unwrap());
        stats.record(&workflow.process(2, short, Some(short)).unwrap());
        stats.record(&workflow.process(3, long, Some(long)).unwrap());

        assert_eq!(stats.records_in, 4);
        assert_eq!(stats.records_out, 1);
        assert_eq!(stats.records_rejected, 3);
        assert_eq!(stats.r1_only_failed, 1);
        assert_eq!(stats.r2_only_failed, 1);
        assert_eq!(stats.both_failed, 1);
        assert_eq!(stats.reason_counts[0], 4, "four short mates in total");
        assert_eq!(
            stats
                .reason_combinations
                .get(&FilterReason::TOO_SHORT.bits()),
            Some(&3)
        );
    }

    #[test]
    fn merging_is_pure_addition() {
        let workflow = workflow();
        let sequence = [b"ACGTACGTACGTACGTACGT".as_slice(), ADAPTER].concat();
        let read = ReadView::unchecked(&sequence, None);

        let mut a = Stats::new(1, 1);
        a.record(&workflow.process(0, read, Some(read)).unwrap());
        let mut b = Stats::new(1, 1);
        b.record(&workflow.process(1, read, Some(read)).unwrap());

        let mut merged = Stats::new(1, 1);
        merged.merge(&a);
        merged.merge(&b);

        let mut sequential = Stats::new(1, 1);
        sequential.record(&workflow.process(0, read, Some(read)).unwrap());
        sequential.record(&workflow.process(1, read, Some(read)).unwrap());

        assert_eq!(merged.records_in, sequential.records_in);
        assert_eq!(merged.adapter_hits, sequential.adapter_hits);
        assert_eq!(merged.trim_bases, sequential.trim_bases);
        assert_eq!(merged.bases_out, sequential.bases_out);
        assert_eq!(merged.reason_combinations, sequential.reason_combinations);
    }
}
