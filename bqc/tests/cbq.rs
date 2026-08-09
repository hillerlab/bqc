// Copyright (c) 2026 The Hiller Lab at the Senckenberg Gesellschaft für Naturforschung
// Distributed under the terms of the GNU General Public License, Version 3.0.

//! End-to-end tests over real CBQ files.

mod common;

use std::path::{Path, PathBuf};

use bqc::io::Schema;
use common::{
    Planted, Record, Sequences, all_schemas, bqc, bqc_ok, complement, count_blocks, describe,
    overlapping_pairs, read_cbq, write_cbq,
};
use tempfile::TempDir;

const ADAPTER_R1: &str = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
const ADAPTER_R2: &str = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT";
const BLOCK: usize = 1 << 20;

struct Fixture {
    dir: TempDir,
}

impl Fixture {
    fn new() -> Self {
        Self {
            dir: tempfile::tempdir().expect("tempdir"),
        }
    }

    fn path(&self, name: &str) -> PathBuf {
        self.dir.path().join(name)
    }

    /// Writes an input file and returns its path.
    fn input(&self, schema: Schema, records: &[Record], block_size: usize) -> PathBuf {
        let path = self.path("input.cbq");
        write_cbq(&path, schema, records, block_size);
        path
    }
}

/// Parses a JSON report.
fn json(path: &Path) -> serde_json::Value {
    serde_json::from_str(&text(path)).expect("parse report")
}

fn text(path: &Path) -> String {
    std::fs::read_to_string(path).expect("read text output")
}

fn simple_records(schema: Schema, count: usize, length: usize) -> Vec<Record> {
    let mut sequences = Sequences::new(0xB0_1234);
    (0..count)
        .map(|index| {
            let r1 = sequences.sequence(length);
            let r2 = sequences.sequence(length);
            let q1 = sequences.quality(length);
            let q2 = sequences.quality(length);
            Record::new(schema, index, &r1, &r2).with_quality(&q1, &q2)
        })
        .collect()
}

/// Records carrying an adapter on both mates after a 20 base insert.
fn adapter_records(schema: Schema, count: usize) -> Vec<Record> {
    let mut sequences = Sequences::new(0x5EED);
    (0..count)
        .map(|index| {
            let insert = sequences.sequence(20);
            let r1 = [insert.as_slice(), ADAPTER_R1.as_bytes()].concat();
            let r2 = [insert.as_slice(), ADAPTER_R2.as_bytes()].concat();
            let q1 = vec![b'I'; r1.len()];
            let q2 = vec![b'I'; r2.len()];
            Record::new(schema, index, &r1, &r2).with_quality(&q1, &q2)
        })
        .collect()
}

// ---------------------------------------------------------------- schema matrix

#[test]
fn every_schema_round_trips_through_a_trim() {
    for schema in all_schemas() {
        let fixture = Fixture::new();
        let records = simple_records(schema, 12, 40);
        let input = fixture.input(schema, &records, BLOCK);
        let output = fixture.path("out.cbq");

        bqc_ok(&[
            "trim",
            input.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
            "--front",
            "3",
            "--tail",
            "5",
        ]);

        let (out_schema, out_records) = read_cbq(&output);
        assert_eq!(
            out_schema,
            schema,
            "schema changed for {}",
            describe(schema)
        );
        assert_eq!(out_records.len(), records.len(), "{}", describe(schema));

        for (index, (before, after)) in records.iter().zip(&out_records).enumerate() {
            let expected = &before.s_seq[3..before.s_seq.len() - 5];
            assert_eq!(after.s_seq, expected, "R1 {index} of {}", describe(schema));
            if schema.quality {
                let quality = before.s_qual.as_ref().unwrap();
                assert_eq!(
                    after.s_qual.as_ref().unwrap(),
                    &quality[3..quality.len() - 5],
                    "R1 quality {index} of {}",
                    describe(schema)
                );
            }
            if schema.headers {
                assert_eq!(after.s_header, before.s_header, "header {index}");
            }
            if schema.flags {
                assert_eq!(after.flag, before.flag, "flag {index}");
            }
            if schema.paired {
                let sequence = before.x_seq.as_ref().unwrap();
                assert_eq!(
                    after.x_seq.as_ref().unwrap(),
                    &sequence[3..sequence.len() - 5],
                    "R2 {index}"
                );
                if schema.quality {
                    let quality = before.x_qual.as_ref().unwrap();
                    assert_eq!(
                        after.x_qual.as_ref().unwrap(),
                        &quality[3..quality.len() - 5],
                        "R2 quality {index}"
                    );
                }
                if schema.headers {
                    assert_eq!(after.x_header, before.x_header, "R2 header {index}");
                }
            }
        }
    }
}

#[test]
fn flags_survive_a_filtering_run() {
    // Guards the metadata-propagation trap called out in the design: it is easy
    // to copy sequence, header and quality while dropping the flag column.
    let schema = Schema {
        paired: true,
        quality: true,
        headers: false,
        flags: true,
    };
    let fixture = Fixture::new();
    let mut records = simple_records(schema, 20, 40);
    records[3].s_seq.truncate(5);
    records[3].s_qual.as_mut().unwrap().truncate(5);
    let input = fixture.input(schema, &records, BLOCK);
    let output = fixture.path("out.cbq");

    bqc_ok(&[
        "filter",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--min-length",
        "10",
    ]);

    let (out_schema, out_records) = read_cbq(&output);
    assert!(out_schema.flags);
    let expected: Vec<Option<u64>> = records
        .iter()
        .enumerate()
        .filter(|(index, _)| *index != 3)
        .map(|(_, record)| record.flag)
        .collect();
    let actual: Vec<Option<u64>> = out_records.iter().map(|record| record.flag).collect();
    assert_eq!(actual, expected, "flags must be preserved and stay aligned");
}

#[test]
fn empty_and_single_record_inputs_are_handled() {
    let schema = Schema {
        paired: false,
        quality: true,
        headers: true,
        flags: false,
    };
    for count in [0usize, 1] {
        let fixture = Fixture::new();
        let records = simple_records(schema, count, 30);
        let input = fixture.input(schema, &records, BLOCK);
        let output = fixture.path("out.cbq");
        bqc_ok(&[
            "trim",
            input.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
            "--tail",
            "1",
        ]);
        let (out_schema, out_records) = read_cbq(&output);
        assert_eq!(out_schema, schema);
        assert_eq!(out_records.len(), count);
    }
}

#[test]
fn variable_length_reads_and_ambiguous_bases_are_preserved() {
    let schema = Schema {
        paired: true,
        quality: true,
        headers: true,
        flags: false,
    };
    let fixture = Fixture::new();
    let mut sequences = Sequences::new(0xFACE);
    let records: Vec<Record> = (0..30)
        .map(|index| {
            let length = 10 + index * 3;
            let r1 = sequences.sequence(length);
            let r2 = sequences.sequence(length + 1);
            let q1 = sequences.quality(length);
            let q2 = sequences.quality(length + 1);
            Record::new(schema, index, &r1, &r2).with_quality(&q1, &q2)
        })
        .collect();
    let input = fixture.input(schema, &records, 4096);
    let output = fixture.path("out.cbq");

    bqc_ok(&[
        "trim",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--front",
        "2",
    ]);

    let (_, out_records) = read_cbq(&output);
    assert_eq!(out_records.len(), records.len());
    for (before, after) in records.iter().zip(&out_records) {
        assert_eq!(after.s_seq, before.s_seq[2..]);
        assert_eq!(
            after.x_seq.as_ref().unwrap(),
            &before.x_seq.as_ref().unwrap()[2..]
        );
    }
}

// ------------------------------------------------------------------- ordering

#[test]
fn record_order_and_output_bytes_are_independent_of_thread_count() {
    let schema = Schema {
        paired: true,
        quality: true,
        headers: true,
        flags: true,
    };
    let fixture = Fixture::new();
    // Small blocks so the file has many independently scheduled chunks.
    let records = simple_records(schema, 600, 60);
    let input = fixture.input(schema, &records, 4096);
    assert!(count_blocks(&input) > 4, "fixture must span several blocks");

    let mut outputs = Vec::new();
    for threads in ["1", "2", "3", "8"] {
        let output = fixture.path(&format!("out{threads}.cbq"));
        bqc_ok(&[
            "workflow",
            input.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
            "-T",
            threads,
            "--quality-tail",
            "20",
            "--min-length",
            "20",
        ]);
        outputs.push(output);
    }

    let (_, reference) = read_cbq(&outputs[0]);
    assert!(!reference.is_empty());
    let reference_bytes = std::fs::read(&outputs[0]).unwrap();
    for output in &outputs[1..] {
        let (_, records) = read_cbq(output);
        assert_eq!(
            records, reference,
            "decoded output differs between thread counts"
        );
        assert_eq!(
            std::fs::read(output).unwrap(),
            reference_bytes,
            "output bytes differ between thread counts"
        );
    }

    // Order is the input order, not merely the same multiset.
    let expected: Vec<Vec<u8>> = records
        .iter()
        .map(|record| record.s_header.clone().unwrap())
        .filter(|header| {
            reference
                .iter()
                .any(|out| out.s_header.as_ref() == Some(header))
        })
        .collect();
    let actual: Vec<Vec<u8>> = reference
        .iter()
        .map(|record| record.s_header.clone().unwrap())
        .collect();
    assert_eq!(actual, expected, "records must stay in input order");
}

#[test]
fn output_block_layout_follows_the_input() {
    let schema = Schema {
        paired: false,
        quality: true,
        headers: false,
        flags: false,
    };
    let fixture = Fixture::new();
    let records = simple_records(schema, 400, 50);
    let input = fixture.input(schema, &records, 4096);
    let output = fixture.path("out.cbq");
    bqc_ok(&[
        "trim",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--tail",
        "1",
        "-T",
        "4",
    ]);
    assert_eq!(
        count_blocks(&output),
        count_blocks(&input),
        "one input block should produce one output block"
    );
}

// ---------------------------------------------------------- command equivalence

fn assert_same_file(a: &Path, b: &Path) {
    assert_eq!(read_cbq(a).1, read_cbq(b).1, "decoded records differ");
    assert_eq!(
        std::fs::read(a).unwrap(),
        std::fs::read(b).unwrap(),
        "bytes differ"
    );
}

#[test]
fn standalone_commands_match_the_equivalent_workflow() {
    let schema = Schema {
        paired: true,
        quality: true,
        headers: true,
        flags: false,
    };
    let fixture = Fixture::new();
    let records = adapter_records(schema, 40);
    let input = fixture.input(schema, &records, 8192);

    let cases: Vec<(&str, Vec<&str>)> = vec![
        (
            "adapter",
            vec![
                "--adapter-r1",
                ADAPTER_R1,
                "--adapter-r2",
                ADAPTER_R2,
                "--min-overlap",
                "8",
            ],
        ),
        (
            "trim",
            vec!["--front", "2", "--quality-tail", "20", "--trim-terminal-n"],
        ),
        (
            "filter",
            vec![
                "--min-length",
                "20",
                "--max-n",
                "2",
                "--min-mean-quality",
                "20",
            ],
        ),
    ];

    for (step, args) in cases {
        let standalone = fixture.path(&format!("{step}-standalone.cbq"));
        let mut invocation = vec![
            step,
            input.to_str().unwrap(),
            "-o",
            standalone.to_str().unwrap(),
        ];
        invocation.extend(args.iter().copied());
        bqc_ok(&invocation);

        let fused = fixture.path(&format!("{step}-workflow.cbq"));
        let mut invocation = vec![
            "workflow",
            input.to_str().unwrap(),
            "-o",
            fused.to_str().unwrap(),
            "--steps",
            step,
        ];
        invocation.extend(args.iter().copied());
        bqc_ok(&invocation);

        assert_same_file(&standalone, &fused);
    }
}

#[test]
fn fusing_stages_matches_running_them_in_sequence() {
    let schema = Schema {
        paired: true,
        quality: true,
        headers: true,
        flags: true,
    };
    let fixture = Fixture::new();
    let records = adapter_records(schema, 60);
    let input = fixture.input(schema, &records, 8192);

    let after_adapter = fixture.path("step1.cbq");
    bqc_ok(&[
        "adapter",
        input.to_str().unwrap(),
        "-o",
        after_adapter.to_str().unwrap(),
        "--adapter-r1",
        ADAPTER_R1,
        "--adapter-r2",
        ADAPTER_R2,
    ]);
    let after_trim = fixture.path("step2.cbq");
    bqc_ok(&[
        "trim",
        after_adapter.to_str().unwrap(),
        "-o",
        after_trim.to_str().unwrap(),
        "--tail",
        "2",
    ]);
    let staged = fixture.path("step3.cbq");
    bqc_ok(&[
        "filter",
        after_trim.to_str().unwrap(),
        "-o",
        staged.to_str().unwrap(),
        "--min-length",
        "10",
    ]);

    let fused = fixture.path("fused.cbq");
    bqc_ok(&[
        "workflow",
        input.to_str().unwrap(),
        "-o",
        fused.to_str().unwrap(),
        "--adapter-r1",
        ADAPTER_R1,
        "--adapter-r2",
        ADAPTER_R2,
        "--tail",
        "2",
        "--min-length",
        "10",
        "-T",
        "3",
    ]);

    assert_eq!(read_cbq(&staged).1, read_cbq(&fused).1);
}

// -------------------------------------------------------------------- adapters

#[test]
fn adapters_are_trimmed_and_counted() {
    let schema = Schema {
        paired: true,
        quality: true,
        headers: true,
        flags: false,
    };
    let fixture = Fixture::new();
    let records = adapter_records(schema, 25);
    let input = fixture.input(schema, &records, BLOCK);
    let output = fixture.path("out.cbq");
    let report = fixture.path("report.json");

    bqc_ok(&[
        "adapter",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--adapter-r1",
        ADAPTER_R1,
        "--adapter-r2",
        ADAPTER_R2,
        "--report",
        report.to_str().unwrap(),
    ]);

    let (_, out_records) = read_cbq(&output);
    for (index, record) in out_records.iter().enumerate() {
        assert_eq!(record.s_seq, records[index].s_seq[..20], "R1 {index}");
        assert_eq!(
            record.x_seq.as_ref().unwrap(),
            &records[index].x_seq.as_ref().unwrap()[..20],
            "R2 {index}"
        );
    }

    let report: serde_json::Value = serde_json::from_str(&text(&report)).unwrap();
    assert_eq!(report["adapter"]["r1_reads_trimmed"], 25);
    assert_eq!(report["adapter"]["r2_reads_trimmed"], 25);
    assert_eq!(report["adapter"]["r1_bases_removed"], 25 * 33);
    assert_eq!(report["adapter"]["per_adapter"][0]["name"], "r1");
    assert_eq!(report["adapter"]["per_adapter"][0]["hits"], 25);
    assert_eq!(report["configuration"]["stage_order"][0], "adapter");
}

#[test]
fn adapters_can_be_supplied_as_a_fasta_file() {
    let schema = Schema {
        paired: false,
        quality: true,
        headers: true,
        flags: false,
    };
    let fixture = Fixture::new();
    let records = adapter_records(schema, 10);
    let input = fixture.input(schema, &records, BLOCK);
    let fasta = fixture.path("adapters.fa");
    std::fs::write(
        &fasta,
        format!(">truseq\n{ADAPTER_R1}\n>nextera\nCTGTCTCTTATACACATCTCCGAGCCCACGAGAC\n"),
    )
    .unwrap();
    let output = fixture.path("out.cbq");
    let report = fixture.path("report.json");

    bqc_ok(&[
        "adapter",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--adapter-fasta",
        fasta.to_str().unwrap(),
        "--report",
        report.to_str().unwrap(),
    ]);

    let (_, out_records) = read_cbq(&output);
    for (index, record) in out_records.iter().enumerate() {
        assert_eq!(record.s_seq, records[index].s_seq[..20]);
    }
    let report: serde_json::Value = serde_json::from_str(&text(&report)).unwrap();
    let per_adapter = report["adapter"]["per_adapter"].as_array().unwrap();
    assert_eq!(per_adapter[0]["name"], "truseq");
    assert_eq!(per_adapter[0]["hits"], 10);
    assert_eq!(per_adapter[1]["name"], "nextera");
    assert_eq!(per_adapter[1]["hits"], 0);
}

// --------------------------------------------------------------------- filters

#[test]
fn strict_pairing_routes_rejected_pairs_and_records_every_reason() {
    let schema = Schema {
        paired: true,
        quality: true,
        headers: true,
        flags: false,
    };
    let fixture = Fixture::new();
    let mut records = simple_records(schema, 6, 40);
    // Record 1: R1 too short. Record 2: R2 too short. Record 3: both.
    records[1].s_seq.truncate(5);
    records[1].s_qual.as_mut().unwrap().truncate(5);
    records[2].x_seq.as_mut().unwrap().truncate(5);
    records[2].x_qual.as_mut().unwrap().truncate(5);
    records[3].s_seq.truncate(4);
    records[3].s_qual.as_mut().unwrap().truncate(4);
    records[3].x_seq.as_mut().unwrap().truncate(4);
    records[3].x_qual.as_mut().unwrap().truncate(4);

    let input = fixture.input(schema, &records, BLOCK);
    let output = fixture.path("passed.cbq");
    let failed = fixture.path("failed.cbq");
    let reasons = fixture.path("reasons.tsv");
    let report = fixture.path("report.json");

    bqc_ok(&[
        "filter",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--failed",
        failed.to_str().unwrap(),
        "--failed-reasons",
        reasons.to_str().unwrap(),
        "--report",
        report.to_str().unwrap(),
        "--min-length",
        "10",
    ]);

    let (_, accepted) = read_cbq(&output);
    let (failed_schema, rejected) = read_cbq(&failed);
    assert_eq!(accepted.len(), 3);
    assert_eq!(rejected.len(), 3);
    assert_eq!(
        failed_schema, schema,
        "the failed output keeps the input schema"
    );
    // Rejected pairs are written unmodified by default.
    assert_eq!(rejected[0], records[1]);
    assert_eq!(rejected[1], records[2]);
    assert_eq!(rejected[2], records[3]);

    let sidecar = text(&reasons);
    let lines: Vec<&str> = sidecar.lines().collect();
    assert!(lines[0].starts_with("record_index\tmate\tstatus\treasons"));
    assert_eq!(lines.len(), 1 + 6, "two rows per rejected pair");
    assert_eq!(lines[1], "1\tR1\tFAIL\tTOO_SHORT\t5\t5\t5\t.\t.");
    assert_eq!(lines[2], "1\tR2\tPASS\tPASS\t40\t40\t40\t.\t.");
    assert_eq!(lines[3], "2\tR1\tPASS\tPASS\t40\t40\t40\t.\t.");
    assert_eq!(lines[4], "2\tR2\tFAIL\tTOO_SHORT\t5\t5\t5\t.\t.");
    assert_eq!(lines[5], "3\tR1\tFAIL\tTOO_SHORT\t4\t4\t4\t.\t.");
    assert_eq!(lines[6], "3\tR2\tFAIL\tTOO_SHORT\t4\t4\t4\t.\t.");

    let report: serde_json::Value = serde_json::from_str(&text(&report)).unwrap();
    assert_eq!(report["filter"]["accepted_records"], 3);
    assert_eq!(report["filter"]["rejected_records"], 3);
    assert_eq!(report["filter"]["r1_only_failed"], 1);
    assert_eq!(report["filter"]["r2_only_failed"], 1);
    assert_eq!(report["filter"]["both_failed"], 1);
    assert_eq!(report["filter"]["per_reason"][0]["reasons"], "TOO_SHORT");
    assert_eq!(report["filter"]["per_reason"][0]["records"], 4);
    assert_eq!(
        report["filter"]["per_combination"][0]["reasons"],
        "TOO_SHORT"
    );
    assert_eq!(report["filter"]["per_combination"][0]["records"], 3);
    // Accounting reconciles with the input.
    assert_eq!(report["counts"]["records_in"], 6);
    assert_eq!(report["counts"]["records_out"], 3);
    assert_eq!(report["counts"]["records_rejected"], 3);
}

#[test]
fn failed_mode_processed_writes_the_transformed_record() {
    let schema = Schema {
        paired: false,
        quality: true,
        headers: true,
        flags: false,
    };
    let fixture = Fixture::new();
    let records = simple_records(schema, 4, 40);
    let input = fixture.input(schema, &records, BLOCK);

    for (mode, expected_length) in [("original", 40), ("processed", 20)] {
        let output = fixture.path(&format!("passed-{mode}.cbq"));
        let failed = fixture.path(&format!("failed-{mode}.cbq"));
        bqc_ok(&[
            "workflow",
            input.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
            "--failed",
            failed.to_str().unwrap(),
            "--failed-mode",
            mode,
            "--tail",
            "20",
            "--min-length",
            "30",
        ]);
        let (_, rejected) = read_cbq(&failed);
        assert_eq!(rejected.len(), 4, "every record fails --min-length 30");
        assert_eq!(rejected[0].s_seq.len(), expected_length, "mode {mode}");
        assert_eq!(rejected[0].s_qual.as_ref().unwrap().len(), expected_length);
        assert_eq!(read_cbq(&output).1.len(), 0);
    }
}

#[test]
fn the_sidecar_reports_the_adapter_that_caused_the_trim() {
    let schema = Schema {
        paired: false,
        quality: true,
        headers: true,
        flags: false,
    };
    let fixture = Fixture::new();
    let records = adapter_records(schema, 3);
    let input = fixture.input(schema, &records, BLOCK);
    let output = fixture.path("out.cbq");
    let reasons = fixture.path("reasons.tsv");
    let fasta = fixture.path("adapters.fa");
    std::fs::write(&fasta, format!(">truseq\n{ADAPTER_R1}\n")).unwrap();

    // The insert is 20 bases, so every read fails --min-length 30 after the
    // adapter is removed: each rejected row must name the adapter and its
    // coordinate in the original read.
    bqc_ok(&[
        "workflow",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--adapter-fasta",
        fasta.to_str().unwrap(),
        "--failed-reasons",
        reasons.to_str().unwrap(),
        "--min-length",
        "30",
    ]);

    let sidecar = text(&reasons);
    let rows: Vec<&str> = sidecar.lines().skip(1).collect();
    assert_eq!(rows.len(), 3);
    for (index, row) in rows.iter().enumerate() {
        assert_eq!(
            *row,
            format!("{index}\tR1\tFAIL\tTOO_SHORT\t53\t20\t20\ttruseq\t20")
        );
    }
    assert_eq!(read_cbq(&output).1.len(), 0);
}

#[test]
fn multiple_failure_reasons_are_all_reported() {
    let schema = Schema {
        paired: false,
        quality: true,
        headers: true,
        flags: false,
    };
    let fixture = Fixture::new();
    let records = vec![Record::new(schema, 0, b"NNNNNNNN", b"").with_quality(b"!!!!!!!!", b"")];
    let input = fixture.input(schema, &records, BLOCK);
    let output = fixture.path("out.cbq");
    let reasons = fixture.path("reasons.tsv");

    bqc_ok(&[
        "filter",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--failed-reasons",
        reasons.to_str().unwrap(),
        "--min-length",
        "10",
        "--max-n",
        "2",
        "--min-mean-quality",
        "20",
        "--min-complexity",
        "0.3",
        "--qualified-quality",
        "15",
        "--max-unqualified-fraction",
        "0.4",
    ]);

    let sidecar = text(&reasons);
    let row = sidecar.lines().nth(1).unwrap();
    assert!(row.contains("TOO_SHORT/TOO_MANY_N/TOO_MANY_LOW_QUAL/LOW_MEAN_QUAL/LOW_COMPLEXITY"));
    assert_eq!(read_cbq(&output).1.len(), 0);
}

// ------------------------------------------------------------------------ span

#[test]
fn span_selects_original_record_indices() {
    let schema = Schema {
        paired: false,
        quality: true,
        headers: true,
        flags: false,
    };
    let fixture = Fixture::new();
    let records = simple_records(schema, 50, 30);
    let input = fixture.input(schema, &records, 1024);
    assert!(count_blocks(&input) > 1);
    let output = fixture.path("out.cbq");

    bqc_ok(&[
        "trim",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--span",
        "10..20",
        "--tail",
        "1",
        "-T",
        "4",
    ]);

    let (_, out_records) = read_cbq(&output);
    assert_eq!(out_records.len(), 10);
    for (offset, record) in out_records.iter().enumerate() {
        assert_eq!(record.s_header, records[10 + offset].s_header);
    }
}

#[test]
fn span_bounds_are_validated_and_clamped() {
    let schema = Schema {
        paired: false,
        quality: true,
        headers: true,
        flags: false,
    };
    let fixture = Fixture::new();
    let records = simple_records(schema, 10, 30);
    let input = fixture.input(schema, &records, BLOCK);
    let output = fixture.path("out.cbq");

    // An end past the last record is clamped.
    bqc_ok(&[
        "trim",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--span",
        "5..",
        "--tail",
        "1",
    ]);
    assert_eq!(read_cbq(&output).1.len(), 5);

    // A start past the last record is an error.
    let other = fixture.path("other.cbq");
    let error = bqc(&[
        "trim",
        input.to_str().unwrap(),
        "-o",
        other.to_str().unwrap(),
        "--span",
        "50..60",
        "--tail",
        "1",
    ])
    .unwrap_err();
    assert!(
        format!("{error}").contains("past the last record"),
        "{error}"
    );
    assert!(!other.exists(), "no partial output is left behind");
}

// ------------------------------------------------------------- error behaviour

#[test]
fn quality_operations_require_a_quality_column() {
    let schema = Schema {
        paired: false,
        quality: false,
        headers: true,
        flags: false,
    };
    let fixture = Fixture::new();
    let records = simple_records(schema, 5, 30);
    let input = fixture.input(schema, &records, BLOCK);
    let output = fixture.path("out.cbq");

    let error = bqc(&[
        "trim",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--quality-tail",
        "20",
    ])
    .unwrap_err();
    assert!(format!("{error}").contains("no quality column"), "{error}");
    assert!(!output.exists());

    // Sequence-only operations still work on the same file.
    bqc_ok(&[
        "trim",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--front",
        "2",
    ]);
    assert_eq!(read_cbq(&output).1.len(), 5);
}

#[test]
fn existing_outputs_are_protected_unless_forced() {
    let schema = Schema {
        paired: false,
        quality: true,
        headers: true,
        flags: false,
    };
    let fixture = Fixture::new();
    let records = simple_records(schema, 5, 30);
    let input = fixture.input(schema, &records, BLOCK);
    let output = fixture.path("out.cbq");
    std::fs::write(&output, b"existing").unwrap();

    let error = bqc(&[
        "trim",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--front",
        "1",
    ])
    .unwrap_err();
    assert!(
        format!("{error}").contains("refusing to overwrite"),
        "{error}"
    );
    assert_eq!(std::fs::read(&output).unwrap(), b"existing");

    bqc_ok(&[
        "trim",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--front",
        "1",
        "--force",
    ]);
    assert_eq!(read_cbq(&output).1.len(), 5);
}

#[test]
fn the_input_cannot_be_used_as_an_output() {
    let schema = Schema {
        paired: false,
        quality: true,
        headers: true,
        flags: false,
    };
    let fixture = Fixture::new();
    let records = simple_records(schema, 5, 30);
    let input = fixture.input(schema, &records, BLOCK);

    let error = bqc(&[
        "trim",
        input.to_str().unwrap(),
        "-o",
        input.to_str().unwrap(),
        "--front",
        "1",
        "--force",
    ])
    .unwrap_err();
    assert!(format!("{error}").contains("same file"), "{error}");
    assert_eq!(read_cbq(&input).1.len(), 5, "the input is untouched");
}

#[test]
fn non_cbq_input_is_rejected_with_guidance() {
    let fixture = Fixture::new();
    let input = fixture.path("reads.fastq");
    std::fs::write(&input, b"@read\nACGT\n+\nIIII\n").unwrap();
    let output = fixture.path("out.cbq");
    let error = bqc(&[
        "trim",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--front",
        "1",
    ])
    .unwrap_err();
    let message = format!("{error}");
    assert!(message.contains("not a readable CBQ file"), "{message}");
    assert!(message.contains("bqtools encode"), "{message}");
}

#[test]
fn no_operation_configured_is_an_error() {
    let schema = Schema {
        paired: false,
        quality: true,
        headers: true,
        flags: false,
    };
    let fixture = Fixture::new();
    let records = simple_records(schema, 3, 30);
    let input = fixture.input(schema, &records, BLOCK);
    let output = fixture.path("out.cbq");

    for command in ["adapter", "trim", "filter", "workflow"] {
        let error = bqc(&[
            command,
            input.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
        ])
        .unwrap_err();
        let message = format!("{error}");
        assert!(
            message.contains("requires") || message.contains("no operation configured"),
            "{command}: {message}"
        );
        assert!(!output.exists(), "{command} left an output behind");
    }
}

// ---------------------------------------------------------------- config files

#[test]
fn a_toml_configuration_drives_a_workflow_and_the_cli_overrides_it() {
    let schema = Schema {
        paired: true,
        quality: true,
        headers: true,
        flags: false,
    };
    let fixture = Fixture::new();
    let records = adapter_records(schema, 30);
    let input = fixture.input(schema, &records, BLOCK);
    let config = fixture.path("illumina.toml");
    let output = fixture.path("clean.cbq");
    let report = fixture.path("report.json");
    std::fs::write(
        &config,
        format!(
            r#"
threads = 2

[adapter]
r1 = "{ADAPTER_R1}"
r2 = "{ADAPTER_R2}"
min_overlap = 8
max_error_rate = 0.10

[trim.quality_tail]
minimum_phred = 20
window = 4

[trim]
terminal_n = true

[filter]
min_length = 15
max_n = 5
"#
        ),
    )
    .unwrap();

    bqc_ok(&[
        "workflow",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--config",
        config.to_str().unwrap(),
        "--report",
        report.to_str().unwrap(),
        "--min-length",
        "25",
    ]);

    let report: serde_json::Value = serde_json::from_str(&text(&report)).unwrap();
    assert_eq!(
        report["configuration"]["threads"], 2,
        "threads come from the file"
    );
    assert_eq!(
        report["configuration"]["workflow"]["filter"]["min_length"], 25,
        "the command line overrides the file"
    );
    assert_eq!(
        report["configuration"]["workflow"]["trim"]["r1"]["terminal_n"],
        true
    );
    assert_eq!(
        report["configuration"]["stage_order"],
        serde_json::json!(["adapter", "trim", "filter"])
    );
    // Every read is 20 bases after adapter removal, so --min-length 25 rejects
    // everything: the override really took effect.
    assert_eq!(read_cbq(&output).1.len(), 0);
}

#[test]
fn reports_can_be_written_as_tsv() {
    let schema = Schema {
        paired: false,
        quality: true,
        headers: true,
        flags: false,
    };
    let fixture = Fixture::new();
    let records = simple_records(schema, 5, 30);
    let input = fixture.input(schema, &records, BLOCK);
    let output = fixture.path("out.cbq");
    let report = fixture.path("report.tsv");

    bqc_ok(&[
        "trim",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--front",
        "1",
        "--report",
        report.to_str().unwrap(),
        "--report-format",
        "tsv",
    ]);

    let rendered = text(&report);
    assert!(rendered.starts_with("key\tvalue\n"));
    assert!(rendered.contains("counts.records_in\t5\n"), "{rendered}");
    assert!(rendered.contains("tool\tbqc\n"));
    assert!(
        rendered.contains("configuration.workflow.trim.r1.front\t1\n"),
        "{rendered}"
    );
}

#[test]
fn reports_describe_the_resolved_configuration() {
    let schema = Schema {
        paired: false,
        quality: true,
        headers: true,
        flags: false,
    };
    let fixture = Fixture::new();
    let records = simple_records(schema, 5, 30);
    let input = fixture.input(schema, &records, BLOCK);
    let output = fixture.path("out.cbq");
    let report = fixture.path("report.json");

    bqc_ok(&[
        "trim",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--poly-g",
        "--report",
        report.to_str().unwrap(),
    ]);

    let report: serde_json::Value = serde_json::from_str(&text(&report)).unwrap();
    assert_eq!(report["tool"], "bqc");
    assert_eq!(report["version"], env!("CARGO_PKG_VERSION"));
    assert_eq!(report["binseq_version"], "0.9.4");
    assert_eq!(report["command"], "trim");
    assert_eq!(report["input"]["records"], 5);
    assert_eq!(report["input"]["schema"]["quality"], true);
    assert_eq!(report["input"]["blocks"], 1);
    // Defaults that the user never typed are still recorded.
    assert_eq!(
        report["configuration"]["workflow"]["trim"]["r1"]["poly_g"]["min_length"],
        10
    );
    assert_eq!(report["configuration"]["workflow"]["pair_policy"], "strict");
    assert_eq!(report["configuration"]["span"], serde_json::Value::Null);
    assert_eq!(report["outputs"]["accepted"], output.to_str().unwrap());
    assert!(report["performance"]["elapsed_seconds"].as_f64().unwrap() >= 0.0);
    assert_eq!(report["trim"][0]["operation"], "poly_g");
}

// ---------------------------------------------------- adversarial review fixes

#[test]
fn report_conflicts_abort_the_run_before_any_output_is_written() {
    let schema = Schema {
        paired: false,
        quality: true,
        headers: true,
        flags: false,
    };
    let fixture = Fixture::new();
    let records = simple_records(schema, 8, 40);
    let input = fixture.input(schema, &records, BLOCK);
    let output = fixture.path("out.cbq");
    let report = fixture.path("report.json");
    std::fs::write(&report, b"PRECIOUS").unwrap();

    // An existing report path without --force must fail up front.
    let error = bqc(&[
        "trim",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--front",
        "2",
        "--report",
        report.to_str().unwrap(),
    ])
    .unwrap_err();
    assert!(
        format!("{error}").contains("refusing to overwrite"),
        "{error}"
    );
    assert!(
        !output.exists(),
        "the accepted output must not exist after an aborted run"
    );
    assert_eq!(
        std::fs::read(&report).unwrap(),
        b"PRECIOUS",
        "the existing report is untouched"
    );

    // The report pointing at the input must fail up front, even with --force.
    let error = bqc(&[
        "trim",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--front",
        "2",
        "--report",
        input.to_str().unwrap(),
        "--force",
    ])
    .unwrap_err();
    assert!(format!("{error}").contains("same file"), "{error}");
    assert!(!output.exists());
    assert_eq!(read_cbq(&input).1.len(), 8, "the input is untouched");
}

#[test]
fn max_length_zero_is_rejected_for_either_mate() {
    let schema = Schema {
        paired: true,
        quality: true,
        headers: true,
        flags: false,
    };
    let fixture = Fixture::new();
    let records = simple_records(schema, 4, 40);
    let input = fixture.input(schema, &records, BLOCK);

    for (r1, r2) in [("0", "5"), ("5", "0"), ("0", "0")] {
        let output = fixture.path(&format!("out-{r1}-{r2}.cbq"));
        let error = bqc(&[
            "trim",
            input.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
            "--max-length-r1",
            r1,
            "--max-length-r2",
            r2,
        ])
        .unwrap_err();
        assert!(
            format!("{error}").contains("--max-length must be at least 1"),
            "r1={r1} r2={r2}: {error}"
        );
        assert!(!output.exists());
    }
}

#[test]
fn truncated_input_is_rejected_instead_of_silently_accepted() {
    let schema = Schema {
        paired: true,
        quality: true,
        headers: true,
        flags: true,
    };
    let fixture = Fixture::new();
    // Small blocks, so the file has several of them.
    let records = simple_records(schema, 200, 40);
    let input = fixture.input(schema, &records, 512);
    assert!(count_blocks(&input) > 2, "fixture must have several blocks");
    let full = std::fs::read(&input).unwrap();

    // Locate the index trailer: exactly two CBQINDEX magics exist, the index
    // header and the index footer.
    let positions: Vec<usize> = full
        .windows(8)
        .enumerate()
        .filter_map(|(i, w)| (w == b"CBQINDEX").then_some(i))
        .collect();
    assert_eq!(positions.len(), 2, "index header and footer magics");
    let index_start = positions[0];

    // Dropping the index and the trailing block must not be accepted.
    let truncated = fixture.path("cut-at-boundary.cbq");
    std::fs::write(&truncated, &full[..index_start]).unwrap();
    let output = fixture.path("out.cbq");
    let error = bqc(&[
        "trim",
        truncated.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--front",
        "1",
    ])
    .unwrap_err();
    assert!(format!("{error}").contains("truncated"), "{error}");
    assert!(!output.exists());

    // Cutting through the middle of a block must not be accepted either.
    let mid_block = fixture.path("cut-mid-block.cbq");
    std::fs::write(&mid_block, &full[..index_start / 2]).unwrap();
    let error = bqc(&[
        "trim",
        mid_block.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--front",
        "1",
    ])
    .unwrap_err();
    assert!(format!("{error}").contains("corrupt"), "{error}");

    // The untouched file processes normally.
    bqc_ok(&[
        "trim",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--front",
        "1",
    ]);
    assert_eq!(read_cbq(&output).1.len(), 200);
}

// ---------------------------------------------------------------- orphan outputs

/// Pairs with per-mate pass/fail quality patterns.
///
/// The pattern cycles every four records: both pass, R2 fails, R1 fails, both
/// fail. Passing mates are Q40, failing mates Q0, and every sequence is
/// distinct so records can be told apart.
fn orphan_records(schema: Schema, count: usize) -> Vec<Record> {
    (0..count)
        .map(|i| {
            let r1: Vec<u8> = (0..30).map(|j| b"ACGT"[(i + j) % 4]).collect();
            let r2: Vec<u8> = (0..30).map(|j| b"TGCA"[(i + j) % 4]).collect();
            // cycle: 0 = both pass, 1 = R2 fails, 2 = R1 fails, 3 = both fail
            let (p1, p2) = match i % 4 {
                0 => (true, true),
                1 => (true, false),
                2 => (false, true),
                _ => (false, false),
            };
            let q1 = vec![if p1 { b'I' } else { b'!' }; 30];
            let q2 = vec![if p2 { b'I' } else { b'!' }; 30];
            Record::new(schema, i, &r1, &r2).with_quality(&q1, &q2)
        })
        .collect()
}

#[test]
// One run writes five files; checking them together is the point of the test.
#[allow(clippy::too_many_lines)]
fn orphan_outputs_route_surviving_mates_to_single_end_files() {
    let schema = Schema {
        paired: true,
        quality: true,
        headers: true,
        flags: true,
    };
    let fixture = Fixture::new();
    let records = orphan_records(schema, 6);
    let input = fixture.input(schema, &records, BLOCK);
    let accepted = fixture.path("clean.cbq");
    let rejected = fixture.path("rejected.cbq");
    let reasons = fixture.path("reasons.tsv");
    let prefix = fixture.path("surviving");
    let report = fixture.path("report.json");

    bqc_ok(&[
        "filter",
        input.to_str().unwrap(),
        "-o",
        accepted.to_str().unwrap(),
        "--failed",
        rejected.to_str().unwrap(),
        "--failed-reasons",
        reasons.to_str().unwrap(),
        "--pair-policy",
        "orphan",
        "--orphan-prefix",
        prefix.to_str().unwrap(),
        "--min-mean-quality",
        "20",
        "--report",
        report.to_str().unwrap(),
        "-T",
        "2",
    ]);

    // The accepted output holds the two pairs where both mates pass.
    let (accepted_schema, accepted_records) = read_cbq(&accepted);
    assert_eq!(accepted_schema, schema, "the accepted output stays paired");
    let accepted_headers: Vec<_> = accepted_records
        .iter()
        .map(|r| r.s_header.clone().unwrap())
        .collect();
    assert_eq!(
        accepted_headers,
        [b"read_0/1".to_vec(), b"read_4/1".to_vec()]
    );

    // The R1 orphan file holds the R1 mates of records 1 and 5, single-end.
    let orphan_r1_path = fixture.path("surviving.R1.cbq");
    let (orphan_r1_schema, orphan_r1) = read_cbq(&orphan_r1_path);
    assert_eq!(
        orphan_r1_schema,
        schema.unpaired(),
        "orphan files are single-end, schema otherwise preserved"
    );
    assert_eq!(orphan_r1.len(), 2);
    assert_eq!(orphan_r1[0].s_seq, records[1].s_seq);
    assert_eq!(orphan_r1[1].s_seq, records[5].s_seq);
    assert_eq!(orphan_r1[0].s_header, records[1].s_header);
    assert_eq!(orphan_r1[0].flag, records[1].flag, "flags are preserved");
    assert_eq!(orphan_r1[0].x_seq, None, "no mate column in an orphan file");

    // The R2 orphan file holds the R2 mate of record 2 — not its R1.
    let orphan_r2_path = fixture.path("surviving.R2.cbq");
    let (_, orphan_r2) = read_cbq(&orphan_r2_path);
    assert_eq!(orphan_r2.len(), 1);
    assert_eq!(orphan_r2[0].s_seq, records[2].x_seq.clone().unwrap());
    assert_eq!(orphan_r2[0].s_header, records[2].x_header);
    assert_eq!(orphan_r2[0].flag, records[2].flag);

    // The rejected output holds record 3, both mates, in original form.
    let (_, rejected_records) = read_cbq(&rejected);
    assert_eq!(rejected_records.len(), 1);
    assert_eq!(rejected_records[0], records[3]);

    // Every input record is accounted for exactly once.
    assert_eq!(
        accepted_records.len() + orphan_r1.len() + orphan_r2.len() + rejected_records.len(),
        records.len()
    );

    // The sidecar explains every non-accepted record, per mate.
    let sidecar = text(&reasons);
    let rows: Vec<&str> = sidecar.lines().skip(1).collect();
    assert_eq!(rows.len(), 8, "four broken records, two mates each");
    assert_eq!(
        rows[0], "1\tR1\tPASS\tPASS\t30\t30\t30\t.\t.",
        "the surviving mate is recorded as PASS"
    );
    assert_eq!(rows[1], "1\tR2\tFAIL\tLOW_MEAN_QUAL\t30\t30\t30\t.\t.");
    assert!(
        rows.iter()
            .any(|row| row.starts_with("3\tR1\tFAIL\tLOW_MEAN_QUAL"))
    );

    // The report carries the orphan counts and paths.
    let report: serde_json::Value = serde_json::from_str(&text(&report)).unwrap();
    assert_eq!(report["filter"]["accepted_records"], 2);
    assert_eq!(report["filter"]["rejected_records"], 1);
    assert_eq!(report["filter"]["orphan_r1_records"], 2);
    assert_eq!(report["filter"]["orphan_r2_records"], 1);
    assert_eq!(report["filter"]["r2_only_failed"], 2);
    assert_eq!(report["filter"]["r1_only_failed"], 1);
    assert_eq!(report["filter"]["both_failed"], 1);
    assert_eq!(
        report["outputs"]["orphan_r1"],
        orphan_r1_path.to_str().unwrap()
    );
    assert_eq!(
        report["outputs"]["orphan_r2"],
        orphan_r2_path.to_str().unwrap()
    );
    assert_eq!(report["configuration"]["workflow"]["pair_policy"], "orphan");
    assert_eq!(report["configuration"]["output"]["pair_policy"], "orphan");

    // The accepted record and base counts must describe the same records:
    // orphan bases belong to the orphan files and are counted separately.
    let counts = &report["counts"];
    assert_eq!(counts["records_out"], 2);
    let accepted_bases = accepted_records
        .iter()
        .map(|record| record.s_seq.len() + record.x_seq.as_ref().map_or(0, Vec::len))
        .sum::<usize>();
    assert_eq!(
        counts["r1_bases_out"].as_u64().unwrap() + counts["r2_bases_out"].as_u64().unwrap(),
        accepted_bases as u64,
        "accepted bases must match the accepted output only"
    );
    let orphan_bases: usize = orphan_r1
        .iter()
        .chain(&orphan_r2)
        .map(|record| record.s_seq.len())
        .sum();
    assert_eq!(
        counts["r1_bases_orphaned"].as_u64().unwrap()
            + counts["r2_bases_orphaned"].as_u64().unwrap(),
        orphan_bases as u64,
        "orphan bases must match the orphan outputs"
    );
}

#[test]
fn orphan_outputs_are_identical_at_any_thread_count() {
    let schema = Schema {
        paired: true,
        quality: true,
        headers: true,
        flags: true,
    };
    let fixture = Fixture::new();
    let records = orphan_records(schema, 400);
    // Small blocks so the parallel engine has several chunks to reorder.
    let input = fixture.input(schema, &records, 512);
    assert!(count_blocks(&input) > 2);

    let mut runs = Vec::new();
    for threads in [1, 4] {
        let accepted = fixture.path(&format!("clean-{threads}.cbq"));
        let prefix = fixture.path(&format!("surviving-{threads}"));
        bqc_ok(&[
            "workflow",
            input.to_str().unwrap(),
            "-o",
            accepted.to_str().unwrap(),
            "--pair-policy",
            "orphan",
            "--orphan-prefix",
            prefix.to_str().unwrap(),
            "--min-mean-quality",
            "20",
            "-T",
            &threads.to_string(),
        ]);
        runs.push((
            accepted,
            fixture.path(&format!("surviving-{threads}.R1.cbq")),
            fixture.path(&format!("surviving-{threads}.R2.cbq")),
        ));
    }
    let [(a1, r1_1, r2_1), (a2, r1_2, r2_2)] = runs.as_slice() else {
        unreachable!()
    };
    for (one, four) in [(a1, a2), (r1_1, r1_2), (r2_1, r2_2)] {
        assert_same_file(one, four);
    }
    // 100 of each class: both pass, R2 fails, R1 fails, both fail.
    assert_eq!(read_cbq(a1).1.len(), 100);
    assert_eq!(read_cbq(r1_1).1.len(), 100);
    assert_eq!(read_cbq(r2_1).1.len(), 100);
}

#[test]
fn the_orphan_policy_is_validated() {
    let paired = Schema {
        paired: true,
        quality: true,
        headers: true,
        flags: false,
    };
    let single_end = Schema {
        paired: false,
        ..paired
    };
    let fixture = Fixture::new();
    let paired_input = fixture.input(paired, &orphan_records(paired, 4), BLOCK);
    let se_input = fixture.path("se.cbq");
    write_cbq(
        &se_input,
        single_end,
        &simple_records(single_end, 4, 30),
        BLOCK,
    );
    let output = fixture.path("out.cbq");
    let prefix = fixture.path("surviving");

    // Single-end input cannot produce orphans.
    let error = bqc(&[
        "filter",
        se_input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--pair-policy",
        "orphan",
        "--orphan-prefix",
        prefix.to_str().unwrap(),
        "--min-length",
        "10",
    ])
    .unwrap_err();
    assert!(format!("{error}").contains("paired input"), "{error}");

    // The policy needs a prefix.
    let error = bqc(&[
        "filter",
        paired_input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--pair-policy",
        "orphan",
        "--min-length",
        "10",
    ])
    .unwrap_err();
    assert!(format!("{error}").contains("--orphan-prefix"), "{error}");

    // The prefix needs the policy.
    let error = bqc(&[
        "filter",
        paired_input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--orphan-prefix",
        prefix.to_str().unwrap(),
        "--min-length",
        "10",
    ])
    .unwrap_err();
    assert!(
        format!("{error}").contains("--pair-policy orphan"),
        "{error}"
    );

    // Nothing can fail without a filter stage.
    let error = bqc(&[
        "trim",
        paired_input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--pair-policy",
        "orphan",
        "--orphan-prefix",
        prefix.to_str().unwrap(),
        "--front",
        "1",
    ])
    .unwrap_err();
    assert!(format!("{error}").contains("requires filtering"), "{error}");

    // An orphan path colliding with another destination is rejected.
    let error = bqc(&[
        "filter",
        paired_input.to_str().unwrap(),
        "-o",
        fixture.path("surviving.R1.cbq").to_str().unwrap(),
        "--pair-policy",
        "orphan",
        "--orphan-prefix",
        prefix.to_str().unwrap(),
        "--min-length",
        "10",
    ])
    .unwrap_err();
    assert!(
        format!("{error}").contains("more than one output"),
        "{error}"
    );
    assert!(!output.exists());
}

// ---------------------------------------------------------------- paired overlap

fn revcomp(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => b'N',
        })
        .collect()
}

/// Pairs with a 60-base insert and 30 bases of adapter read-through per mate.
fn read_through_records(schema: Schema, count: usize) -> Vec<Record> {
    let mut sequences = Sequences::new(0x000E_41A9);
    (0..count)
        .map(|index| {
            let mut insert = sequences.sequence(60);
            for base in &mut insert {
                if *base == b'N' {
                    *base = b'A';
                }
            }
            let tail_r1 = sequences.sequence(30);
            let tail_r2 = sequences.sequence(30);
            let r1 = [insert.as_slice(), &tail_r1].concat();
            let r2 = [revcomp(&insert).as_slice(), &tail_r2].concat();
            let q1 = vec![b'I'; r1.len()];
            let q2 = vec![b'I'; r2.len()];
            Record::new(schema, index, &r1, &r2).with_quality(&q1, &q2)
        })
        .collect()
}

#[test]
fn paired_overlap_trims_read_through_without_adapter_sequences() {
    let schema = Schema {
        paired: true,
        quality: true,
        headers: true,
        flags: true,
    };
    let fixture = Fixture::new();
    let records = read_through_records(schema, 20);
    let input = fixture.input(schema, &records, BLOCK);
    let output = fixture.path("out.cbq");
    let report = fixture.path("report.json");

    bqc_ok(&[
        "adapter",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--paired-overlap",
        "--report",
        report.to_str().unwrap(),
    ]);

    let (out_schema, out_records) = read_cbq(&output);
    assert_eq!(out_schema, schema);
    assert_eq!(out_records.len(), 20);
    for (original, trimmed) in records.iter().zip(&out_records) {
        assert_eq!(trimmed.s_seq, &original.s_seq[..60], "R1 cut at the insert");
        assert_eq!(
            trimmed.x_seq.as_deref(),
            Some(&original.x_seq.as_ref().unwrap()[..60]),
            "R2 cut at the insert"
        );
        assert_eq!(trimmed.s_header, original.s_header);
        assert_eq!(trimmed.flag, original.flag);
    }

    let report: serde_json::Value = serde_json::from_str(&text(&report)).unwrap();
    assert_eq!(report["adapter"]["r1_overlap_reads_trimmed"], 20);
    assert_eq!(report["adapter"]["r2_overlap_reads_trimmed"], 20);
    assert_eq!(report["adapter"]["r1_overlap_bases_removed"], 600);
    assert_eq!(
        report["configuration"]["workflow"]["adapter"]["paired_overlap"]["min_overlap"],
        30
    );
}

#[test]
fn indel_aware_matching_trims_adapters_with_insertions() {
    let schema = Schema {
        paired: false,
        quality: true,
        headers: true,
        flags: false,
    };
    let fixture = Fixture::new();
    // Three records: adapter with an insertion, adapter with a deletion, and
    // a clean read, each after a 20-base insert.
    let mut sequences = Sequences::new(0x1DE1);
    let mut with_insertion = ADAPTER_R1.as_bytes().to_vec();
    with_insertion.insert(10, b'C');
    let mut with_deletion = ADAPTER_R1.as_bytes().to_vec();
    with_deletion.remove(10);
    let inserts: Vec<Vec<u8>> = (0..3).map(|_| sequences.sequence(20)).collect();
    let tails: Vec<&[u8]> = vec![&with_insertion, &with_deletion, b"CCCCCCCCCCCC"];
    let records: Vec<Record> = (0..3)
        .map(|i| {
            let seq = [inserts[i].as_slice(), tails[i]].concat();
            let qual = vec![b'I'; seq.len()];
            Record::new(schema, i, &seq, b"").with_quality(&qual, b"")
        })
        .collect();
    let input = fixture.input(schema, &records, BLOCK);

    // Substitutions alone cannot trim the indel-containing adapters.
    let plain = fixture.path("plain.cbq");
    bqc_ok(&[
        "adapter",
        input.to_str().unwrap(),
        "-o",
        plain.to_str().unwrap(),
        "--adapter-r1",
        ADAPTER_R1,
    ]);
    let plain_records = read_cbq(&plain).1;
    assert_eq!(plain_records[0].s_seq.len(), 20 + with_insertion.len());
    assert_eq!(plain_records[1].s_seq.len(), 20 + with_deletion.len());

    // With --allow-indels both are cut at the insert boundary.
    let output = fixture.path("indel.cbq");
    bqc_ok(&[
        "adapter",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--adapter-r1",
        ADAPTER_R1,
        "--allow-indels",
        "-T",
        "2",
    ]);
    let out_records = read_cbq(&output).1;
    for (i, record) in out_records.iter().take(2).enumerate() {
        assert_eq!(
            record.s_seq.len(),
            20,
            "record {i} must be cut at the insert boundary"
        );
        assert_eq!(record.s_seq, inserts[i]);
        assert_eq!(record.s_header, records[i].s_header);
    }
    // The clean read is untouched.
    assert_eq!(out_records[2].s_seq.len(), 32);
}

#[test]
fn paired_overlap_is_validated() {
    let single_end = Schema {
        paired: false,
        quality: true,
        headers: true,
        flags: false,
    };
    let paired = Schema {
        paired: true,
        ..single_end
    };
    let fixture = Fixture::new();
    let se_input = fixture.input(single_end, &simple_records(single_end, 4, 40), BLOCK);
    let pe_input = fixture.path("paired.cbq");
    write_cbq(&pe_input, paired, &simple_records(paired, 4, 40), BLOCK);
    let output = fixture.path("out.cbq");

    let error = bqc(&[
        "adapter",
        se_input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--paired-overlap",
    ])
    .unwrap_err();
    assert!(format!("{error}").contains("paired input"), "{error}");

    let error = bqc(&[
        "adapter",
        pe_input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--paired-overlap-min-overlap",
        "20",
    ])
    .unwrap_err();
    assert!(
        format!("{error}").contains("requires --paired-overlap"),
        "{error}"
    );
    assert!(!output.exists());
}

#[test]
fn failed_detection_falls_back_to_paired_overlap() {
    // Detection finds nothing on clean data, but paired-overlap inference needs
    // no adapter sequence: asking for both still trims, and asking for
    // auto-detection alone passes the file through untrimmed.
    let schema = Schema {
        paired: true,
        quality: true,
        headers: true,
        flags: false,
    };
    let fixture = Fixture::new();
    let records = simple_records(schema, 40, 50);
    let input = fixture.input(schema, &records, BLOCK);
    let output = fixture.path("out.cbq");
    let report = fixture.path("report.json");

    bqc_ok(&[
        "adapter",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--auto-detect",
        "--paired-overlap",
        "--report",
        report.to_str().unwrap(),
    ]);
    assert_eq!(read_cbq(&output).1.len(), 40);
    let report: serde_json::Value = serde_json::from_str(&text(&report)).unwrap();
    assert_eq!(
        report["adapter"]["detection"]["r1"]["decision"],
        "inconclusive"
    );
    assert!(report["adapter"]["detection"]["r1"]["recommended_sequence"].is_null());
    assert!(
        !report["configuration"]["workflow"]["adapter"]["paired_overlap"].is_null(),
        "overlap inference stays configured"
    );

    // Without a fallback the file passes through untrimmed.
    let other = fixture.path("other.cbq");
    let other_report = fixture.path("other-report.json");
    bqc_ok(&[
        "adapter",
        input.to_str().unwrap(),
        "-o",
        other.to_str().unwrap(),
        "--auto-detect",
        "--report",
        other_report.to_str().unwrap(),
    ]);
    assert_eq!(read_cbq(&other).1.len(), 40);
    let report: serde_json::Value = serde_json::from_str(&text(&other_report)).unwrap();
    assert_eq!(
        report["adapter"]["detection"]["r1"]["decision"],
        "inconclusive"
    );
    assert!(report["adapter"]["detection"]["r1"]["recommended_sequence"].is_null());
    assert_eq!(report["adapter"]["r1_reads_trimmed"], 0);
}

#[test]
fn detection_samples_only_the_requested_span() {
    // Records outside --span must not contribute evidence: the contaminated
    // half of the file is excluded, so detection sees clean reads only.
    let schema = Schema {
        paired: false,
        quality: true,
        headers: true,
        flags: false,
    };
    let fixture = Fixture::new();
    let mut sequences = Sequences::new(0x5A11_0F5E);
    let records: Vec<Record> = (0..800)
        .map(|index| {
            let sequence = if index < 400 {
                // The first half is heavily contaminated.
                let insert = sequences.sequence(20);
                [insert.as_slice(), ADAPTER_R1.as_bytes()].concat()
            } else {
                sequences.sequence(53)
            };
            let quality = vec![b'I'; sequence.len()];
            Record::new(schema, index, &sequence, b"").with_quality(&quality, b"")
        })
        .collect();
    // Small blocks so the span really excludes whole blocks.
    let input = fixture.input(schema, &records, 1024);
    let output = fixture.path("out.cbq");
    let report = fixture.path("report.json");

    // Sampling the contaminated half detects the adapter.
    bqc_ok(&[
        "adapter",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--span",
        "0..400",
        "--auto-detect",
        "--report",
        report.to_str().unwrap(),
    ]);
    let parsed: serde_json::Value = serde_json::from_str(&text(&report)).unwrap();
    let detection = &parsed["adapter"]["detection"];
    assert_eq!(
        detection["r1"]["recommended_name"], "illumina-truseq",
        "{detection}"
    );
    assert_eq!(detection["r1"]["sampled_reads"], 400);

    // Sampling the clean half finds nothing and passes the file through.
    let clean = fixture.path("clean.cbq");
    let clean_report = fixture.path("clean-report.json");
    bqc_ok(&[
        "adapter",
        input.to_str().unwrap(),
        "-o",
        clean.to_str().unwrap(),
        "--span",
        "400..800",
        "--auto-detect",
        "--report",
        clean_report.to_str().unwrap(),
    ]);
    assert_eq!(read_cbq(&clean).1.len(), 400);
    let report: serde_json::Value = serde_json::from_str(&text(&clean_report)).unwrap();
    assert_eq!(
        report["adapter"]["detection"]["r1"]["decision"],
        "inconclusive"
    );
    assert_eq!(report["adapter"]["r1_reads_trimmed"], 0);
}

// ---------------------------------------------------------------- auto-detection

/// Pairs where `fraction` of records carry the given adapter prefixes after a
/// random insert, on both mates.
fn contaminated_records(
    schema: Schema,
    count: usize,
    r1_adapter: &[u8],
    r2_adapter: &[u8],
) -> Vec<Record> {
    let mut sequences = Sequences::new(0x00DE_7EC7);
    (0..count)
        .map(|index| {
            let insert = sequences.sequence(20 + index % 30);
            let mut r1 = sequences.sequence(40 + index % 30);
            let mut r2 = sequences.sequence(40 + index % 30);
            if index % 5 == 0 {
                // 20% contamination: random insert, then the adapter prefix.
                r1 = [insert.as_slice(), &r1_adapter[..20]].concat();
                r2 = [insert.as_slice(), &r2_adapter[..20]].concat();
            }
            let q1 = vec![b'I'; r1.len()];
            let q2 = vec![b'I'; r2.len()];
            Record::new(schema, index, &r1, &r2).with_quality(&q1, &q2)
        })
        .collect()
}

#[test]
fn auto_detection_finds_a_known_adapter() {
    let schema = Schema {
        paired: true,
        quality: true,
        headers: true,
        flags: false,
    };
    let fixture = Fixture::new();
    let records = contaminated_records(schema, 200, ADAPTER_R1.as_bytes(), ADAPTER_R2.as_bytes());
    let input = fixture.input(schema, &records, BLOCK);
    let output = fixture.path("out.cbq");
    let report = fixture.path("report.json");

    bqc_ok(&[
        "adapter",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--auto-detect",
        "--report",
        report.to_str().unwrap(),
    ]);

    // The 40 contaminated pairs are cut to their inserts.
    let (_, out_records) = read_cbq(&output);
    assert_eq!(out_records.len(), 200);
    for (index, record) in out_records.iter().enumerate() {
        if index % 5 == 0 {
            assert!(
                record.s_seq.len() < records[index].s_seq.len(),
                "record {index} was not trimmed"
            );
        }
    }

    let report: serde_json::Value = serde_json::from_str(&text(&report)).unwrap();
    let detection = &report["adapter"]["detection"];
    assert_eq!(detection["database_version"], 1);
    assert_eq!(detection["r1"]["decision"], "confident", "{detection}");
    assert_eq!(detection["r1"]["recommended_name"], "illumina-truseq");
    assert_eq!(
        detection["r1"]["recommended_sequence"].as_str().unwrap(),
        ADAPTER_R1
    );
    assert_eq!(detection["r2"]["decision"], "confident");
    assert_eq!(
        detection["r2"]["recommended_sequence"].as_str().unwrap(),
        ADAPTER_R2
    );
    // The report carries the evidence, not just the answer.
    let leader = &detection["r1"]["candidates"][0];
    assert_eq!(leader["confidence"], "high");
    assert!(
        leader["support_fraction"].as_f64().unwrap() >= 0.01,
        "{leader}"
    );
    // The detected sequences are part of the resolved configuration.
    let adapters = report["configuration"]["workflow"]["adapter"]["r1"]
        .as_array()
        .unwrap();
    assert!(
        adapters
            .iter()
            .any(|a| a["sequence"].as_str() == Some(ADAPTER_R1)),
        "{adapters:?}"
    );
}

#[test]
fn auto_detection_with_indels_is_not_lost_by_candidate_seeding() {
    // Every adapter has an insertion inside the fixed 5' seed used by the old
    // shortlist. The indel-aware matcher accepts it, so candidate generation
    // must still deliver the known entry to verification.
    let fixture = Fixture::new();
    let mut sequences = Sequences::new(0x1ADE_15E1);
    let records: Vec<Record> = (0..400)
        .map(|index| {
            let insert = sequences.sequence(30 + index % 20);
            let mut damaged = ADAPTER_R1.as_bytes().to_vec();
            damaged.insert(6, b'T');
            let read = [insert.as_slice(), damaged.as_slice()].concat();
            let quality = vec![b'I'; read.len()];
            Record::new(SINGLE, index, &read, b"").with_quality(&quality, b"")
        })
        .collect();
    let input = fixture.input(SINGLE, &records, BLOCK);
    let output = fixture.path("trimmed.cbq");
    let report = fixture.path("report.json");

    bqc_ok(&[
        "adapter",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--auto-detect",
        "--allow-indels",
        "--report",
        report.to_str().unwrap(),
    ]);

    let report = json(&report);
    let detection = &report["adapter"]["detection"]["r1"];
    assert_eq!(detection["decision"], "confident", "{detection}");
    assert_eq!(detection["recommended_sequence"], ADAPTER_R1);
    assert!(
        detection["candidates"][0]["indel_matches"]
            .as_u64()
            .unwrap()
            > 0,
        "{detection}"
    );

    let (_, trimmed) = read_cbq(&output);
    for (index, record) in trimmed.iter().enumerate() {
        assert_eq!(record.s_seq.len(), 30 + index % 20, "record {index}");
    }
}

#[test]
fn auto_detection_completes_the_known_partner_mate() {
    let schema = Schema {
        paired: true,
        quality: true,
        headers: true,
        flags: false,
    };
    let fixture = Fixture::new();
    // Only R1 is contaminated; R2 reads are clean throughout.
    let mut sequences = Sequences::new(0x00DE_7EC7);
    let records: Vec<Record> = (0..200)
        .map(|index| {
            let r1 = if index % 5 == 0 {
                let insert = sequences.sequence(20 + index % 30);
                [insert.as_slice(), &ADAPTER_R1.as_bytes()[..20]].concat()
            } else {
                sequences.sequence(40 + index % 30)
            };
            let r2 = sequences.sequence(40 + index % 30);
            let q1 = vec![b'I'; r1.len()];
            let q2 = vec![b'I'; r2.len()];
            Record::new(schema, index, &r1, &r2).with_quality(&q1, &q2)
        })
        .collect();
    let input = fixture.input(schema, &records, BLOCK);
    let output = fixture.path("out.cbq");
    let report = fixture.path("report.json");

    bqc_ok(&[
        "adapter",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--auto-detect",
        "--report",
        report.to_str().unwrap(),
    ]);

    let report: serde_json::Value = serde_json::from_str(&text(&report)).unwrap();
    let detection = &report["adapter"]["detection"];
    assert_eq!(detection["r1"]["decision"], "confident");
    assert_eq!(detection["r1"]["recommended_name"], "illumina-truseq");
    // The mates are inferred independently, and R2 carries no adapter, so no
    // sequence is invented for it. Assuming both mates share a chemistry would
    // configure R2 with a sequence nothing in the data supports.
    assert_eq!(detection["r2"]["decision"], "inconclusive", "{detection}");
    assert!(detection["r2"]["recommended_sequence"].is_null());
    assert!(
        report["configuration"]["workflow"]["adapter"]["r2"]
            .as_array()
            .unwrap()
            .is_empty(),
        "no adapter may be configured for a mate with no evidence"
    );

    // R2 is therefore untouched, and every record still passes through in order.
    let (_, out_records) = read_cbq(&output);
    assert_eq!(out_records.len(), 200);
    for (index, record) in out_records.iter().enumerate() {
        assert_eq!(
            record.x_seq.as_ref().unwrap().len(),
            records[index].x_seq.as_ref().unwrap().len(),
            "R2 of record {index} must be untouched"
        );
    }
}

#[test]
fn auto_detection_assembles_a_consensus_for_an_unknown_adapter() {
    let schema = Schema {
        paired: false,
        quality: true,
        headers: true,
        flags: false,
    };
    let fixture = Fixture::new();
    // An adapter that is not in the known library.
    let unknown = b"GTCGATCGTACGGCATCCGATCGTACG";
    let mut sequences = Sequences::new(0xC0C0A);
    let records: Vec<Record> = (0..200)
        .map(|index| {
            let seq = if index % 3 == 0 {
                // 33% contamination.
                let insert = sequences.sequence(20 + index % 30);
                [insert.as_slice(), &unknown[..20]].concat()
            } else {
                sequences.sequence(40 + index % 30)
            };
            let qual = vec![b'I'; seq.len()];
            Record::new(schema, index, &seq, b"").with_quality(&qual, b"")
        })
        .collect();
    let input = fixture.input(schema, &records, BLOCK);
    let output = fixture.path("out.cbq");
    let report = fixture.path("report.json");

    bqc_ok(&[
        "adapter",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--auto-detect",
        "--report",
        report.to_str().unwrap(),
    ]);

    let report: serde_json::Value = serde_json::from_str(&text(&report)).unwrap();
    let detection = &report["adapter"]["detection"];
    assert_eq!(detection["r1"]["decision"], "confident", "{detection}");
    let leader = &detection["r1"]["candidates"][0];
    // Nothing in the library matches, so the evidence is k-mer assembly alone.
    assert_eq!(leader["evidence_sources"][0], "kmer_consensus");
    assert!(leader["known_name"].is_null());
    assert_eq!(
        detection["r1"]["recommended_sequence"].as_str().unwrap(),
        String::from_utf8_lossy(&unknown[..20]),
        "the consensus reproduces the unknown adapter from its 5' start"
    );

    let (_, out_records) = read_cbq(&output);
    for (index, record) in out_records.iter().enumerate() {
        if index % 3 == 0 {
            assert!(
                record.s_seq.len() < records[index].s_seq.len(),
                "record {index} was not trimmed"
            );
        }
    }
}

#[test]
fn auto_detection_passes_through_when_nothing_clears_the_gates() {
    let schema = Schema {
        paired: false,
        quality: true,
        headers: true,
        flags: false,
    };
    let fixture = Fixture::new();
    let records = simple_records(schema, 100, 50);
    let input = fixture.input(schema, &records, BLOCK);
    let output = fixture.path("out.cbq");
    let report = fixture.path("report.json");

    bqc_ok(&[
        "adapter",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--auto-detect",
        "--report",
        report.to_str().unwrap(),
    ]);
    assert_eq!(read_cbq(&output).1.len(), 100);
    let report: serde_json::Value = serde_json::from_str(&text(&report)).unwrap();
    assert_eq!(
        report["adapter"]["detection"]["r1"]["decision"],
        "inconclusive"
    );
    assert_eq!(report["adapter"]["r1_reads_trimmed"], 0);

    // Detection thresholds need the flag itself.
    let error = bqc(&[
        "adapter",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--adapter-r1",
        ADAPTER_R1,
        "--detect-sample-size",
        "10",
    ])
    .unwrap_err();
    assert!(
        format!("{error}").contains("require --auto-detect"),
        "{error}"
    );
}

// ----------------------------------------------------------------- correction

/// The overlap geometry `overlapping_pairs` builds, for a 60 base read: R1
/// covers insert[0..60], R2 covers insert[30..90], so the overlap is R1[30..60]
/// and R2[30..60].
const CORRECTION_READ: usize = 60;

fn correction_plan() -> Vec<Planted> {
    // One of each class, repeated so several blocks see every case.
    (0..40)
        .map(|index| match index % 4 {
            0 => Planted::Clean,
            1 => Planted::R2Doubtful(40),
            2 => Planted::R1Doubtful(45),
            _ => Planted::R1Confident(50),
        })
        .collect()
}

#[test]
// One run's every output is checked together: sequences, qualities, metadata,
// report and log. Splitting it would only scatter the assertions.
#[allow(clippy::too_many_lines)]
fn correction_fixes_doubtful_bases_from_the_confident_mate() {
    let schema = Schema {
        paired: true,
        quality: true,
        headers: true,
        flags: true,
    };
    let fixture = Fixture::new();
    let (records, offset) = overlapping_pairs(schema, CORRECTION_READ, &correction_plan());
    let input = fixture.input(schema, &records, BLOCK);
    let output = fixture.path("corrected.cbq");
    let log = fixture.path("corrections.tsv");
    let report = fixture.path("report.json");

    bqc_ok(&[
        "correct",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--correction-log",
        log.to_str().unwrap(),
        "--report",
        report.to_str().unwrap(),
        "-T",
        "2",
    ]);

    let (out_schema, corrected) = read_cbq(&output);
    assert_eq!(out_schema, schema, "the schema is preserved");
    assert_eq!(corrected.len(), records.len(), "every pair is written");

    for (index, (before, after)) in records.iter().zip(&corrected).enumerate() {
        assert_eq!(after.s_header, before.s_header, "headers preserved");
        assert_eq!(after.flag, before.flag, "flags preserved");
        assert_eq!(after.s_seq.len(), before.s_seq.len(), "R1 length unchanged");
        let r2_before = before.x_seq.as_ref().unwrap();
        let r2_after = after.x_seq.as_ref().unwrap();
        assert_eq!(r2_after.len(), r2_before.len(), "R2 length unchanged");
        assert_eq!(after.s_qual.as_ref().unwrap().len(), after.s_seq.len());

        match index % 4 {
            // Clean and confident-error pairs must come back byte-identical.
            0 | 3 => assert_eq!(after, before, "pair {index} must not change"),
            // R2 was doubtful: R1 donates, so R2[40] becomes the complement of
            // the aligned R1 base and takes R1's quality byte.
            1 => {
                let position2 = 40;
                let position1 = r2_before.len() - 1 - position2 + offset;
                assert_eq!(r2_after[position2], complement(after.s_seq[position1]));
                assert_ne!(r2_after[position2], r2_before[position2], "R2 corrected");
                assert_eq!(
                    after.x_qual.as_ref().unwrap()[position2],
                    after.s_qual.as_ref().unwrap()[position1],
                    "donor quality byte copied exactly"
                );
                assert_eq!(after.s_seq, before.s_seq, "R1 untouched");
            }
            // R1 was doubtful at 45, so R2 donates.
            _ => {
                let position1 = 45;
                let position2 = r2_before.len() - 1 - (position1 - offset);
                assert_eq!(after.s_seq[position1], complement(r2_after[position2]));
                assert_ne!(after.s_seq[position1], before.s_seq[position1]);
                assert_eq!(
                    after.s_qual.as_ref().unwrap()[position1],
                    after.x_qual.as_ref().unwrap()[position2]
                );
                assert_eq!(r2_after, r2_before, "R2 untouched");
            }
        }
    }

    // Aggregate statistics distinguish pairs, reads and bases.
    let report: serde_json::Value = serde_json::from_str(&text(&report)).unwrap();
    let correction = &report["correction"];
    assert_eq!(correction["pairs_examined"], 40);
    assert_eq!(correction["pairs_with_overlap"], 40);
    assert_eq!(correction["pairs_with_mismatches"], 30);
    assert_eq!(correction["corrected_pairs"], 20);
    assert_eq!(correction["corrected_r1_bases"], 10);
    assert_eq!(correction["corrected_r2_bases"], 10);
    assert_eq!(correction["corrected_r1_reads"], 10);
    assert_eq!(correction["corrected_r2_reads"], 10);
    assert_eq!(correction["corrected_pairs_both_mates"], 0);
    assert_eq!(correction["unresolved_mismatches"], 10);
    assert_eq!(correction["noncanonical_donors_skipped"], 0);
    // The two large tables live in JSON only.
    assert_eq!(correction["corrections_per_pair"][0]["corrections"], 1);
    assert_eq!(correction["corrections_per_pair"][0]["pairs"], 20);
    assert!(
        !correction["substitutions"].as_array().unwrap().is_empty(),
        "substitutions are captured before mutation"
    );
    let substituted: u64 = correction["substitutions"]
        .as_array()
        .unwrap()
        .iter()
        .map(|entry| entry["bases"].as_u64().unwrap())
        .sum();
    assert_eq!(
        substituted, 20,
        "the matrix reconciles with corrected bases"
    );
    assert_eq!(
        correction["corrected_pairs_by_disposition"][0]["disposition"],
        "accepted"
    );
    assert_eq!(
        report["configuration"]["stage_order"],
        serde_json::json!(["correct"])
    );
    assert_eq!(
        report["configuration"]["workflow"]["correction"]["donor_quality"],
        30
    );
    assert_eq!(
        report["configuration"]["workflow"]["correction"]["recipient_quality"],
        14
    );

    // The read-level log has one row per corrected pair, in record order.
    let rendered = text(&log);
    let rows: Vec<&str> = rendered.lines().skip(1).map(str::trim_end).collect();
    assert_eq!(rows.len(), 20);
    let indices: Vec<u64> = rows
        .iter()
        .map(|row| row.split('\t').next().unwrap().parse().unwrap())
        .collect();
    let mut sorted = indices.clone();
    sorted.sort_unstable();
    assert_eq!(indices, sorted, "log rows are ordered by record index");
    assert!(rows[0].contains("\taccepted"), "{}", rows[0]);
    assert!(
        rows[0].contains("read_1/1"),
        "headers are recorded: {}",
        rows[0]
    );
}

#[test]
fn correction_requires_paired_input_with_stored_qualities() {
    let fixture = Fixture::new();
    let single = Schema {
        paired: false,
        quality: true,
        headers: true,
        flags: false,
    };
    let input = fixture.input(single, &simple_records(single, 10, 60), BLOCK);
    let output = fixture.path("out.cbq");
    let error = bqc(&[
        "correct",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
    ])
    .unwrap_err();
    assert!(
        format!("{error}").contains("requires a paired input"),
        "{error}"
    );
    assert!(!output.exists());

    let unqualified = Schema {
        paired: true,
        quality: false,
        headers: true,
        flags: false,
    };
    let fixture = Fixture::new();
    let input = fixture.input(unqualified, &simple_records(unqualified, 10, 60), BLOCK);
    let output = fixture.path("out.cbq");
    let error = bqc(&[
        "correct",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
    ])
    .unwrap_err();
    assert!(format!("{error}").contains("no quality column"), "{error}");
    assert!(!output.exists());
}

#[test]
fn correction_thresholds_are_validated() {
    let schema = Schema {
        paired: true,
        quality: true,
        headers: true,
        flags: false,
    };
    let fixture = Fixture::new();
    let (records, _) = overlapping_pairs(schema, CORRECTION_READ, &correction_plan());
    let input = fixture.input(schema, &records, BLOCK);
    let output = fixture.path("out.cbq");

    let error = bqc(&[
        "correct",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--donor-quality",
        "14",
        "--recipient-quality",
        "14",
    ])
    .unwrap_err();
    assert!(
        format!("{error}").contains("must be greater than --recipient-quality"),
        "{error}"
    );

    // Threshold options without the stage are a mistake, not a silent no-op.
    let error = bqc(&[
        "workflow",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--donor-quality",
        "30",
        "--min-length",
        "10",
    ])
    .unwrap_err();
    assert!(
        format!("{error}").contains("require --correction"),
        "{error}"
    );
    assert!(!output.exists());
}

#[test]
fn corrected_qualities_are_visible_to_filtering() {
    // A doubtful base drags a read's mean quality below the threshold; once the
    // donor's quality is copied across, the same read passes. This is the reason
    // correction has to run before trimming and filtering.
    let schema = Schema {
        paired: true,
        quality: true,
        headers: true,
        flags: false,
    };
    let fixture = Fixture::new();
    let plan: Vec<Planted> = (0..20)
        .map(|index| {
            if index % 2 == 0 {
                Planted::Clean
            } else {
                Planted::R2Doubtful(40)
            }
        })
        .collect();
    let (records, _) = overlapping_pairs(schema, CORRECTION_READ, &plan);
    let input = fixture.input(schema, &records, BLOCK);

    let without = fixture.path("without.cbq");
    bqc_ok(&[
        "filter",
        input.to_str().unwrap(),
        "-o",
        without.to_str().unwrap(),
        "--min-mean-quality",
        "35",
    ]);
    assert_eq!(
        read_cbq(&without).1.len(),
        10,
        "doubtful pairs are rejected"
    );

    let with = fixture.path("with.cbq");
    bqc_ok(&[
        "workflow",
        input.to_str().unwrap(),
        "-o",
        with.to_str().unwrap(),
        "--correction",
        "--min-mean-quality",
        "35",
    ]);
    assert_eq!(
        read_cbq(&with).1.len(),
        20,
        "correction rescues them before the filter runs"
    );
}

#[test]
fn standalone_correction_matches_the_workflow_step() {
    let schema = Schema {
        paired: true,
        quality: true,
        headers: true,
        flags: true,
    };
    let fixture = Fixture::new();
    let (records, _) = overlapping_pairs(schema, CORRECTION_READ, &correction_plan());
    let input = fixture.input(schema, &records, 4096);

    let standalone = fixture.path("standalone.cbq");
    let standalone_log = fixture.path("standalone.tsv");
    bqc_ok(&[
        "correct",
        input.to_str().unwrap(),
        "-o",
        standalone.to_str().unwrap(),
        "--correction-log",
        standalone_log.to_str().unwrap(),
        "--donor-quality",
        "25",
        "--recipient-quality",
        "10",
    ]);

    let fused = fixture.path("fused.cbq");
    let fused_log = fixture.path("fused.tsv");
    bqc_ok(&[
        "workflow",
        input.to_str().unwrap(),
        "-o",
        fused.to_str().unwrap(),
        "--steps",
        "correct",
        "--correction",
        "--correction-log",
        fused_log.to_str().unwrap(),
        "--donor-quality",
        "25",
        "--recipient-quality",
        "10",
    ]);

    assert_same_file(&standalone, &fused);
    assert_eq!(text(&standalone_log), text(&fused_log));
}

#[test]
fn correction_output_is_identical_at_any_thread_count() {
    let schema = Schema {
        paired: true,
        quality: true,
        headers: true,
        flags: true,
    };
    let fixture = Fixture::new();
    let plan: Vec<Planted> = (0..400)
        .map(|index| match index % 4 {
            0 => Planted::Clean,
            1 => Planted::R2Doubtful(40),
            2 => Planted::R1Doubtful(45),
            _ => Planted::R1Confident(50),
        })
        .collect();
    let (records, _) = overlapping_pairs(schema, CORRECTION_READ, &plan);
    // Small blocks so the parallel engine has several chunks to reorder.
    let input = fixture.input(schema, &records, 4096);
    assert!(count_blocks(&input) > 2);

    let mut runs = Vec::new();
    for threads in [1, 3, 8] {
        for detail in ["reads", "bases"] {
            let output = fixture.path(&format!("out-{threads}-{detail}.cbq"));
            let log = fixture.path(&format!("log-{threads}-{detail}.tsv"));
            bqc_ok(&[
                "correct",
                input.to_str().unwrap(),
                "-o",
                output.to_str().unwrap(),
                "--correction-log",
                log.to_str().unwrap(),
                "--correction-log-detail",
                detail,
                "-T",
                &threads.to_string(),
            ]);
            runs.push((threads, detail.to_string(), output, log));
        }
    }
    for detail in ["reads", "bases"] {
        let reference = runs
            .iter()
            .find(|(threads, kind, ..)| *threads == 1 && kind == detail)
            .unwrap();
        for candidate in runs
            .iter()
            .filter(|(threads, kind, ..)| *threads != 1 && kind == detail)
        {
            assert_same_file(&reference.2, &candidate.2);
            assert_eq!(
                text(&reference.3),
                text(&candidate.3),
                "{detail} log differs at T={}",
                candidate.0
            );
        }
    }
    // The base-level log has one row per corrected base, ordered by record then
    // mate then position.
    let base_log = &runs
        .iter()
        .find(|(threads, kind, ..)| *threads == 1 && kind == "bases")
        .unwrap()
        .3;
    let rendered = text(base_log);
    let rows: Vec<&str> = rendered.lines().skip(1).collect();
    assert_eq!(rows.len(), 200, "one row per corrected base");
    let keys: Vec<(u64, String, usize)> = rows
        .iter()
        .map(|row| {
            let mut fields = row.split('\t');
            let index = fields.next().unwrap().parse().unwrap();
            let mate = fields.next().unwrap().to_string();
            let position = fields.next().unwrap().parse().unwrap();
            (index, mate, position)
        })
        .collect();
    let mut sorted = keys.clone();
    sorted.sort();
    assert_eq!(keys, sorted, "base log order is record, mate, position");
}

#[test]
fn failed_mode_original_writes_uncorrected_pairs() {
    let schema = Schema {
        paired: true,
        quality: true,
        headers: true,
        flags: false,
    };
    let fixture = Fixture::new();
    let plan: Vec<Planted> = (0..20).map(|_| Planted::R2Doubtful(40)).collect();
    let (records, _) = overlapping_pairs(schema, CORRECTION_READ, &plan);
    let input = fixture.input(schema, &records, BLOCK);

    // No pair can reach mean Q36, so everything is rejected either way.
    for (mode, expect_corrections) in [("original", false), ("processed", true)] {
        let accepted = fixture.path(&format!("acc-{mode}.cbq"));
        let failed = fixture.path(&format!("failed-{mode}.cbq"));
        bqc_ok(&[
            "workflow",
            input.to_str().unwrap(),
            "-o",
            accepted.to_str().unwrap(),
            "--correction",
            "--min-mean-quality",
            "36",
            "--failed",
            failed.to_str().unwrap(),
            "--failed-mode",
            mode,
        ]);
        let (_, rejected) = read_cbq(&failed);
        assert_eq!(rejected.len(), 20, "mode {mode}");
        let changed = rejected
            .iter()
            .zip(&records)
            .any(|(after, before)| after.x_seq != before.x_seq);
        assert_eq!(
            changed,
            expect_corrections,
            "failed-mode {mode} must {} corrections",
            if expect_corrections { "keep" } else { "drop" }
        );
    }
}

#[test]
fn correction_logs_escape_hostile_headers() {
    let schema = Schema {
        paired: true,
        quality: true,
        headers: true,
        flags: false,
    };
    let fixture = Fixture::new();
    let (mut records, _) = overlapping_pairs(schema, CORRECTION_READ, &[Planted::R2Doubtful(40)]);
    records[0].s_header = Some(b"tab\there\nnewline\\slash".to_vec());
    records[0].x_header = Some(b"carriage\rreturn".to_vec());
    let input = fixture.input(schema, &records, BLOCK);
    let output = fixture.path("out.cbq");
    let log = fixture.path("log.tsv");

    bqc_ok(&[
        "correct",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--correction-log",
        log.to_str().unwrap(),
    ]);

    let rendered = text(&log);
    let rows: Vec<&str> = rendered.lines().collect();
    assert_eq!(rows.len(), 2, "a header must not spill into extra rows");
    let fields: Vec<&str> = rows[1].split('\t').collect();
    assert_eq!(
        fields.len(),
        10,
        "a header must not spill into extra columns"
    );
    assert_eq!(fields[1], "tab\\there\\nnewline\\\\slash");
    assert_eq!(fields[2], "carriage\\rreturn");
}

#[test]
fn correction_applies_to_surviving_orphan_mates() {
    let schema = Schema {
        paired: true,
        quality: true,
        headers: true,
        flags: false,
    };
    let fixture = Fixture::new();
    // R1 is doubtful inside the overlap and R2 is short enough to fail a length
    // filter, so R1 survives alone — carrying its correction.
    let plan: Vec<Planted> = (0..8).map(|_| Planted::R1Doubtful(45)).collect();
    let (mut records, _) = overlapping_pairs(schema, CORRECTION_READ, &plan);
    for record in &mut records {
        // The overlap is R2[30..60]; ambiguous bases before it fail an N filter
        // without disturbing the alignment or the donor evidence.
        let sequence = record.x_seq.as_mut().unwrap();
        sequence[..6].fill(b'N');
    }
    let input = fixture.input(schema, &records, BLOCK);
    let accepted = fixture.path("clean.cbq");
    let prefix = fixture.path("surviving");

    bqc_ok(&[
        "workflow",
        input.to_str().unwrap(),
        "-o",
        accepted.to_str().unwrap(),
        "--correction",
        "--max-n",
        "2",
        "--pair-policy",
        "orphan",
        "--orphan-prefix",
        prefix.to_str().unwrap(),
    ]);

    let (orphan_schema, orphans) = read_cbq(&fixture.path("surviving.R1.cbq"));
    assert_eq!(orphan_schema, schema.unpaired());
    assert_eq!(orphans.len(), 8, "every R1 survives alone");
    for (orphan, before) in orphans.iter().zip(&records) {
        assert_ne!(
            orphan.s_seq[45], before.s_seq[45],
            "the surviving mate keeps its correction"
        );
    }
}

// ------------------------------------------------------ linked segmentation

// Deliberately non-repetitive: a self-similar flank can legitimately match
// earlier inside the insert, which is correct behaviour but makes a fixture's
// expected geometry ambiguous.
const LINKED_5P: &str = "AGGTCAGTCTAC";
const LINKED_3P: &str = "CTTACGGATCCA";

/// Reads shaped `prefix + 5' flank + insert + 3' flank + suffix`, with a few
/// reads deliberately missing a flank.
fn linked_records(schema: Schema, count: usize) -> Vec<Record> {
    let mut sequences = Sequences::new(0x0114_ED01);
    (0..count)
        .map(|index| {
            let mut insert: Vec<u8> = sequences
                .sequence(30)
                .into_iter()
                .map(|base| if base == b'N' { b'A' } else { base })
                .collect();
            // An accidental flank prefix inside the insert would move the true
            // boundary, so perturb one base if that happens.
            for flank in [LINKED_5P, LINKED_3P] {
                let seed = &flank.as_bytes()[..8];
                while let Some(at) = insert.windows(8).position(|window| window == seed) {
                    insert[at] = if insert[at] == b'A' { b'C' } else { b'A' };
                }
            }
            let read = match index % 4 {
                // Both flanks, with two leading bases before the 5' flank.
                0 | 1 => [
                    b"AT".as_slice(),
                    LINKED_5P.as_bytes(),
                    &insert,
                    LINKED_3P.as_bytes(),
                    b"CCC",
                ]
                .concat(),
                // 3' flank missing.
                2 => [LINKED_5P.as_bytes(), &insert].concat(),
                // Neither flank.
                _ => insert.clone(),
            };
            let quality = vec![b'I'; read.len()];
            Record::new(schema, index, &read, &read).with_quality(&quality, &quality)
        })
        .collect()
}

#[test]
fn linked_adapters_retain_the_insert_between_the_flanks() {
    let schema = Schema {
        paired: false,
        quality: true,
        headers: true,
        flags: true,
    };
    let fixture = Fixture::new();
    let records = linked_records(schema, 40);
    let input = fixture.input(schema, &records, BLOCK);
    let output = fixture.path("inserts.cbq");
    let report = fixture.path("report.json");

    bqc_ok(&[
        "adapter",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--linked-5p-r1",
        LINKED_5P,
        "--linked-3p-r1",
        LINKED_3P,
        "--linked-min-insert-length",
        "10",
        "--report",
        report.to_str().unwrap(),
        "-T",
        "2",
    ]);

    let (out_schema, out) = read_cbq(&output);
    assert_eq!(out_schema, schema, "the schema is preserved");
    assert_eq!(out.len(), records.len(), "one record in, one record out");
    for (index, (before, after)) in records.iter().zip(&out).enumerate() {
        assert_eq!(after.s_header, before.s_header, "headers preserved");
        assert_eq!(after.flag, before.flag, "flags preserved");
        assert_eq!(
            after.s_qual.as_ref().unwrap().len(),
            after.s_seq.len(),
            "quality tracks sequence"
        );
        match index % 4 {
            // Both flanks: exactly the 30 base insert survives.
            0 | 1 => {
                assert_eq!(after.s_seq.len(), 30, "record {index}");
                let start = 2 + LINKED_5P.len();
                assert_eq!(
                    after.s_seq,
                    before.s_seq[start..start + 30],
                    "record {index}"
                );
                assert_eq!(
                    after.s_qual.as_ref().unwrap(),
                    &before.s_qual.as_ref().unwrap()[start..start + 30]
                );
            }
            // No 3' flank, `both` required: the read falls through unchanged,
            // because the default unmatched policy is to continue.
            _ => assert_eq!(after.s_seq, before.s_seq, "record {index} must not change"),
        }
    }

    let report: serde_json::Value = serde_json::from_str(&text(&report)).unwrap();
    let linked = &report["linked"];
    assert_eq!(linked["both_adapters"], 20);
    assert_eq!(linked["five_prime_only"], 10, "the 3' flank was missing");
    assert_eq!(linked["neither_adapter"], 10);
    assert_eq!(linked["invalid_order"], 0);
    assert_eq!(
        linked["leading_bases_removed"],
        20 * (2 + LINKED_5P.len()) as u64
    );
    assert_eq!(
        linked["trailing_bases_removed"],
        20 * (LINKED_3P.len() + 3) as u64
    );
    assert_eq!(linked["per_definition"][0]["name"], "linked1");
    assert_eq!(linked["per_definition"][0]["mate"], "R1");
    assert_eq!(linked["per_definition"][0]["matches"], 20);
    assert_eq!(linked["unmatched_policy"], "continue");
    assert_eq!(
        report["configuration"]["stage_order"],
        serde_json::json!(["linked"]),
        "no ordinary adapter sequence was configured, so only linked ran"
    );
}

#[test]
fn the_unmatched_policy_decides_what_happens_to_the_rest() {
    let schema = Schema {
        paired: false,
        quality: true,
        headers: true,
        flags: false,
    };
    let fixture = Fixture::new();
    let records = linked_records(schema, 40);
    let input = fixture.input(schema, &records, BLOCK);

    // `fail` routes unmatched reads through the ordinary rejection path.
    let accepted = fixture.path("accepted.cbq");
    let failed = fixture.path("failed.cbq");
    let reasons = fixture.path("reasons.tsv");
    bqc_ok(&[
        "workflow",
        input.to_str().unwrap(),
        "-o",
        accepted.to_str().unwrap(),
        "--linked-5p-r1",
        LINKED_5P,
        "--linked-3p-r1",
        LINKED_3P,
        "--linked-unmatched",
        "fail",
        "--min-length",
        "1",
        "--failed",
        failed.to_str().unwrap(),
        "--failed-reasons",
        reasons.to_str().unwrap(),
    ]);
    assert_eq!(
        read_cbq(&accepted).1.len(),
        20,
        "only matched reads survive"
    );
    assert_eq!(read_cbq(&failed).1.len(), 20);
    let sidecar = text(&reasons);
    assert!(
        sidecar.contains("LINKED_UNMATCHED"),
        "the dedicated reason is recorded: {}",
        sidecar.lines().nth(1).unwrap_or_default()
    );

    // `keep` leaves unmatched reads untouched and skips ordinary adapter work.
    let kept = fixture.path("kept.cbq");
    bqc_ok(&[
        "adapter",
        input.to_str().unwrap(),
        "-o",
        kept.to_str().unwrap(),
        "--linked-5p-r1",
        LINKED_5P,
        "--linked-3p-r1",
        LINKED_3P,
        "--linked-unmatched",
        "keep",
    ]);
    let (_, out) = read_cbq(&kept);
    for (index, (before, after)) in records.iter().zip(&out).enumerate() {
        if index % 4 >= 2 {
            assert_eq!(after.s_seq, before.s_seq, "record {index} kept verbatim");
        }
    }
}

#[test]
fn linked_definitions_come_from_a_configuration_file() {
    let schema = Schema {
        paired: true,
        quality: true,
        headers: true,
        flags: false,
    };
    let fixture = Fixture::new();
    let records = linked_records(schema, 20);
    let input = fixture.input(schema, &records, BLOCK);
    let config = fixture.path("amplicon.toml");
    std::fs::write(
        &config,
        format!(
            r#"
[[adapter.linked_r1]]
name = "amplicon"
five_prime = "{LINKED_5P}"
three_prime = "{LINKED_3P}"
require = "both"
max_five_prime_offset = 3
minimum_insert_length = 20

[[adapter.linked_r2]]
name = "amplicon"
five_prime = "{LINKED_5P}"
three_prime = "{LINKED_3P}"
"#
        ),
    )
    .unwrap();
    let output = fixture.path("out.cbq");
    let report = fixture.path("report.json");
    bqc_ok(&[
        "workflow",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--config",
        config.to_str().unwrap(),
        "--report",
        report.to_str().unwrap(),
    ]);
    let (_, out) = read_cbq(&output);
    assert_eq!(out.len(), 20);
    assert_eq!(out[0].s_seq.len(), 30, "the insert survives on R1");
    assert_eq!(
        out[0].x_seq.as_ref().unwrap().len(),
        30,
        "and independently on R2"
    );
    let report: serde_json::Value = serde_json::from_str(&text(&report)).unwrap();
    assert_eq!(report["linked"]["per_definition"][0]["name"], "amplicon");
    assert_eq!(
        report["configuration"]["workflow"]["linked"]["r1"][0]["minimum_insert_length"],
        20
    );
}

#[test]
fn linked_output_is_identical_at_any_thread_count() {
    let schema = Schema {
        paired: true,
        quality: true,
        headers: true,
        flags: true,
    };
    let fixture = Fixture::new();
    let records = linked_records(schema, 400);
    let input = fixture.input(schema, &records, 4096);
    assert!(count_blocks(&input) > 2);
    let mut outputs = Vec::new();
    for threads in ["1", "3", "8"] {
        let output = fixture.path(&format!("out{threads}.cbq"));
        bqc_ok(&[
            "workflow",
            input.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
            "--linked-5p-r1",
            LINKED_5P,
            "--linked-3p-r1",
            LINKED_3P,
            "--linked-5p-r2",
            LINKED_5P,
            "--linked-3p-r2",
            LINKED_3P,
            "--linked-unmatched",
            "keep",
            "--quality-tail",
            "20",
            "-T",
            threads,
        ]);
        outputs.push(output);
    }
    for other in &outputs[1..] {
        assert_same_file(&outputs[0], other);
    }
}

#[test]
fn mate_two_linked_definitions_require_a_paired_input() {
    let schema = Schema {
        paired: false,
        quality: true,
        headers: true,
        flags: false,
    };
    let fixture = Fixture::new();
    let records = linked_records(schema, 8);
    let input = fixture.input(schema, &records, BLOCK);
    let output = fixture.path("out.cbq");
    let error = bqc(&[
        "adapter",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--linked-5p-r2",
        LINKED_5P,
        "--linked-3p-r2",
        LINKED_3P,
    ])
    .unwrap_err();
    assert!(
        format!("{error}").contains("require a paired input"),
        "{error}"
    );
    assert!(!output.exists());

    // A definition needs both sides.
    let error = bqc(&[
        "adapter",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--linked-5p-r1",
        LINKED_5P,
    ])
    .unwrap_err();
    assert!(
        format!("{error}").contains("needs both a 5' and a 3' sequence"),
        "{error}"
    );
}

// ------------------------------------------------------------------ segmenting

/// Two delimiters that share no prefix with each other or with random sequence.
const DELIMITER_A: &str = "AGGTCAGTCTACAT";
const DELIMITER_B: &str = "CTTACGGATCCAGG";

/// Reads carrying zero, one or two internal delimiters, plus one that begins
/// with a delimiter — the case that produces an empty prefix fragment.
fn segment_records(schema: Schema, count: usize) -> Vec<Record> {
    let mut sequences = Sequences::new(0x5E6_1234);
    (0..count)
        .map(|index| {
            let mut piece = |length: usize| {
                let mut piece = sequences.sequence(length);
                // An accidental delimiter prefix inside a piece would add a cut
                // the test did not plant, so perturb one base if that happens.
                for delimiter in [DELIMITER_A, DELIMITER_B] {
                    let seed = &delimiter.as_bytes()[..8];
                    while let Some(at) = piece.windows(8).position(|window| window == seed) {
                        piece[at] = if piece[at] == b'A' { b'C' } else { b'A' };
                    }
                }
                piece
            };
            let read = match index % 4 {
                // No delimiter: the whole read is one fragment.
                0 => piece(30),
                // One delimiter: a prefix and a suffix.
                1 => [piece(30).as_slice(), DELIMITER_A.as_bytes(), &piece(20)].concat(),
                // Two delimiters: prefix, middle, suffix.
                2 => [
                    piece(30).as_slice(),
                    DELIMITER_A.as_bytes(),
                    &piece(20),
                    DELIMITER_B.as_bytes(),
                    &piece(25),
                ]
                .concat(),
                // A leading delimiter, so the prefix fragment is empty.
                _ => [DELIMITER_A.as_bytes(), piece(30).as_slice()].concat(),
            };
            let quality = sequences.quality(read.len());
            Record::new(schema, index, &read, &read).with_quality(&quality, &quality)
        })
        .collect()
}

/// The fragments each planted read should produce, as `(start, end)` pairs.
fn expected_fragments(record: &Record, index: usize) -> Vec<(usize, usize)> {
    let length = record.s_seq.len();
    let a = DELIMITER_A.len();
    match index % 4 {
        0 => vec![(0, length)],
        1 => vec![(0, 30), (30 + a, length)],
        2 => vec![
            (0, 30),
            (30 + a, 50 + a),
            (50 + a + DELIMITER_B.len(), length),
        ],
        _ => vec![(a, length)],
    }
}

fn delimiter_fasta(path: &Path) {
    std::fs::write(
        path,
        format!(">left\n{DELIMITER_A}\n>right\n{DELIMITER_B}\n"),
    )
    .expect("write delimiter fasta");
}

/// Checks every emitted fragment against the source record it came from, and
/// against its own provenance row.
fn verify_fragments(
    records: &[Record],
    expected: &[(usize, (usize, usize))],
    out: &[Record],
    rows: &[Vec<String>],
) {
    // Fragments arrive ordered by source record, then by segment index.
    let mut segment_index = 0usize;
    let mut previous = 0usize;
    for (position, ((source, (start, end)), fragment)) in expected.iter().zip(out).enumerate() {
        let record = &records[*source];
        segment_index = if *source == previous {
            segment_index
        } else {
            0
        };
        previous = *source;
        assert_eq!(
            fragment.s_seq,
            &record.s_seq[*start..*end],
            "fragment {segment_index} of record {source}"
        );
        assert_eq!(
            fragment.s_qual.as_ref().unwrap(),
            &record.s_qual.as_ref().unwrap()[*start..*end],
            "quality is sliced with the sequence"
        );
        assert_eq!(fragment.flag, record.flag, "flags describe the source read");
        // The header carries the provenance of the fragment.
        let header = String::from_utf8(fragment.s_header.clone().unwrap()).unwrap();
        let source_header = String::from_utf8(record.s_header.clone().unwrap()).unwrap();
        assert_eq!(
            header,
            format!("{source_header}|segment={segment_index}|span={start}-{end}"),
        );
        // And the sidecar says the same thing in columns.
        let row = &rows[position];
        assert_eq!(row[0], source.to_string(), "source_record_index");
        assert_eq!(row[1], segment_index.to_string(), "segment_index");
        assert_eq!(row[2], "R1");
        assert_eq!(row[3], start.to_string(), "start");
        assert_eq!(row[4], end.to_string(), "end");
        assert_eq!(row[5], (end - start).to_string(), "length");
        assert_eq!(row[8], source_header, "original_header");
        assert_eq!(row[9], "PASS");
        // Only the middle fragment of a two-delimiter read has both flanks.
        let both = row[6] != "." && row[7] != ".";
        assert_eq!(
            both,
            source % 4 == 2 && segment_index == 1,
            "flank columns of fragment {segment_index} of record {source}"
        );
        segment_index += 1;
    }
}

#[test]
fn segmentation_splits_reads_at_every_internal_delimiter() {
    let schema = Schema {
        paired: false,
        quality: true,
        headers: true,
        flags: true,
    };
    let fixture = Fixture::new();
    let records = segment_records(schema, 40);
    let input = fixture.input(schema, &records, BLOCK);
    let output = fixture.path("fragments.cbq");
    let fasta = fixture.path("delimiters.fa");
    let sidecar = fixture.path("segments.tsv");
    let report = fixture.path("report.json");
    delimiter_fasta(&fasta);

    bqc_ok(&[
        "segment",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--adapter-fasta",
        fasta.to_str().unwrap(),
        "--segments",
        sidecar.to_str().unwrap(),
        "--report",
        report.to_str().unwrap(),
        "-T",
        "2",
    ]);

    let (out_schema, out) = read_cbq(&output);
    assert_eq!(out_schema, schema, "the schema is preserved");
    let expected: Vec<(usize, (usize, usize))> = records
        .iter()
        .enumerate()
        .flat_map(|(index, record)| {
            expected_fragments(record, index)
                .into_iter()
                .map(move |span| (index, span))
        })
        .collect();
    assert_eq!(
        out.len(),
        expected.len(),
        "one output record per expected fragment"
    );

    let rows: Vec<Vec<String>> = text(&sidecar)
        .lines()
        .skip(1)
        .map(|line| line.split('\t').map(str::to_string).collect())
        .collect();
    assert_eq!(rows.len(), expected.len(), "one sidecar row per fragment");

    verify_fragments(&records, &expected, &out, &rows);

    // No delimiter may survive inside any fragment.
    for fragment in &out {
        for delimiter in [DELIMITER_A, DELIMITER_B] {
            assert!(
                !fragment
                    .s_seq
                    .windows(delimiter.len())
                    .any(|window| window == delimiter.as_bytes()),
                "a delimiter survived segmentation"
            );
        }
    }

    let report = json(&report);
    let segment = &report["segment"];
    assert_eq!(segment["source_records"], 40);
    assert_eq!(segment["fragments_emitted"], out.len());
    // Ten reads begin with a delimiter, so ten empty prefixes are discarded.
    assert_eq!(segment["discarded_empty"], 10);
    assert_eq!(
        segment["internal_fragments"], 10,
        "the two-delimiter middles"
    );
    assert_eq!(report["counts"]["records_in"], 40, "source records");
    assert_eq!(report["counts"]["records_out"], out.len(), "fragments");
}

#[test]
fn terminal_fragments_and_short_fragments_can_be_dropped() {
    let schema = Schema {
        paired: false,
        quality: true,
        headers: true,
        flags: false,
    };
    let fixture = Fixture::new();
    let records = segment_records(schema, 40);
    let input = fixture.input(schema, &records, BLOCK);
    let fasta = fixture.path("delimiters.fa");
    delimiter_fasta(&fasta);

    // Only fragments with a delimiter on both sides survive: the middles of the
    // two-delimiter reads, and nothing else.
    let output = fixture.path("internal.cbq");
    bqc_ok(&[
        "segment",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--adapter-fasta",
        fasta.to_str().unwrap(),
        "--terminal-fragments",
        "discard",
    ]);
    let (_, internal) = read_cbq(&output);
    assert_eq!(internal.len(), 10, "one middle per two-delimiter read");
    for fragment in &internal {
        assert_eq!(fragment.s_seq.len(), 20, "the planted middle length");
    }

    // A minimum length above the middles removes them too, and the report
    // accounts for them separately from the empty ones.
    let output = fixture.path("long.cbq");
    let report = fixture.path("report.json");
    bqc_ok(&[
        "segment",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--adapter-fasta",
        fasta.to_str().unwrap(),
        "--min-segment-length",
        "21",
        "--report",
        report.to_str().unwrap(),
    ]);
    let (_, long) = read_cbq(&output);
    assert!(
        long.iter().all(|fragment| fragment.s_seq.len() >= 21),
        "no fragment shorter than the minimum survives"
    );
    let report = json(&report);
    assert_eq!(
        report["segment"]["discarded_too_short"], 20,
        "ten middles and ten 20 base suffixes"
    );

    // The safety limit truncates a read to its first fragments.
    let output = fixture.path("capped.cbq");
    let report = fixture.path("capped.json");
    bqc_ok(&[
        "segment",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--adapter-fasta",
        fasta.to_str().unwrap(),
        "--max-segments-per-read",
        "1",
        "--report",
        report.to_str().unwrap(),
    ]);
    let (_, capped) = read_cbq(&output);
    assert_eq!(capped.len(), records.len(), "exactly one fragment per read");
    let report = json(&report);
    assert_eq!(report["segment"]["discarded_over_limit"], 30);
}

#[test]
fn fragments_are_trimmed_and_filtered_individually() {
    let schema = Schema {
        paired: false,
        quality: true,
        headers: true,
        flags: false,
    };
    let fixture = Fixture::new();
    let records = segment_records(schema, 40);
    let input = fixture.input(schema, &records, BLOCK);
    let fasta = fixture.path("delimiters.fa");
    let output = fixture.path("fragments.cbq");
    let failed = fixture.path("failed.cbq");
    let reasons = fixture.path("reasons.tsv");
    delimiter_fasta(&fasta);

    // Five bases off the front of every fragment, then a length filter that only
    // the longest fragments survive.
    bqc_ok(&[
        "segment",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--adapter-fasta",
        fasta.to_str().unwrap(),
        "--front",
        "5",
        "--min-length",
        "25",
        "--failed",
        failed.to_str().unwrap(),
        "--failed-reasons",
        reasons.to_str().unwrap(),
        "-T",
        "3",
    ]);

    let (_, accepted) = read_cbq(&output);
    let (_, rejected) = read_cbq(&failed);
    assert!(!accepted.is_empty() && !rejected.is_empty());
    for fragment in &accepted {
        assert!(
            fragment.s_seq.len() >= 25,
            "the filter admitted a short one"
        );
        // Trimming is visible in the sequence, but the header still describes the
        // fragment's position in the source read.
        let header = String::from_utf8(fragment.s_header.clone().unwrap()).unwrap();
        assert!(header.contains("|segment="), "{header}");
    }
    // Every rejected fragment is explained, one row per fragment.
    let text = text(&reasons);
    let rows: Vec<&str> = text.lines().skip(1).collect();
    assert_eq!(rows.len(), rejected.len());
    for row in rows {
        let columns: Vec<&str> = row.split('\t').collect();
        assert!(columns[1].starts_with("segment"), "{row}");
        assert_eq!(columns[2], "FAIL");
        assert!(columns[3].contains("TOO_SHORT"), "{row}");
    }
}

#[test]
fn segment_output_is_identical_at_any_thread_count() {
    let schema = Schema {
        paired: false,
        quality: true,
        headers: true,
        flags: true,
    };
    let fixture = Fixture::new();
    // Several blocks, so the committer has to reorder chunks.
    let records = segment_records(schema, 400);
    let input = fixture.input(schema, &records, 4096);
    let fasta = fixture.path("delimiters.fa");
    delimiter_fasta(&fasta);
    assert!(count_blocks(&input) > 1, "the fixture must span blocks");

    let mut reference: Option<(Vec<u8>, String)> = None;
    for threads in ["1", "2", "5"] {
        let output = fixture.path(&format!("fragments_{threads}.cbq"));
        let sidecar = fixture.path(&format!("segments_{threads}.tsv"));
        bqc_ok(&[
            "segment",
            input.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
            "--adapter-fasta",
            fasta.to_str().unwrap(),
            "--segments",
            sidecar.to_str().unwrap(),
            "-T",
            threads,
        ]);
        let bytes = std::fs::read(&output).expect("read output");
        let rows = text(&sidecar);
        match &reference {
            None => reference = Some((bytes, rows)),
            Some((expected_bytes, expected_rows)) => {
                assert_eq!(&bytes, expected_bytes, "output differs at -T {threads}");
                assert_eq!(&rows, expected_rows, "sidecar differs at -T {threads}");
            }
        }
    }
}

#[test]
fn segmentation_requires_a_single_end_input_and_a_delimiter() {
    let schema = Schema {
        paired: true,
        quality: true,
        headers: true,
        flags: false,
    };
    let fixture = Fixture::new();
    let records = segment_records(schema, 8);
    let input = fixture.input(schema, &records, BLOCK);
    let output = fixture.path("fragments.cbq");

    let error = bqc(&[
        "segment",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--adapter-r1",
        DELIMITER_A,
    ])
    .unwrap_err();
    assert!(
        format!("{error}").contains("requires a single-end input file"),
        "{error}"
    );
    assert!(!output.exists(), "nothing is written for a rejected run");

    // And a delimiter is not optional: there is nothing to split on without one.
    let single = Schema {
        paired: false,
        ..schema
    };
    let fixture = Fixture::new();
    let input = fixture.input(single, &segment_records(single, 8), BLOCK);
    let output = fixture.path("fragments.cbq");
    let error = bqc(&[
        "segment",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
    ])
    .unwrap_err();
    assert!(
        format!("{error}").contains("requires --adapter-r1 or --adapter-fasta"),
        "{error}"
    );
}

#[test]
fn a_failed_segment_run_leaves_no_output_behind() {
    let schema = Schema {
        paired: false,
        quality: true,
        headers: true,
        flags: false,
    };
    let fixture = Fixture::new();
    let records = segment_records(schema, 12);
    let input = fixture.input(schema, &records, BLOCK);
    let output = fixture.path("fragments.cbq");
    // The sidecar destination is a directory, so creating it fails after the
    // output file has already been opened.
    let sidecar = fixture.path("occupied");
    std::fs::create_dir(&sidecar).expect("create directory");

    let error = bqc(&[
        "segment",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--adapter-r1",
        DELIMITER_A,
        "--segments",
        sidecar.to_str().unwrap(),
    ])
    .unwrap_err();
    assert!(format!("{error}").contains("occupied"), "{error}");
    assert!(!output.exists(), "no output is renamed into place");
    // And no temporary file survives either.
    let leftovers: Vec<String> = std::fs::read_dir(fixture.dir.path())
        .expect("list directory")
        .map(|entry| entry.unwrap().file_name().to_string_lossy().into_owned())
        .filter(|name| name.contains("bqc-tmp"))
        .collect();
    assert!(leftovers.is_empty(), "temporary files left: {leftovers:?}");
}

#[test]
fn indel_aware_segmentation_cuts_at_a_damaged_delimiter() {
    let schema = Schema {
        paired: false,
        quality: true,
        headers: true,
        flags: false,
    };
    let fixture = Fixture::new();
    let mut sequences = Sequences::new(0x1_1DE1);
    // One delimiter copy per read, missing its fifth base.
    let mut damaged = DELIMITER_A.as_bytes().to_vec();
    damaged.remove(4);
    let records: Vec<Record> = (0..16)
        .map(|index| {
            let left = sequences.sequence(25);
            let right = sequences.sequence(25);
            let read = [left.as_slice(), &damaged, &right].concat();
            let quality = sequences.quality(read.len());
            Record::new(schema, index, &read, &read).with_quality(&quality, &quality)
        })
        .collect();
    let input = fixture.input(schema, &records, BLOCK);
    let output = fixture.path("fragments.cbq");

    bqc_ok(&[
        "segment",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--adapter-r1",
        DELIMITER_A,
        "--allow-indels",
    ]);

    let (_, out) = read_cbq(&output);
    assert_eq!(out.len(), records.len() * 2, "every read is split in two");
    for (index, record) in records.iter().enumerate() {
        // The cut uses the read bases the alignment consumed, so no base of the
        // damaged delimiter survives at either fragment boundary.
        assert_eq!(out[index * 2].s_seq, record.s_seq[..25], "prefix {index}");
        assert_eq!(
            out[index * 2 + 1].s_seq,
            record.s_seq[25 + damaged.len()..],
            "suffix {index}"
        );
    }
}

#[test]
fn a_header_free_input_keeps_its_schema_and_requires_the_sidecar() {
    let schema = Schema {
        paired: false,
        quality: true,
        headers: false,
        flags: false,
    };
    let fixture = Fixture::new();
    let records = segment_records(schema, 12);
    let input = fixture.input(schema, &records, BLOCK);
    let output = fixture.path("fragments.cbq");

    // Without headers there is nowhere to put the `|segment=` suffix, so the
    // sidecar is the only provenance and is therefore mandatory.
    let error = bqc(&[
        "segment",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--adapter-r1",
        DELIMITER_A,
    ])
    .unwrap_err();
    assert!(
        format!("{error}").contains("--segments is required"),
        "{error}"
    );

    let sidecar = fixture.path("segments.tsv");
    bqc_ok(&[
        "segment",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--adapter-r1",
        DELIMITER_A,
        "--segments",
        sidecar.to_str().unwrap(),
    ]);
    let (out_schema, out) = read_cbq(&output);
    assert_eq!(out_schema, schema, "a header-free input stays header-free");
    assert!(out.iter().all(|fragment| fragment.s_header.is_none()));
    let sidecar = text(&sidecar);
    let rows: Vec<&str> = sidecar.lines().skip(1).collect();
    assert_eq!(rows.len(), out.len(), "every fragment is in the sidecar");
    for row in rows {
        let columns: Vec<&str> = row.split('\t').collect();
        assert_eq!(columns[8], ".", "there is no original header to report");
    }
}

// ---------------------------------------------------------------------- sniff

/// Records whose contamination sits only in the final `tail` of the file.
///
/// A sampler that reads the leading records sees a clean library here; one that
/// spreads its sample across the file sees the adapter.
fn contaminated_tail(schema: Schema, count: usize, tail: usize) -> Vec<Record> {
    let mut sequences = Sequences::new(0x5C1F_F1A5);
    (0..count)
        .map(|index| {
            let insert = sequences.sequence(30 + index % 20);
            let clean_r1 = sequences.sequence(60);
            let clean_r2 = sequences.sequence(60);
            let (r1, r2) = if index >= count - tail {
                (
                    [insert.as_slice(), ADAPTER_R1.as_bytes()].concat(),
                    [insert.as_slice(), ADAPTER_R2.as_bytes()].concat(),
                )
            } else {
                (clean_r1, clean_r2)
            };
            let q1 = vec![b'I'; r1.len()];
            let q2 = vec![b'I'; r2.len()];
            Record::new(schema, index, &r1, &r2).with_quality(&q1, &q2)
        })
        .collect()
}

/// Pairs where `every`th record carries the full adapter on both mates.
fn sniff_records(schema: Schema, count: usize, every: usize) -> Vec<Record> {
    let mut sequences = Sequences::new(0xF00D_5EED);
    (0..count)
        .map(|index| {
            let insert = sequences.sequence(30 + index % 20);
            let (r1, r2) = if index % every == 0 {
                (
                    [insert.as_slice(), ADAPTER_R1.as_bytes()].concat(),
                    [insert.as_slice(), ADAPTER_R2.as_bytes()].concat(),
                )
            } else {
                (sequences.sequence(80), sequences.sequence(80))
            };
            let q1 = vec![b'I'; r1.len()];
            let q2 = vec![b'I'; r2.len()];
            Record::new(schema, index, &r1, &r2).with_quality(&q1, &q2)
        })
        .collect()
}

const SINGLE: Schema = Schema {
    paired: false,
    quality: true,
    headers: true,
    flags: false,
};

const PAIRED: Schema = Schema {
    paired: true,
    quality: true,
    headers: true,
    flags: false,
};

#[test]
fn sniff_adapters_identifies_the_library_in_both_mates() {
    let fixture = Fixture::new();
    let records = sniff_records(PAIRED, 4000, 2);
    let input = fixture.input(PAIRED, &records, BLOCK);
    let out = fixture.path("adapters.json");

    bqc_ok(&[
        "sniff",
        "adapters",
        input.to_str().unwrap(),
        "--format",
        "json",
        "-o",
        out.to_str().unwrap(),
    ]);

    let report = json(&out);
    assert_eq!(report["schema_version"], 1);
    assert_eq!(report["command"], "sniff adapters");
    assert_eq!(report["sample"]["method"], "deterministic-distributed");
    assert_eq!(report["result"]["r1"]["decision"], "confident");
    assert_eq!(report["result"]["r2"]["decision"], "confident");
    assert_eq!(
        report["result"]["r1"]["recommended_sequence"]
            .as_str()
            .unwrap(),
        ADAPTER_R1
    );
    assert_eq!(
        report["result"]["r2"]["recommended_sequence"]
            .as_str()
            .unwrap(),
        ADAPTER_R2
    );

    let leader = &report["result"]["r1"]["candidates"][0];
    assert_eq!(leader["confidence"], "high");
    assert_eq!(leader["known_name"], "illumina-truseq");
    assert_eq!(leader["known_category"], "adapter");
    assert_eq!(leader["evidence_sources"][0], "known_database");
    // Half the records are contaminated, and support is measured once.
    let support = leader["support_fraction"].as_f64().unwrap();
    assert!((0.4..=0.6).contains(&support), "support fraction {support}");
    assert!(leader["supporting_reads"].as_u64().unwrap() > 100);
}

#[test]
fn an_internal_known_adapter_is_not_recommended() {
    // Regression: identity with a known adapter is not enough to justify
    // trimming. Every read carries the exact sequence internally, surrounded by
    // independently generated biological sequence on both sides.
    let fixture = Fixture::new();
    let mut sequences = Sequences::new(0x1A7E_4A11);
    let records: Vec<Record> = (0..4000)
        .map(|index| {
            let prefix = sequences.sequence(120);
            let suffix = sequences.sequence(120);
            let read = [prefix.as_slice(), ADAPTER_R1.as_bytes(), suffix.as_slice()].concat();
            let quality = vec![b'I'; read.len()];
            Record::new(SINGLE, index, &read, b"").with_quality(&quality, b"")
        })
        .collect();
    let input = fixture.input(SINGLE, &records, BLOCK);
    let out = fixture.path("internal.json");

    bqc_ok(&[
        "sniff",
        "adapters",
        input.to_str().unwrap(),
        "--format",
        "json",
        "-o",
        out.to_str().unwrap(),
    ]);

    let report = json(&out);
    let result = &report["result"]["r1"];
    assert_eq!(result["decision"], "inconclusive", "{report}");
    assert!(result["recommended_sequence"].is_null());
    let known = result["candidates"]
        .as_array()
        .unwrap()
        .iter()
        .find(|candidate| candidate["known_name"] == "illumina-truseq")
        .expect("known internal candidate is still reported");
    assert_eq!(known["confidence"], "medium", "{known}");
    assert_eq!(known["median_start"], 120);
}

#[test]
fn a_start_only_known_primer_is_reported_but_never_auto_selected() {
    const PRIMER: &[u8] = b"AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA";

    // : sequence identity at the 5' boundary identifies this primer,
    // but it does not establish a safe 3' trimming coordinate.
    let fixture = Fixture::new();
    let mut sequences = Sequences::new(0x5A17_0A11);
    let records: Vec<Record> = (0..4000)
        .map(|index| {
            let biological = sequences.sequence(160);
            let read = [PRIMER, biological.as_slice()].concat();
            let quality = vec![b'I'; read.len()];
            Record::new(SINGLE, index, &read, b"").with_quality(&quality, b"")
        })
        .collect();
    let input = fixture.input(SINGLE, &records, BLOCK);
    let sniffed = fixture.path("start-only.json");

    bqc_ok(&[
        "sniff",
        "adapters",
        input.to_str().unwrap(),
        "--format",
        "json",
        "-o",
        sniffed.to_str().unwrap(),
    ]);

    let report = json(&sniffed);
    let result = &report["result"]["r1"];
    assert_eq!(result["decision"], "inconclusive", "{report}");
    assert!(result["recommended_sequence"].is_null());
    let primer = result["candidates"]
        .as_array()
        .unwrap()
        .iter()
        .find(|candidate| candidate["known_name"] == "illumina-expression-pcr-primer-2")
        .expect("the start-only primer is still reported");
    assert_eq!(primer["known_category"], "primer");
    assert_eq!(primer["confidence"], "medium", "{primer}");
    assert_eq!(primer["median_start"], 0);

    let output = fixture.path("auto.cbq");
    let report_path = fixture.path("auto-report.json");
    bqc_ok(&[
        "adapter",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--auto-detect",
        "--report",
        report_path.to_str().unwrap(),
    ]);
    assert_eq!(read_cbq(&output).1.len(), 4000);
    let report: serde_json::Value = serde_json::from_str(&text(&report_path)).unwrap();
    assert_eq!(
        report["adapter"]["detection"]["r1"]["decision"],
        "inconclusive"
    );
    assert_eq!(report["adapter"]["r1_reads_trimmed"], 0);
}

#[test]
fn a_terminal_minimum_overlap_recovers_the_full_known_adapter() {
    let fixture = Fixture::new();
    let mut sequences = Sequences::new(0x8BA5_EED5);
    let records: Vec<Record> = (0..400)
        .map(|index| {
            let insert = sequences.sequence(70 + index % 20);
            let read = [insert.as_slice(), &ADAPTER_R1.as_bytes()[..8]].concat();
            let quality = vec![b'I'; read.len()];
            Record::new(SINGLE, index, &read, b"").with_quality(&quality, b"")
        })
        .collect();
    let input = fixture.input(SINGLE, &records, BLOCK);
    let out = fixture.path("partial.json");

    bqc_ok(&[
        "sniff",
        "adapters",
        input.to_str().unwrap(),
        "--format",
        "json",
        "-o",
        out.to_str().unwrap(),
    ]);

    let report = json(&out);
    let result = &report["result"]["r1"];
    assert_eq!(result["decision"], "confident", "{report}");
    assert_eq!(result["recommended_sequence"], ADAPTER_R1);
    let leader = &result["candidates"][0];
    assert_eq!(leader["known_category"], "adapter");
    assert!(
        leader["tail_enrichment"].is_null(),
        "tail-only enrichment has explicit null semantics: {leader}"
    );
}

#[test]
fn sniff_reads_the_whole_file_not_just_the_first_records() {
    // Contamination lives only in the last 20% of a multi-block file. A
    // leading-prefix sampler cannot see it; the distributed sample must.
    let fixture = Fixture::new();
    let records = contaminated_tail(PAIRED, 5000, 1000);
    let input = fixture.input(PAIRED, &records, 1 << 14);
    assert!(count_blocks(&input) > 1, "fixture must span several blocks");
    let out = fixture.path("adapters.json");

    bqc_ok(&[
        "sniff",
        "adapters",
        input.to_str().unwrap(),
        "--sample-size",
        "500",
        "--format",
        "json",
        "-o",
        out.to_str().unwrap(),
    ]);

    let report = json(&out);
    let leader = &report["result"]["r1"]["candidates"][0];
    assert_eq!(leader["known_name"], "illumina-truseq", "{report}");
    assert!(
        leader["supporting_reads"].as_u64().unwrap() > 50,
        "the tail contamination was missed: {leader}"
    );
    assert_eq!(report["sample"]["selected"], 500);
}

#[test]
fn sniff_never_modifies_its_input() {
    let fixture = Fixture::new();
    let records = sniff_records(PAIRED, 2000, 3);
    let input = fixture.input(PAIRED, &records, BLOCK);
    let before = std::fs::read(&input).unwrap();

    bqc_ok(&[
        "sniff",
        "adapters",
        input.to_str().unwrap(),
        "--format",
        "json",
    ]);

    let after = std::fs::read(&input).unwrap();
    assert_eq!(before, after, "sniffing rewrote its input");
    // And no stray output was created next to it.
    let (_, unchanged) = read_cbq(&input);
    assert_eq!(unchanged.len(), records.len());
}

#[test]
fn sniff_adapters_is_identical_at_any_thread_count() {
    let fixture = Fixture::new();
    let records = sniff_records(PAIRED, 6000, 2);
    let input = fixture.input(PAIRED, &records, 1 << 14);
    assert!(count_blocks(&input) > 2, "fixture must span several blocks");

    let mut rendered = Vec::new();
    for threads in ["1", "3", "8"] {
        let out = fixture.path(&format!("t{threads}.json"));
        bqc_ok(&[
            "sniff",
            "adapters",
            input.to_str().unwrap(),
            "--sample-size",
            "1500",
            "--format",
            "json",
            "-T",
            threads,
            "-o",
            out.to_str().unwrap(),
        ]);
        rendered.push(text(&out));
    }
    assert_eq!(rendered[0], rendered[1], "-T 1 and -T 3 disagree");
    assert_eq!(rendered[0], rendered[2], "-T 1 and -T 8 disagree");
}

#[test]
fn a_clean_library_is_inconclusive_and_that_is_not_a_failure() {
    let fixture = Fixture::new();
    let records = simple_records(PAIRED, 2000, 100);
    let input = fixture.input(PAIRED, &records, BLOCK);
    let out = fixture.path("adapters.json");

    // Without --require-confident an inconclusive answer still succeeds.
    bqc_ok(&[
        "sniff",
        "adapters",
        input.to_str().unwrap(),
        "--format",
        "json",
        "-o",
        out.to_str().unwrap(),
    ]);
    let report = json(&out);
    assert_eq!(report["result"]["r1"]["decision"], "inconclusive");
    assert!(report["result"]["r1"]["recommended_sequence"].is_null());

    // With it, the run reports a distinct outcome rather than an error.
    let outcome = common::bqc_outcome(&[
        "sniff",
        "adapters",
        input.to_str().unwrap(),
        "--require-confident",
        "--format",
        "json",
        "-o",
        fixture.path("strict.json").to_str().unwrap(),
    ])
    .expect("an inconclusive result is not an error");
    assert_eq!(outcome, bqc::cli::Outcome::NotConfident);
}

#[test]
fn a_confident_result_satisfies_require_confident() {
    let fixture = Fixture::new();
    let records = sniff_records(PAIRED, 4000, 2);
    let input = fixture.input(PAIRED, &records, BLOCK);
    let outcome = common::bqc_outcome(&[
        "sniff",
        "adapters",
        input.to_str().unwrap(),
        "--require-confident",
        "--format",
        "json",
        "-o",
        fixture.path("out.json").to_str().unwrap(),
    ])
    .expect("run succeeds");
    assert_eq!(outcome, bqc::cli::Outcome::Success);
}

#[test]
fn every_projection_describes_the_same_result() {
    let fixture = Fixture::new();
    let records = sniff_records(PAIRED, 3000, 2);
    let input = fixture.input(PAIRED, &records, BLOCK);

    for format in ["text", "json", "tsv"] {
        let out = fixture.path(&format!("out.{format}"));
        bqc_ok(&[
            "sniff",
            "adapters",
            input.to_str().unwrap(),
            "--format",
            format,
            "-o",
            out.to_str().unwrap(),
        ]);
        let rendered = text(&out);
        assert!(
            rendered.contains(ADAPTER_R1),
            "{format} output omits the adapter: {rendered}"
        );
    }

    // TSV is one header plus one row per reported candidate, both mates.
    let tsv = text(&fixture.path("out.tsv"));
    let mut lines = tsv.lines();
    let header = lines.next().unwrap();
    assert!(header.starts_with("input\tmate\tdecision\tsequence\t"));
    let rows: Vec<&str> = lines.collect();
    assert!(rows.iter().any(|row| row.contains("\tR1\tconfident\t")));
    assert!(rows.iter().any(|row| row.contains("\tR2\tconfident\t")));

    // Text names the decision in words.
    let rendered = text(&fixture.path("out.text"));
    assert!(
        rendered.contains("R1 adapter result: confident"),
        "{rendered}"
    );
    assert!(rendered.contains("Sampling: deterministic-distributed"));
}

#[test]
fn a_single_end_file_reports_one_mate() {
    let schema = Schema {
        paired: false,
        quality: true,
        headers: true,
        flags: false,
    };
    let fixture = Fixture::new();
    let records = sniff_records(schema, 3000, 2);
    let input = fixture.input(schema, &records, BLOCK);
    let out = fixture.path("out.json");

    bqc_ok(&[
        "sniff",
        "adapters",
        input.to_str().unwrap(),
        "--format",
        "json",
        "-o",
        out.to_str().unwrap(),
    ]);
    let report = json(&out);
    assert_eq!(report["result"]["r1"]["decision"], "confident");
    assert!(report["result"]["r2"].is_null(), "single-end has no R2");
}

#[test]
fn sniffing_works_without_qualities_or_headers() {
    // Adapter discovery reads sequence only, so the leanest schema must work.
    let schema = Schema {
        paired: true,
        quality: false,
        headers: false,
        flags: false,
    };
    let fixture = Fixture::new();
    let records = sniff_records(schema, 3000, 2);
    let input = fixture.input(schema, &records, BLOCK);
    let out = fixture.path("out.json");

    bqc_ok(&[
        "sniff",
        "adapters",
        input.to_str().unwrap(),
        "--format",
        "json",
        "-o",
        out.to_str().unwrap(),
    ]);
    let report = json(&out);
    assert_eq!(report["input"]["quality"], false);
    assert_eq!(report["input"]["headers"], false);
    assert_eq!(report["result"]["r1"]["decision"], "confident");
}

#[test]
fn a_span_restricts_sniffing_to_those_records() {
    let fixture = Fixture::new();
    // Only the first 1000 records are contaminated.
    let records = contaminated_tail(PAIRED, 4000, 4000)
        .into_iter()
        .take(1000)
        .chain(simple_records(PAIRED, 3000, 100))
        .collect::<Vec<_>>();
    let input = fixture.input(PAIRED, &records, 1 << 14);

    let clean = fixture.path("clean.json");
    bqc_ok(&[
        "sniff",
        "adapters",
        input.to_str().unwrap(),
        "--span",
        "1000..4000",
        "--format",
        "json",
        "-o",
        clean.to_str().unwrap(),
    ]);
    let report = json(&clean);
    assert_eq!(report["sample"]["range_start"], 1000);
    assert_eq!(report["sample"]["range_end"], 4000);
    assert_eq!(report["result"]["r1"]["decision"], "inconclusive");

    let dirty = fixture.path("dirty.json");
    bqc_ok(&[
        "sniff",
        "adapters",
        input.to_str().unwrap(),
        "--span",
        "0..1000",
        "--format",
        "json",
        "-o",
        dirty.to_str().unwrap(),
    ]);
    assert_eq!(json(&dirty)["result"]["r1"]["decision"], "confident");
}

/// An adapter absent from the bundled library, planted after a variable insert.
const NOVEL_ADAPTER: &[u8] = b"GTCGATCGTACGGCATCCGATCGTACGATCGGCATTAGCGCATTAGCGG";

#[test]
fn de_novo_discovery_finds_an_adapter_that_is_not_in_the_library() {
    let fixture = Fixture::new();
    let mut sequences = Sequences::new(0xDE_1204);
    let records: Vec<Record> = (0..6000)
        .map(|index| {
            let insert = sequences.sequence(20 + index % 25);
            let (r1, r2) = if index % 2 == 0 {
                (
                    [insert.as_slice(), NOVEL_ADAPTER].concat(),
                    [insert.as_slice(), NOVEL_ADAPTER].concat(),
                )
            } else {
                (sequences.sequence(90), sequences.sequence(90))
            };
            let q1 = vec![b'I'; r1.len()];
            let q2 = vec![b'I'; r2.len()];
            Record::new(PAIRED, index, &r1, &r2).with_quality(&q1, &q2)
        })
        .collect();
    let input = fixture.input(PAIRED, &records, BLOCK);
    let out = fixture.path("out.json");

    bqc_ok(&[
        "sniff",
        "adapters",
        input.to_str().unwrap(),
        "--format",
        "json",
        "-o",
        out.to_str().unwrap(),
    ]);

    let report = json(&out);
    let leader = &report["result"]["r1"]["candidates"][0];
    assert_eq!(report["result"]["r1"]["decision"], "confident", "{report}");
    assert_eq!(leader["confidence"], "high");
    // Assembled from k-mer evidence alone: nothing in the library matches.
    assert_eq!(leader["evidence_sources"][0], "kmer_consensus");
    assert!(leader["known_name"].is_null());

    // The consensus must be a genuine stretch of the planted adapter, not a
    // seed-length fragment.
    let found = leader["sequence"].as_str().unwrap().as_bytes();
    assert!(found.len() >= 20, "consensus is only {} bases", found.len());
    let novel = String::from_utf8_lossy(NOVEL_ADAPTER);
    assert!(
        novel.contains(leader["sequence"].as_str().unwrap()),
        "assembled {} which is not part of the planted adapter",
        leader["sequence"]
    );
}

#[test]
fn known_and_de_novo_evidence_for_one_adapter_become_one_candidate() {
    let fixture = Fixture::new();
    let records = sniff_records(PAIRED, 6000, 2);
    let input = fixture.input(PAIRED, &records, BLOCK);
    let out = fixture.path("out.json");

    bqc_ok(&[
        "sniff",
        "adapters",
        input.to_str().unwrap(),
        "--format",
        "json",
        "-o",
        out.to_str().unwrap(),
    ]);

    let report = json(&out);
    let candidates = report["result"]["r1"]["candidates"].as_array().unwrap();
    let leader = &candidates[0];
    // One candidate carrying both provenances, not two candidates.
    let sources: Vec<&str> = leader["evidence_sources"]
        .as_array()
        .unwrap()
        .iter()
        .map(|s| s.as_str().unwrap())
        .collect();
    assert_eq!(
        sources,
        vec!["known_database", "kmer_consensus"],
        "{leader}"
    );
    assert_eq!(leader["known_name"], "illumina-truseq");

    // Support is measured once: it cannot exceed the reads that were sampled.
    let sampled = report["result"]["r1"]["sampled_reads"].as_u64().unwrap();
    for candidate in candidates {
        assert!(
            candidate["supporting_reads"].as_u64().unwrap() <= sampled,
            "support was summed across evidence sources: {candidate}"
        );
    }
}

#[test]
fn an_abundant_repeat_is_not_reported_as_an_adapter() {
    // A motif at the read's 3' end in every read, but followed by different
    // sequence each time — a biological repeat, not read-through. It is
    // enriched, so it survives k-mer counting; only the consensus extension
    // distinguishes it.
    let fixture = Fixture::new();
    let mut sequences = Sequences::new(0x8EEF_00D1);
    let records: Vec<Record> = (0..6000)
        .map(|index| {
            let mut r1 = sequences.sequence(40);
            r1.extend_from_slice(b"CCTGAGCTAAGCTT");
            r1.extend(sequences.sequence(30));
            let mut r2 = sequences.sequence(40);
            r2.extend_from_slice(b"CCTGAGCTAAGCTT");
            r2.extend(sequences.sequence(30));
            let q1 = vec![b'I'; r1.len()];
            let q2 = vec![b'I'; r2.len()];
            Record::new(PAIRED, index, &r1, &r2).with_quality(&q1, &q2)
        })
        .collect();
    let input = fixture.input(PAIRED, &records, BLOCK);
    let out = fixture.path("out.json");

    bqc_ok(&[
        "sniff",
        "adapters",
        input.to_str().unwrap(),
        "--format",
        "json",
        "-o",
        out.to_str().unwrap(),
    ]);

    let report = json(&out);
    assert_eq!(
        report["result"]["r1"]["decision"], "inconclusive",
        "an abundant repeat was recommended as an adapter: {report}"
    );
    assert!(report["result"]["r1"]["recommended_sequence"].is_null());
}

#[test]
fn a_poly_a_tail_is_reported_as_an_artifact_not_an_adapter() {
    let fixture = Fixture::new();
    let mut sequences = Sequences::new(0xA11A_AA01);
    let records: Vec<Record> = (0..4000)
        .map(|index| {
            let mut r1 = sequences.sequence(40);
            r1.extend(std::iter::repeat_n(b'A', 70));
            let mut r2 = sequences.sequence(40);
            r2.extend(std::iter::repeat_n(b'A', 70));
            let q1 = vec![b'I'; r1.len()];
            let q2 = vec![b'I'; r2.len()];
            Record::new(PAIRED, index, &r1, &r2).with_quality(&q1, &q2)
        })
        .collect();
    let input = fixture.input(PAIRED, &records, BLOCK);
    let out = fixture.path("out.json");

    bqc_ok(&[
        "sniff",
        "adapters",
        input.to_str().unwrap(),
        "--format",
        "json",
        "-o",
        out.to_str().unwrap(),
    ]);

    let report = json(&out);
    let mate = &report["result"]["r1"];
    assert_eq!(mate["decision"], "inconclusive", "{report}");
    assert_eq!(mate["poly_a_signal"], mate["sampled_reads"]);
    for candidate in mate["candidates"].as_array().unwrap() {
        let sequence = candidate["sequence"].as_str().unwrap();
        assert!(
            !sequence.contains("AAAAAAAAAA"),
            "a poly-A run was proposed as an adapter: {sequence}"
        );
    }
}

#[test]
fn de_novo_discovery_is_identical_at_any_thread_count() {
    let fixture = Fixture::new();
    let mut sequences = Sequences::new(0x7EAD_1234);
    let records: Vec<Record> = (0..8000)
        .map(|index| {
            let insert = sequences.sequence(20 + index % 25);
            let (r1, r2) = if index % 2 == 0 {
                (
                    [insert.as_slice(), NOVEL_ADAPTER].concat(),
                    [insert.as_slice(), NOVEL_ADAPTER].concat(),
                )
            } else {
                (sequences.sequence(90), sequences.sequence(90))
            };
            let q1 = vec![b'I'; r1.len()];
            let q2 = vec![b'I'; r2.len()];
            Record::new(PAIRED, index, &r1, &r2).with_quality(&q1, &q2)
        })
        .collect();
    let input = fixture.input(PAIRED, &records, 1 << 14);
    assert!(count_blocks(&input) > 2, "fixture must span several blocks");

    let mut rendered = Vec::new();
    for threads in ["1", "3", "8"] {
        let out = fixture.path(&format!("t{threads}.json"));
        bqc_ok(&[
            "sniff",
            "adapters",
            input.to_str().unwrap(),
            "--sample-size",
            "2000",
            "--format",
            "json",
            "-T",
            threads,
            "-o",
            out.to_str().unwrap(),
        ]);
        rendered.push(text(&out));
    }
    assert_eq!(rendered[0], rendered[1], "-T 1 and -T 3 disagree");
    assert_eq!(rendered[0], rendered[2], "-T 1 and -T 8 disagree");
}

/// A second novel adapter, so the two mates can carry different ones.
const NOVEL_ADAPTER_R2: &[u8] = b"TTGCAGCTAGCTGGATCCTTAGCAGGCTTAAGCCTGATCGATCGGTACC";

/// Pairs whose insert is shorter than the read, so both mates read through into
/// adapter and the inferred insert boundary exposes the overhang directly.
fn short_insert_pairs(count: usize) -> Vec<Record> {
    let mut sequences = Sequences::new(0x0FF5_E701);
    (0..count)
        .map(|index| {
            let insert: Vec<u8> = sequences
                .sequence(35 + index % 10)
                .into_iter()
                .map(|base| if base == b'N' { b'A' } else { base })
                .collect();
            let reversed: Vec<u8> = insert.iter().rev().map(|&b| complement(b)).collect();
            let r1 = [insert.as_slice(), NOVEL_ADAPTER].concat();
            let r2 = [reversed.as_slice(), NOVEL_ADAPTER_R2].concat();
            let q1 = vec![b'I'; r1.len()];
            let q2 = vec![b'I'; r2.len()];
            Record::new(PAIRED, index, &r1, &r2).with_quality(&q1, &q2)
        })
        .collect()
}

#[test]
fn overlap_overhangs_are_adapter_evidence_for_both_mates() {
    let fixture = Fixture::new();
    let records = short_insert_pairs(4000);
    let input = fixture.input(PAIRED, &records, BLOCK);
    let out = fixture.path("out.json");

    bqc_ok(&[
        "sniff",
        "adapters",
        input.to_str().unwrap(),
        "--format",
        "json",
        "-o",
        out.to_str().unwrap(),
    ]);

    let report = json(&out);
    let mut support = Vec::new();
    for (mate, planted) in [("r1", NOVEL_ADAPTER), ("r2", NOVEL_ADAPTER_R2)] {
        let result = &report["result"][mate];
        assert_eq!(result["decision"], "confident", "{mate}: {report}");
        let leader = &result["candidates"][0];
        let sources: Vec<&str> = leader["evidence_sources"]
            .as_array()
            .unwrap()
            .iter()
            .map(|s| s.as_str().unwrap())
            .collect();
        assert!(
            sources.contains(&"paired_overlap"),
            "{mate} has no overlap evidence: {sources:?}"
        );
        // The overhang past the inferred boundary is the adapter itself.
        let found = leader["sequence"].as_str().unwrap();
        assert!(
            String::from_utf8_lossy(planted).contains(found),
            "{mate} assembled {found}, which is not the planted adapter"
        );
        support.push(leader["paired_overlap_support"].as_u64().unwrap());
    }

    assert!(support[0] > 0, "no pair contributed overlap evidence");
    // One analysis per pair serves both mates, so the two counts are the same
    // number of pairs — not two independently inferred overlaps.
    assert_eq!(
        support[0], support[1],
        "the mates disagree on how many pairs overlapped, \
         so the overlap was not inferred once per pair"
    );
}

#[test]
fn a_single_end_file_has_no_overlap_evidence() {
    let schema = Schema {
        paired: false,
        quality: true,
        headers: true,
        flags: false,
    };
    let fixture = Fixture::new();
    let records = sniff_records(schema, 3000, 2);
    let input = fixture.input(schema, &records, BLOCK);
    let out = fixture.path("out.json");

    bqc_ok(&[
        "sniff",
        "adapters",
        input.to_str().unwrap(),
        "--format",
        "json",
        "-o",
        out.to_str().unwrap(),
    ]);

    let report = json(&out);
    for candidate in report["result"]["r1"]["candidates"].as_array().unwrap() {
        let sources: Vec<&str> = candidate["evidence_sources"]
            .as_array()
            .unwrap()
            .iter()
            .map(|s| s.as_str().unwrap())
            .collect();
        assert!(
            !sources.contains(&"paired_overlap"),
            "single-end input produced overlap evidence: {sources:?}"
        );
        assert_eq!(candidate["paired_overlap_support"], 0);
    }
}

#[test]
fn pairs_that_do_not_overlap_contribute_no_overhang() {
    // Inserts longer than the reads: there is no read-through and no overhang,
    // so the only evidence available is the library scan.
    let fixture = Fixture::new();
    let mut sequences = Sequences::new(0x100B_EE01);
    let records: Vec<Record> = (0..3000)
        .map(|index| {
            let r1 = sequences.sequence(100);
            let r2 = sequences.sequence(100);
            let q1 = vec![b'I'; r1.len()];
            let q2 = vec![b'I'; r2.len()];
            Record::new(PAIRED, index, &r1, &r2).with_quality(&q1, &q2)
        })
        .collect();
    let input = fixture.input(PAIRED, &records, BLOCK);
    let out = fixture.path("out.json");

    bqc_ok(&[
        "sniff",
        "adapters",
        input.to_str().unwrap(),
        "--format",
        "json",
        "-o",
        out.to_str().unwrap(),
    ]);

    let report = json(&out);
    assert_eq!(report["result"]["r1"]["decision"], "inconclusive");
    for candidate in report["result"]["r1"]["candidates"].as_array().unwrap() {
        assert_eq!(candidate["paired_overlap_support"], 0, "{candidate}");
    }
}

#[test]
fn a_confident_result_can_be_emitted_as_configuration() {
    let fixture = Fixture::new();
    let records = sniff_records(PAIRED, 4000, 2);
    let input = fixture.input(PAIRED, &records, BLOCK);
    let config = fixture.path("adapters.toml");

    bqc_ok(&[
        "sniff",
        "adapters",
        input.to_str().unwrap(),
        "--emit-config",
        config.to_str().unwrap(),
        "--format",
        "json",
        "-o",
        fixture.path("out.json").to_str().unwrap(),
    ]);

    let emitted = text(&config);
    assert!(emitted.contains("[adapter]"), "{emitted}");
    assert!(
        emitted.contains(&format!("r1 = \"{ADAPTER_R1}\"")),
        "{emitted}"
    );
    assert!(
        emitted.contains(&format!("r2 = \"{ADAPTER_R2}\"")),
        "{emitted}"
    );

    // The fragment is a valid bqc configuration: it drives a real run.
    let trimmed = fixture.path("trimmed.cbq");
    bqc_ok(&[
        "workflow",
        input.to_str().unwrap(),
        "-o",
        trimmed.to_str().unwrap(),
        "--config",
        config.to_str().unwrap(),
        "--no-trim",
        "--no-filter",
    ]);
    let (_, out_records) = read_cbq(&trimmed);
    assert_eq!(out_records.len(), records.len());
    for (index, record) in out_records.iter().enumerate() {
        if index % 2 == 0 {
            assert!(
                record.s_seq.len() < records[index].s_seq.len(),
                "record {index} kept its adapter"
            );
        }
    }
}

#[test]
fn configuration_is_never_emitted_without_a_confident_result() {
    let fixture = Fixture::new();
    let records = simple_records(PAIRED, 2000, 100);
    let input = fixture.input(PAIRED, &records, BLOCK);
    let config = fixture.path("adapters.toml");

    let error = bqc(&[
        "sniff",
        "adapters",
        input.to_str().unwrap(),
        "--emit-config",
        config.to_str().unwrap(),
        "--format",
        "json",
        "-o",
        fixture.path("out.json").to_str().unwrap(),
    ])
    .unwrap_err();
    assert!(format!("{error}").contains("uniquely confident"), "{error}");
    assert!(
        !config.exists(),
        "an unusable configuration was left behind"
    );
}

#[test]
fn sniff_report_and_configuration_must_have_distinct_destinations() {
    let fixture = Fixture::new();
    let records = sniff_records(PAIRED, 400, 2);
    let input = fixture.input(PAIRED, &records, BLOCK);
    let shared = fixture.path("shared.out");

    for force in [false, true] {
        let mut arguments = vec![
            "sniff",
            "adapters",
            input.to_str().unwrap(),
            "--format",
            "json",
            "-o",
            shared.to_str().unwrap(),
            "--emit-config",
            shared.to_str().unwrap(),
        ];
        if force {
            arguments.push("--force");
        }
        let error = bqc(&arguments).unwrap_err();
        assert!(
            format!("{error}").contains("more than one output"),
            "{error}"
        );
        assert!(!shared.exists(), "a colliding output was written");
    }
}

#[test]
fn sniffing_and_auto_detection_reach_the_same_conclusion() {
    // There is one detector. `sniff adapters` reports what it found and
    // `--auto-detect` consumes the same recommendation, so the two can never
    // disagree about what is in a file.
    let fixture = Fixture::new();
    let records = sniff_records(PAIRED, 4000, 2);
    let input = fixture.input(PAIRED, &records, BLOCK);

    let sniffed = fixture.path("sniff.json");
    bqc_ok(&[
        "sniff",
        "adapters",
        input.to_str().unwrap(),
        "--format",
        "json",
        "-o",
        sniffed.to_str().unwrap(),
    ]);

    let report = fixture.path("run.json");
    bqc_ok(&[
        "adapter",
        input.to_str().unwrap(),
        "-o",
        fixture.path("out.cbq").to_str().unwrap(),
        "--auto-detect",
        "--report",
        report.to_str().unwrap(),
    ]);

    let sniffed = json(&sniffed);
    let detected = json(&report);
    for mate in ["r1", "r2"] {
        assert_eq!(
            sniffed["result"][mate]["recommended_sequence"],
            detected["adapter"]["detection"][mate]["recommended_sequence"],
            "{mate}: the two paths disagree"
        );
        assert_eq!(
            sniffed["result"][mate]["decision"],
            detected["adapter"]["detection"][mate]["decision"]
        );
    }

    // And the recommendation is what the run actually configured.
    let configured = detected["configuration"]["workflow"]["adapter"]["r1"]
        .as_array()
        .unwrap();
    assert!(
        configured.iter().any(|adapter| {
            adapter["sequence"] == sniffed["result"]["r1"]["recommended_sequence"]
        }),
        "{configured:?}"
    );
}

#[test]
fn auto_detection_trims_with_the_strongest_of_two_libraries() {
    // Half the reads carry TruSeq, half carry an unrelated adapter. The
    // strongest candidate is trimmed and the report shows the mixed decision;
    // picking one silently would mistrim the other half, so the run surfaces
    // the choice instead of stopping.
    let fixture = Fixture::new();
    let mut sequences = Sequences::new(0x_ABCD_0001_u64);
    let records: Vec<Record> = (0..6000)
        .map(|index| {
            let insert = sequences.sequence(25 + index % 20);
            let adapter: &[u8] = if index % 2 == 0 {
                ADAPTER_R1.as_bytes()
            } else {
                NOVEL_ADAPTER
            };
            let r1 = [insert.as_slice(), adapter].concat();
            let q1 = vec![b'I'; r1.len()];
            Record::new(SINGLE, index, &r1, b"").with_quality(&q1, b"")
        })
        .collect();
    let input = fixture.input(SINGLE, &records, BLOCK);

    let sniffed = fixture.path("sniff.json");
    bqc_ok(&[
        "sniff",
        "adapters",
        input.to_str().unwrap(),
        "--format",
        "json",
        "-o",
        sniffed.to_str().unwrap(),
    ]);
    assert_eq!(json(&sniffed)["result"]["r1"]["decision"], "mixed");

    // The recommended candidate is exactly what `assemble` orders first.
    let sniffed_json = json(&sniffed);
    let winner = sniffed_json["result"]["r1"]["candidates"][0]["sequence"]
        .as_str()
        .unwrap();
    assert!(
        winner == ADAPTER_R1 || winner == std::str::from_utf8(NOVEL_ADAPTER).unwrap(),
        "{winner}"
    );

    let out = fixture.path("out.cbq");
    let report_path = fixture.path("report.json");
    bqc_ok(&[
        "adapter",
        input.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--auto-detect",
        "--report",
        report_path.to_str().unwrap(),
    ]);
    let report = json(&report_path);
    assert_eq!(
        report["configuration"]["workflow"]["adapter"]["r1"][0]["sequence"],
        winner
    );
    assert_eq!(report["adapter"]["detection"]["r1"]["decision"], "mixed");
    assert_eq!(
        report["adapter"]["detection"]["r1"]["recommended_sequence"],
        winner
    );
    let trimmed = report["adapter"]["r1_reads_trimmed"].as_u64().unwrap();
    assert!(trimmed > 0 && trimmed < 6000, "{trimmed}");
    assert!(
        report["adapter"]["r1_bases_removed"].as_u64().unwrap() > 0
    );
}
