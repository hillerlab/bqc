// Copyright (c) 2026 The Hiller Lab at the Senckenberg Gesellschaft für Naturforschung
// Distributed under the terms of the GNU General Public License, Version 3.0.

//! End-to-end tests for `bqc sniff strand`.
//!
//! These build a real Salmon index from a small synthetic transcriptome and map
//! real reads through it, so they exercise the actual mapper rather than a
//! stand-in. They only exist when the `sniff-strand` feature is enabled.

#![cfg(feature = "sniff-strand")]

mod common;

use std::path::{Path, PathBuf};

use bqc::io::Schema;
use common::{Record, bqc, bqc_outcome, write_cbq};
use tempfile::TempDir;

const BLOCK: usize = 1 << 20;

/// A deterministic transcript sequence with no long repeats.
fn transcript(seed: u64, length: usize) -> Vec<u8> {
    let mut state = seed | 1;
    (0..length)
        .map(|_| {
            state ^= state >> 12;
            state ^= state << 25;
            state ^= state >> 27;
            state = state.wrapping_mul(0x2545_F491_4F6C_DD1D);
            b"ACGT"[(state >> 33) as usize % 4]
        })
        .collect()
}

fn revcomp(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .rev()
        .map(|&base| match base {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            other => other,
        })
        .collect()
}

/// A small transcriptome and an index built from it.
struct Reference {
    dir: TempDir,
    transcripts: Vec<Vec<u8>>,
}

impl Reference {
    /// Builds a Salmon index in-process. No `salmon` binary is involved.
    fn build() -> Self {
        use std::fmt::Write as _;
        let dir = tempfile::tempdir().expect("tempdir");
        let transcripts: Vec<Vec<u8>> = (0..8).map(|i| transcript(0x7A01 + i * 2, 2000)).collect();

        let fasta = dir.path().join("transcripts.fa");
        let mut text = String::new();
        for (index, sequence) in transcripts.iter().enumerate() {
            let _ = writeln!(text, ">tx{index}");
            text.push_str(&String::from_utf8_lossy(sequence));
            text.push('\n');
        }
        std::fs::write(&fasta, text).expect("write transcriptome");

        let index_dir = dir.path().join("index");
        let mut options = salmon_index::IndexBuildOptions::new(vec![fasta], index_dir);
        options.threads = 2;
        salmon_index::build(&options).expect("build a salmon index");
        Self { dir, transcripts }
    }

    fn index(&self) -> PathBuf {
        self.dir.path().join("index")
    }

    fn transcriptome(&self) -> PathBuf {
        self.dir.path().join("transcripts.fa")
    }

    fn path(&self, name: &str) -> PathBuf {
        self.dir.path().join(name)
    }

    /// Pairs drawn from the transcripts in a known orientation.
    ///
    /// `first_strand` reverses the roles: R1 antisense and R2 sense, which is
    /// what a dUTP protocol produces and what Salmon calls `ISR`.
    fn pairs(&self, count: usize, first_strand: bool, schema: Schema) -> Vec<Record> {
        (0..count)
            .map(|index| {
                let transcript = &self.transcripts[index % self.transcripts.len()];
                let start = (index * 37) % (transcript.len() - 400);
                let sense = transcript[start..start + 120].to_vec();
                let downstream = &transcript[start + 200..start + 320];
                let antisense = revcomp(downstream);
                let (r1, r2) = if first_strand {
                    (revcomp(downstream), transcript[start..start + 120].to_vec())
                } else {
                    (sense, antisense)
                };
                let q1 = vec![b'I'; r1.len()];
                let q2 = vec![b'I'; r2.len()];
                Record::new(schema, index, &r1, &r2).with_quality(&q1, &q2)
            })
            .collect()
    }
}

const PAIRED: Schema = Schema {
    paired: true,
    quality: true,
    headers: true,
    flags: false,
};

fn json(path: &Path) -> serde_json::Value {
    let text = std::fs::read_to_string(path).expect("read report");
    serde_json::from_str(&text).expect("parse report")
}

#[test]
fn a_second_strand_library_is_reported_as_forward() {
    let reference = Reference::build();
    let records = reference.pairs(4000, false, PAIRED);
    let input = reference.path("reads.cbq");
    write_cbq(&input, PAIRED, &records, BLOCK);
    let out = reference.path("strand.json");

    bqc(&[
        "sniff",
        "strand",
        input.to_str().unwrap(),
        "--index",
        reference.index().to_str().unwrap(),
        "--min-informative",
        "100",
        "--format",
        "json",
        "-o",
        out.to_str().unwrap(),
    ])
    .expect("sniff strand runs");

    let report = json(&out);
    let result = &report["result"];
    assert_eq!(report["command"], "sniff strand");
    assert_eq!(result["strandedness"], "forward", "{result}");
    assert_eq!(result["salmon_library_type"], "ISF");
    assert_eq!(result["pair_orientation"], "inward");
    assert_eq!(result["featurecounts_strand"], 1);
    assert_eq!(result["htseq_stranded"], "yes");
    assert!(result["informative_records"].as_u64().unwrap() > 100);
}

#[test]
fn a_first_strand_library_is_reported_as_reverse() {
    let reference = Reference::build();
    let records = reference.pairs(4000, true, PAIRED);
    let input = reference.path("reads.cbq");
    write_cbq(&input, PAIRED, &records, BLOCK);
    let out = reference.path("strand.json");

    bqc(&[
        "sniff",
        "strand",
        input.to_str().unwrap(),
        "--index",
        reference.index().to_str().unwrap(),
        "--min-informative",
        "100",
        "--format",
        "json",
        "-o",
        out.to_str().unwrap(),
    ])
    .expect("sniff strand runs");

    let result = &json(&out)["result"];
    assert_eq!(result["strandedness"], "reverse", "{result}");
    assert_eq!(result["salmon_library_type"], "ISR");
    assert_eq!(result["featurecounts_strand"], 2);
    assert_eq!(result["htseq_stranded"], "reverse");
}

#[test]
fn strand_inference_does_not_depend_on_the_thread_count() {
    let reference = Reference::build();
    let records = reference.pairs(6000, true, PAIRED);
    let input = reference.path("reads.cbq");
    write_cbq(&input, PAIRED, &records, 1 << 14);

    let mut rendered = Vec::new();
    for threads in ["1", "3", "8"] {
        let out = reference.path(&format!("t{threads}.json"));
        bqc(&[
            "sniff",
            "strand",
            input.to_str().unwrap(),
            "--index",
            reference.index().to_str().unwrap(),
            "--min-informative",
            "100",
            "--format",
            "json",
            "-T",
            threads,
            "-o",
            out.to_str().unwrap(),
        ])
        .expect("sniff strand runs");
        rendered.push(std::fs::read_to_string(&out).expect("read report"));
    }
    assert_eq!(rendered[0], rendered[1], "-T 1 and -T 3 disagree");
    assert_eq!(rendered[0], rendered[2], "-T 1 and -T 8 disagree");
}

#[test]
fn reads_that_do_not_map_yield_undetermined_rather_than_unstranded() {
    // Salmon's own inference returns `IU` from an all-zero count array. A
    // detector must not present that as a measurement, so the evidence gates
    // run first and the answer is `undetermined`.
    let reference = Reference::build();
    let mut unrelated = Vec::new();
    for index in 0..2000 {
        let r1 = transcript(0xDEAD_0000 + index as u64, 120);
        let r2 = transcript(0xBEEF_0000 + index as u64, 120);
        let q1 = vec![b'I'; r1.len()];
        let q2 = vec![b'I'; r2.len()];
        unrelated.push(Record::new(PAIRED, index, &r1, &r2).with_quality(&q1, &q2));
    }
    let input = reference.path("unrelated.cbq");
    write_cbq(&input, PAIRED, &unrelated, BLOCK);
    let out = reference.path("strand.json");

    let outcome = bqc_outcome(&[
        "sniff",
        "strand",
        input.to_str().unwrap(),
        "--index",
        reference.index().to_str().unwrap(),
        "--require-confident",
        "--format",
        "json",
        "-o",
        out.to_str().unwrap(),
    ])
    .expect("an undetermined result is not an error");

    let result = &json(&out)["result"];
    assert_eq!(result["strandedness"], "undetermined", "{result}");
    assert_eq!(result["decision"], "undetermined");
    assert_eq!(result["failure_reason"], "insufficient_mapping_evidence");
    // No downstream parameter is manufactured from an answer nobody established.
    assert!(result["salmon_library_type"].is_null());
    assert!(result["featurecounts_strand"].is_null());
    assert!(result["htseq_stranded"].is_null());
    assert_eq!(outcome, bqc::cli::Outcome::NotConfident);
}

#[test]
fn an_unreadable_index_is_refused_with_an_actionable_message() {
    let dir = tempfile::tempdir().expect("tempdir");
    let records: Vec<Record> = (0..100)
        .map(|index| {
            let sequence = transcript(index as u64, 120);
            let quality = vec![b'I'; sequence.len()];
            Record::new(PAIRED, index, &sequence, &sequence).with_quality(&quality, &quality)
        })
        .collect();
    let input = dir.path().join("reads.cbq");
    write_cbq(&input, PAIRED, &records, BLOCK);

    // A directory that is not an index at all.
    let empty = dir.path().join("empty");
    std::fs::create_dir(&empty).expect("mkdir");
    let error = bqc(&[
        "sniff",
        "strand",
        input.to_str().unwrap(),
        "--index",
        empty.to_str().unwrap(),
        "--format",
        "json",
        "-o",
        dir.path().join("out.json").to_str().unwrap(),
    ])
    .unwrap_err();
    let message = format!("{error}");
    assert!(
        message.contains("cannot read the Salmon index"),
        "{message}"
    );
    assert!(message.contains("salmon index -t"), "{message}");

    // A path that is not a directory.
    let error = bqc(&[
        "sniff",
        "strand",
        input.to_str().unwrap(),
        "--index",
        input.to_str().unwrap(),
        "--format",
        "json",
        "-o",
        dir.path().join("out2.json").to_str().unwrap(),
    ])
    .unwrap_err();
    assert!(format!("{error}").contains("is not a directory"), "{error}");
}

#[test]
fn the_index_provenance_is_recorded_with_the_result() {
    let reference = Reference::build();
    let records = reference.pairs(2000, false, PAIRED);
    let input = reference.path("reads.cbq");
    write_cbq(&input, PAIRED, &records, BLOCK);
    let out = reference.path("strand.json");

    bqc(&[
        "sniff",
        "strand",
        input.to_str().unwrap(),
        "--index",
        reference.index().to_str().unwrap(),
        "--min-informative",
        "100",
        "--format",
        "json",
        "-o",
        out.to_str().unwrap(),
    ])
    .expect("sniff strand runs");

    let metadata = &json(&out)["result"]["index_metadata"];
    assert_eq!(metadata["num_refs"], 8);
    assert_eq!(metadata["k"], 31);
    assert_eq!(metadata["has_decoys"], false);
    assert!(!metadata["seq_hash"].as_str().unwrap().is_empty());
    assert!(!metadata["salmon_version"].as_str().unwrap().is_empty());
}

#[test]
fn a_summary_row_is_written_for_cohort_aggregation() {
    let reference = Reference::build();
    let records = reference.pairs(2000, true, PAIRED);
    let input = reference.path("reads.cbq");
    write_cbq(&input, PAIRED, &records, BLOCK);
    let out = reference.path("strand.tsv");

    bqc(&[
        "sniff",
        "strand",
        input.to_str().unwrap(),
        "--index",
        reference.index().to_str().unwrap(),
        "--min-informative",
        "100",
        "--format",
        "tsv",
        "-o",
        out.to_str().unwrap(),
    ])
    .expect("sniff strand runs");

    let tsv = std::fs::read_to_string(&out).expect("read tsv");
    let lines: Vec<&str> = tsv.lines().collect();
    assert_eq!(lines.len(), 2, "a header and exactly one summary row");
    assert!(lines[0].starts_with("input\tdecision\tsalmon_library_type\t"));
    assert!(
        lines[1].contains("\tISR\treverse\tinward\t"),
        "{}",
        lines[1]
    );
}

#[test]
fn a_transcriptome_builds_the_index_on_the_fly() {
    let reference = Reference::build();
    let records = reference.pairs(2000, true, PAIRED);
    let input = reference.path("reads.cbq");
    write_cbq(&input, PAIRED, &records, BLOCK);

    bqc(&[
        "sniff",
        "strand",
        input.to_str().unwrap(),
        "--transcriptome",
        reference.transcriptome().to_str().unwrap(),
        "--min-informative",
        "100",
    ])
    .expect("sniff strand runs with a transcriptome");
}
