// Copyright (c) 2026 The Hiller Lab at the Senckenberg Gesellschaft für Naturforschung
// Distributed under the terms of the GNU General Public License, Version 3.0.

//! Exact deduplication, end to end.
//!
//! Every test passes `--memory-mb 0`, which collapses the Bloom filters to a
//! single bit. That forces every record through the exact-equality path, so the
//! tests exercise the collision-tolerant design rather than a lucky hash.

mod common;

use bqc::io::Schema;
use common::{Record, bqc, read_cbq, write_cbq};

fn se_schema() -> Schema {
    Schema {
        paired: false,
        quality: true,
        headers: true,
        flags: false,
    }
}

fn pe_schema() -> Schema {
    Schema {
        paired: true,
        quality: true,
        headers: true,
        flags: false,
    }
}

fn run(input: &std::path::Path, output: &std::path::Path) {
    bqc(&[
        "dedup",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--memory-mb",
        "0",
    ])
    .unwrap();
}

#[test]
fn dedup_se_keeps_first_occurrence_unchanged() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("in.cbq");
    let output = dir.path().join("out.cbq");
    let first = Record {
        s_seq: b"ACGTACGT".to_vec(),
        s_qual: Some(b"IIIIIIII".to_vec()),
        s_header: Some(b"read_0/1".to_vec()),
        x_seq: None,
        x_qual: None,
        x_header: None,
        flag: None,
    };
    let duplicate = Record {
        s_seq: b"ACGTACGT".to_vec(),
        s_qual: Some(b"JJJJJJJJ".to_vec()),
        s_header: Some(b"read_1/1".to_vec()),
        x_seq: None,
        x_qual: None,
        x_header: None,
        flag: None,
    };
    let unique = Record {
        s_seq: b"TTTTTTTT".to_vec(),
        s_qual: Some(b"IIIIIIII".to_vec()),
        s_header: Some(b"read_2/1".to_vec()),
        x_seq: None,
        x_qual: None,
        x_header: None,
        flag: None,
    };
    write_cbq(&input, se_schema(), &[first.clone(), duplicate, unique.clone()], 512);

    run(&input, &output);

    let (_, records) = read_cbq(&output);
    assert_eq!(records.len(), 2);
    // The earliest occurrence survives, byte-for-byte: different qualities and
    // headers on the duplicate must not matter.
    assert_eq!(records[0], first);
    assert_eq!(records[1], unique);
}

#[test]
fn dedup_pe_uses_the_ordered_pair_as_key() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("in.cbq");
    let output = dir.path().join("out.cbq");
    let pair = |r1: &[u8], r2: &[u8]| Record {
        s_seq: r1.to_vec(),
        s_qual: Some(vec![b'I'; r1.len()]),
        s_header: Some(b"read/1".to_vec()),
        x_seq: Some(r2.to_vec()),
        x_qual: Some(vec![b'I'; r2.len()]),
        x_header: Some(b"read/2".to_vec()),
        flag: None,
    };
    let records = [
        pair(b"AAAA", b"CCCC"),
        pair(b"AAAA", b"CCCC"), // duplicate
        pair(b"AAAA", b"GGGG"), // different R2
    ];
    write_cbq(&input, pe_schema(), &records, 512);

    run(&input, &output);

    let (_, out) = read_cbq(&output);
    assert_eq!(out.len(), 2);
    assert_eq!(out[0], records[0]);
    assert_eq!(out[1], records[2]);
}

#[test]
fn dedup_does_not_alias_across_mate_boundaries() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("in.cbq");
    let output = dir.path().join("out.cbq");
    let pair = |r1: &[u8], r2: &[u8]| Record {
        s_seq: r1.to_vec(),
        s_qual: Some(vec![b'I'; r1.len()]),
        s_header: Some(b"read/1".to_vec()),
        x_seq: Some(r2.to_vec()),
        x_qual: Some(vec![b'I'; r2.len()]),
        x_header: Some(b"read/2".to_vec()),
        flag: None,
    };
    // (AC, CGT) and (ACC, GT) concatenate to the same bytes but are different
    // records; both must survive.
    let records = [pair(b"AC", b"CGT"), pair(b"ACC", b"GT")];
    write_cbq(&input, pe_schema(), &records, 512);

    run(&input, &output);

    let (_, out) = read_cbq(&output);
    assert_eq!(out.len(), 2);
    assert_eq!(out[0], records[0]);
    assert_eq!(out[1], records[1]);
}

#[test]
fn dedup_treats_n_as_an_ordinary_base() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("in.cbq");
    let output = dir.path().join("out.cbq");
    let record = |seq: &[u8], name: &str| Record {
        s_seq: seq.to_vec(),
        s_qual: Some(vec![b'I'; seq.len()]),
        s_header: Some(name.as_bytes().to_vec()),
        x_seq: None,
        x_qual: None,
        x_header: None,
        flag: None,
    };
    let records = [
        record(b"ACGTNACG", "read_0/1"),
        record(b"ACGTNACG", "read_1/1"), // duplicate
        record(b"ACGTNACN", "read_2/1"), // differs at one base
    ];
    write_cbq(&input, se_schema(), &records, 512);

    run(&input, &output);

    let (_, out) = read_cbq(&output);
    assert_eq!(out.len(), 2);
    assert_eq!(out[0], records[0]);
    assert_eq!(out[1], records[2]);
}

#[test]
fn dedup_preserves_input_order() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("in.cbq");
    let output = dir.path().join("out.cbq");
    let records: Vec<Record> = (0..40)
        .map(|i| {
            let seq = vec![b'A'; i / 2 + 1];
            Record {
                s_seq: seq.clone(),
                s_qual: Some(vec![b'I'; seq.len()]),
                s_header: Some(format!("read_{i}/1").into_bytes()),
                x_seq: None,
                x_qual: None,
                x_header: None,
                flag: None,
            }
        })
        .collect();
    write_cbq(&input, se_schema(), &records, 64);

    run(&input, &output);

    // Pairs (0,1) share a sequence, (2,3) share another, and so on; the even
    // records survive, in order.
    let (_, out) = read_cbq(&output);
    assert_eq!(out.len(), 20);
    let expected: Vec<&Record> = records.iter().step_by(2).collect();
    for (actual, expected) in out.iter().zip(expected) {
        assert_eq!(actual, expected);
    }
}

#[test]
fn dedup_thread_parity() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("in.cbq");
    let one = dir.path().join("one.cbq");
    let eight = dir.path().join("eight.cbq");
    let records: Vec<Record> = (0..500)
        .map(|i| {
            // Distinct sequences of varying length, each repeated once, so the
            // dedup decision spans every record and block.
            let seq = vec![b"ACGT"[(i / 2) % 4]; (i / 2) % 20 + 1];
            Record {
                s_seq: seq.clone(),
                s_qual: Some(vec![b'I'; seq.len()]),
                s_header: Some(format!("read_{i}/1").into_bytes()),
                x_seq: None,
                x_qual: None,
                x_header: None,
                flag: None,
            }
        })
        .collect();
    write_cbq(&input, se_schema(), &records, 64);

    for (threads, output) in [("1", &one), ("8", &eight)] {
        bqc(&[
            "dedup",
            input.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
            "--memory-mb",
            "0",
            "-T",
            threads,
        ])
        .unwrap();
    }

    assert_eq!(
        std::fs::read(&one).unwrap(),
        std::fs::read(&eight).unwrap(),
        "dedup output must not depend on thread count"
    );
}
