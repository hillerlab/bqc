// Copyright (c) 2026 The Hiller Lab at the Senckenberg Gesellschaft für Naturforschung
// Distributed under the terms of the GNU General Public License, Version 3.0.

//! UMI extraction and relocation, end to end.

mod common;

use bqc::io::Schema;
use common::{Record, bqc, read_cbq, write_cbq};

const UMI: &[u8] = b"ACGTACGT";

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

/// Builds a record with explicit sequence, quality and header.
fn se_record(seq: &[u8], header: &str) -> Record {
    Record {
        s_seq: seq.to_vec(),
        s_qual: Some(vec![b'I'; seq.len()]),
        s_header: Some(header.as_bytes().to_vec()),
        x_seq: None,
        x_qual: None,
        x_header: None,
        flag: None,
    }
}

#[test]
fn umi_read1_clips_sequence_and_tags_the_header() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("in.cbq");
    let output = dir.path().join("out.cbq");
    let mut seq = UMI.to_vec();
    seq.extend_from_slice(b"GGGG");
    write_cbq(&input, se_schema(), &[se_record(&seq, "read_0/1")], 512);

    bqc(&[
        "umi",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--umi-location",
        "read1",
        "--umi-length",
        "8",
    ])
    .unwrap();

    let (_, records) = read_cbq(&output);
    assert_eq!(records[0].s_seq, b"GGGG");
    assert_eq!(records[0].s_qual.as_deref(), Some(&b"IIII"[..]));
    assert_eq!(
        records[0].s_header.as_deref(),
        Some(&b"read_0/1:ACGTACGT"[..])
    );
}

#[test]
fn umi_read2_clips_only_r2_and_tags_both() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("in.cbq");
    let output = dir.path().join("out.cbq");
    let mut r2 = b"TTTTGGGG".to_vec();
    r2.extend_from_slice(b"ACGT");
    let records = [Record {
        s_seq: b"AAAACCCC".to_vec(),
        s_qual: Some(vec![b'I'; 8]),
        s_header: Some(b"read_0/1".to_vec()),
        x_seq: Some(r2),
        x_qual: Some(vec![b'I'; 12]),
        x_header: Some(b"read_0/2".to_vec()),
        flag: None,
    }];
    write_cbq(&input, pe_schema(), &records, 512);

    bqc(&[
        "umi",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--umi-location",
        "read2",
        "--umi-length",
        "8",
    ])
    .unwrap();

    let (_, records) = read_cbq(&output);
    assert_eq!(records[0].s_seq, b"AAAACCCC");
    assert_eq!(records[0].x_seq.as_deref(), Some(&b"ACGT"[..]));
    assert_eq!(
        records[0].s_header.as_deref(),
        Some(&b"read_0/1:TTTTGGGG"[..])
    );
    assert_eq!(
        records[0].x_header.as_deref(),
        Some(&b"read_0/2:TTTTGGGG"[..])
    );
}

#[test]
fn umi_index_and_per_index_parse_the_header() {
    let dir = tempfile::tempdir().unwrap();
    let header = "read_0 1:N:0:TATAGCCT+GGTCCCGA";

    // index1
    let input = dir.path().join("idx1.cbq");
    let output = dir.path().join("idx1_out.cbq");
    write_cbq(&input, se_schema(), &[se_record(b"AAAACCCC", header)], 512);
    bqc(&[
        "umi",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--umi-location",
        "index1",
    ])
    .unwrap();
    let (_, records) = read_cbq(&output);
    assert_eq!(records[0].s_seq, b"AAAACCCC");
    assert_eq!(
        records[0].s_header.as_deref(),
        Some(&b"read_0:TATAGCCT 1:N:0:TATAGCCT+GGTCCCGA"[..])
    );

    // per_index (paired)
    let input = dir.path().join("peridx.cbq");
    let output = dir.path().join("peridx_out.cbq");
    let records = [Record {
        s_seq: b"AAAACCCC".to_vec(),
        s_qual: Some(vec![b'I'; 8]),
        s_header: Some(header.as_bytes().to_vec()),
        x_seq: Some(b"TTTTGGGG".to_vec()),
        x_qual: Some(vec![b'I'; 8]),
        x_header: Some(header.as_bytes().to_vec()),
        flag: None,
    }];
    write_cbq(&input, pe_schema(), &records, 512);
    bqc(&[
        "umi",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--umi-location",
        "per_index",
    ])
    .unwrap();
    let (_, records) = read_cbq(&output);
    let expected = b"read_0:TATAGCCT_GGTCCCGA 1:N:0:TATAGCCT+GGTCCCGA";
    assert_eq!(records[0].s_header.as_deref(), Some(&expected[..]));
    assert_eq!(records[0].x_header.as_deref(), Some(&expected[..]));
}

#[test]
fn umi_per_read_clips_both_and_joins_the_tag() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("in.cbq");
    let output = dir.path().join("out.cbq");
    let mut r1 = b"ACGTACGT".to_vec();
    r1.extend_from_slice(b"GGGG");
    let mut r2 = b"TTTTGGGG".to_vec();
    r2.extend_from_slice(b"CCCC");
    let records = [Record {
        s_seq: r1,
        s_qual: Some(vec![b'I'; 12]),
        s_header: Some(b"read_0/1".to_vec()),
        x_seq: Some(r2),
        x_qual: Some(vec![b'I'; 12]),
        x_header: Some(b"read_0/2".to_vec()),
        flag: None,
    }];
    write_cbq(&input, pe_schema(), &records, 512);

    bqc(&[
        "umi",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--umi-location",
        "per_read",
        "--umi-length",
        "8",
    ])
    .unwrap();

    let (_, records) = read_cbq(&output);
    assert_eq!(records[0].s_seq, b"GGGG");
    assert_eq!(records[0].x_seq.as_deref(), Some(&b"CCCC"[..]));
    assert_eq!(
        records[0].s_header.as_deref(),
        Some(&b"read_0/1:ACGTACGT_TTTTGGGG"[..])
    );
    assert_eq!(
        records[0].x_header.as_deref(),
        Some(&b"read_0/2:ACGTACGT_TTTTGGGG"[..])
    );
}

#[test]
fn umi_prefix_and_delimiter_shape_the_tag() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("in.cbq");
    let output = dir.path().join("out.cbq");
    let mut seq = UMI.to_vec();
    seq.extend_from_slice(b"GGGG");
    write_cbq(&input, se_schema(), &[se_record(&seq, "read_0/1")], 512);

    bqc(&[
        "umi",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--umi-location",
        "read1",
        "--umi-length",
        "8",
        "--umi-prefix",
        "UMI",
    ])
    .unwrap();

    let (_, records) = read_cbq(&output);
    assert_eq!(
        records[0].s_header.as_deref(),
        Some(&b"read_0/1:UMI_ACGTACGT"[..])
    );
}

#[test]
fn umi_short_read_is_an_error() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("in.cbq");
    let output = dir.path().join("out.cbq");
    write_cbq(&input, se_schema(), &[se_record(b"ACGT", "read_0/1")], 512);

    let err = bqc(&[
        "umi",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--umi-location",
        "read1",
        "--umi-length",
        "8",
    ])
    .unwrap_err();
    assert!(
        err.to_string().contains("shorter than"),
        "unexpected error: {err}"
    );
    // No partial output is published.
    assert!(!output.exists());
}

#[test]
fn umi_requires_stored_headers() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("in.cbq");
    let output = dir.path().join("out.cbq");
    let schema = Schema {
        paired: false,
        quality: true,
        headers: false,
        flags: false,
    };
    write_cbq(&input, schema, &[se_record(b"AAAACCCC", "ignored")], 512);

    let err = bqc(&[
        "umi",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--umi-location",
        "index1",
    ])
    .unwrap_err();
    assert!(
        err.to_string().contains("requires stored read headers"),
        "unexpected error: {err}"
    );
}

#[test]
fn umi_thread_parity() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("in.cbq");
    let one = dir.path().join("one.cbq");
    let eight = dir.path().join("eight.cbq");
    let records: Vec<Record> = (0..200)
        .map(|i| {
            let mut seq = UMI.to_vec();
            seq.extend_from_slice(&format!("record_{i}").into_bytes());
            se_record(&seq, &format!("read_{i}/1"))
        })
        .collect();
    write_cbq(&input, se_schema(), &records, 64);

    for (threads, output) in [("1", &one), ("8", &eight)] {
        bqc(&[
            "umi",
            input.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
            "--umi-location",
            "read1",
            "--umi-length",
            "8",
            "-T",
            threads,
        ])
        .unwrap();
    }

    assert_eq!(
        std::fs::read(&one).unwrap(),
        std::fs::read(&eight).unwrap(),
        "UMI output must not depend on thread count"
    );
}
