// Copyright (c) 2026 The Hiller Lab at the Senckenberg Gesellschaft für Naturforschung
// Distributed under the terms of the GNU General Public License, Version 3.0.

//! Fixtures for the integration tests.
//!
//! CBQ files are written and read back with `binseq`'s own streaming APIs, so
//! the tests validate `bqc` against the format rather than against its own
//! reader. The streaming reader visits blocks in file order, which makes it a
//! valid oracle for order preservation.

#![allow(dead_code)]

use std::fs::File;
use std::path::Path;

use binseq::cbq::{BlockRange, ColumnarBlockWriter, FileHeaderBuilder, Reader};
use binseq::{BinseqRecord, SequencingRecordBuilder};
use bqc::cli::{Cli, Outcome, run};
use bqc::io::Schema;
use clap::Parser;

/// An owned record used to build and compare fixtures.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Record {
    pub s_seq: Vec<u8>,
    pub s_qual: Option<Vec<u8>>,
    pub s_header: Option<Vec<u8>>,
    pub x_seq: Option<Vec<u8>>,
    pub x_qual: Option<Vec<u8>>,
    pub x_header: Option<Vec<u8>>,
    pub flag: Option<u64>,
}

impl Record {
    /// Builds a record matching `schema` from sequences and an index.
    pub fn new(schema: Schema, index: usize, s_seq: &[u8], x_seq: &[u8]) -> Self {
        Self {
            s_seq: s_seq.to_vec(),
            s_qual: schema.quality.then(|| vec![b'I'; s_seq.len()]),
            s_header: schema
                .headers
                .then(|| format!("read_{index}/1").into_bytes()),
            x_seq: schema.paired.then(|| x_seq.to_vec()),
            x_qual: (schema.paired && schema.quality).then(|| vec![b'I'; x_seq.len()]),
            x_header: (schema.paired && schema.headers)
                .then(|| format!("read_{index}/2").into_bytes()),
            flag: schema.flags.then_some(index as u64 * 7 + 1),
        }
    }

    /// Replaces the quality string of R1 (and R2, when paired).
    pub fn with_quality(mut self, s_qual: &[u8], x_qual: &[u8]) -> Self {
        if self.s_qual.is_some() {
            self.s_qual = Some(s_qual.to_vec());
        }
        if self.x_qual.is_some() {
            self.x_qual = Some(x_qual.to_vec());
        }
        self
    }
}

/// Writes a CBQ file with the given schema and block size.
pub fn write_cbq(path: &Path, schema: Schema, records: &[Record], block_size: usize) {
    let header = FileHeaderBuilder::default()
        .is_paired(schema.paired)
        .with_qualities(schema.quality)
        .with_headers(schema.headers)
        .with_flags(schema.flags)
        .with_block_size(block_size)
        .with_compression_level(1)
        .build();
    let file = File::create(path).expect("create fixture");
    let mut writer = ColumnarBlockWriter::new(file, header).expect("cbq writer");
    for record in records {
        let built = SequencingRecordBuilder::default()
            .s_seq(&record.s_seq)
            .opt_s_qual(record.s_qual.as_deref())
            .opt_s_header(record.s_header.as_deref())
            .opt_x_seq(record.x_seq.as_deref())
            .opt_x_qual(record.x_qual.as_deref())
            .opt_x_header(record.x_header.as_deref())
            .opt_flag(record.flag)
            .build()
            .expect("build fixture record");
        writer.push(built).expect("push fixture record");
    }
    writer.finish().expect("finish fixture");
}

/// Reads the presence flags straight out of the CBQ file header.
///
/// This deliberately avoids `binseq`'s readers: it is an independent oracle for
/// schema preservation, and it works on files with zero records, which
/// `MmapReader` cannot currently open.
pub fn schema_of(path: &Path) -> Schema {
    let bytes = std::fs::read(path).expect("read cbq");
    assert!(bytes.len() >= 64, "not a CBQ file");
    assert_eq!(&bytes[..7], b"CBQFILE", "not a CBQ file");
    let flags = u64::from_le_bytes(bytes[8..16].try_into().unwrap());
    Schema {
        paired: flags & (1 << 0) != 0,
        quality: flags & (1 << 1) != 0,
        headers: flags & (1 << 2) != 0,
        flags: flags & (1 << 3) != 0,
    }
}

/// Reads a CBQ file back in file order.
pub fn read_cbq(path: &Path) -> (Schema, Vec<Record>) {
    let schema = schema_of(path);
    let file = File::open(path).expect("open cbq");
    let mut reader = Reader::new(file).expect("cbq reader");
    let mut records = Vec::new();
    let mut cumulative = 0u64;
    while let Some(block_header) = reader.read_block().expect("read block") {
        cumulative += block_header.num_records;
        reader.block.decompress_columns().expect("decompress block");
        for record in reader.block.iter_records(BlockRange::new(0, cumulative)) {
            records.push(Record {
                s_seq: record.sseq().to_vec(),
                s_qual: schema.quality.then(|| record.squal().to_vec()),
                s_header: schema.headers.then(|| record.sheader().to_vec()),
                x_seq: schema.paired.then(|| record.xseq().to_vec()),
                x_qual: (schema.paired && schema.quality).then(|| record.xqual().to_vec()),
                x_header: (schema.paired && schema.headers).then(|| record.xheader().to_vec()),
                flag: schema.flags.then(|| record.flag().expect("flag column")),
            });
        }
    }
    (schema, records)
}

/// Number of blocks in a CBQ file.
pub fn count_blocks(path: &Path) -> usize {
    let file = File::open(path).expect("open cbq");
    let mut reader = Reader::new(file).expect("cbq reader");
    let mut blocks = 0;
    while reader.read_block().expect("read block").is_some() {
        blocks += 1;
    }
    blocks
}

/// Runs `bqc` in-process with the given arguments.
pub fn bqc(args: &[&str]) -> bqc::Result<()> {
    bqc_outcome(args).map(|_| ())
}

/// Runs `bqc` in-process, keeping the outcome so exit statuses can be
/// asserted. `--quiet` is not appended: not every command accepts it.
pub fn bqc_outcome(args: &[&str]) -> bqc::Result<Outcome> {
    let mut command = vec!["bqc"];
    command.extend_from_slice(args);
    if !matches!(args.first(), Some(&"sniff")) {
        command.push("--quiet");
    }
    let cli = Cli::try_parse_from(command).expect("arguments parse");
    run(&cli)
}

/// Runs `bqc`, asserting success.
pub fn bqc_ok(args: &[&str]) {
    if let Err(error) = bqc(args) {
        panic!("bqc {args:?} failed: {error}");
    }
}

/// The 16 combinations of CBQ presence flags.
pub fn all_schemas() -> Vec<Schema> {
    let mut schemas = Vec::new();
    for paired in [false, true] {
        for quality in [false, true] {
            for headers in [false, true] {
                for flags in [false, true] {
                    schemas.push(Schema {
                        paired,
                        quality,
                        headers,
                        flags,
                    });
                }
            }
        }
    }
    schemas
}

/// A readable description of a schema, used in assertion messages.
pub fn describe(schema: Schema) -> String {
    format!(
        "paired={} quality={} headers={} flags={}",
        schema.paired, schema.quality, schema.headers, schema.flags
    )
}

/// The Watson-Crick complement of a base.
#[must_use]
pub fn complement(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        other => other,
    }
}

/// One planted overlap error, used to build correction fixtures.
#[derive(Clone, Copy)]
pub enum Planted {
    /// No error: the mates agree everywhere.
    Clean,
    /// R2 is wrong and doubtful at this R2 position; R1 should fix it.
    R2Doubtful(usize),
    /// R1 is wrong and doubtful at this R1 position; R2 should fix it.
    R1Doubtful(usize),
    /// R1 is wrong but confident: nothing may be corrected.
    R1Confident(usize),
}

/// Builds overlapping pairs with planted errors.
///
/// R1 covers the first `read_length` bases of a `1.5 x read_length` insert and R2
/// the last, so the mates overlap on the middle half. Planted positions must lie
/// inside that overlap for the correction stage to see them.
pub fn overlapping_pairs(
    schema: Schema,
    read_length: usize,
    plan: &[Planted],
) -> (Vec<Record>, usize) {
    let mut sequences = Sequences::new(0x0C0_FFEE);
    let insert_length = read_length * 3 / 2;
    let overlap_offset = insert_length - read_length;
    let records = plan
        .iter()
        .enumerate()
        .map(|(index, planted)| {
            // Only ACGT: an unplanned N would be a mismatch of its own.
            let insert: Vec<u8> = sequences
                .sequence(insert_length)
                .into_iter()
                .map(|base| if base == b'N' { b'A' } else { base })
                .collect();
            let mut r1 = insert[..read_length].to_vec();
            let mut r2: Vec<u8> = insert[overlap_offset..]
                .iter()
                .rev()
                .map(|&base| complement(base))
                .collect();
            let mut q1 = vec![b'!' + 35; read_length];
            let mut q2 = vec![b'!' + 35; read_length];
            let flip = |base: u8| if base == b'A' { b'C' } else { b'A' };
            match *planted {
                Planted::Clean => {}
                Planted::R2Doubtful(position) => {
                    r2[position] = flip(r2[position]);
                    q2[position] = b'!' + 5;
                }
                Planted::R1Doubtful(position) => {
                    r1[position] = flip(r1[position]);
                    q1[position] = b'!' + 5;
                }
                Planted::R1Confident(position) => r1[position] = flip(r1[position]),
            }
            Record::new(schema, index, &r1, &r2).with_quality(&q1, &q2)
        })
        .collect();
    (records, overlap_offset)
}

/// Deterministic pseudo-random sequences with a fixed seed.
pub struct Sequences(u64);

impl Sequences {
    pub fn new(seed: u64) -> Self {
        Self(seed | 1)
    }

    fn next_u64(&mut self) -> u64 {
        // xorshift64*: deterministic and dependency free.
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        self.0.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }

    /// A sequence of `length` bases over `ACGTN`.
    pub fn sequence(&mut self, length: usize) -> Vec<u8> {
        const BASES: [u8; 8] = *b"ACGTACGN";
        (0..length)
            .map(|_| BASES[(self.next_u64() % 8) as usize])
            .collect()
    }

    /// A quality string of `length` bytes spanning Q0..Q40.
    pub fn quality(&mut self, length: usize) -> Vec<u8> {
        (0..length)
            .map(|_| b'!' + (self.next_u64() % 41) as u8)
            .collect()
    }
}
