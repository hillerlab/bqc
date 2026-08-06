// Copyright (c) 2026 The Hiller Lab at the Senckenberg Gesellschaft für Naturforschung
// Distributed under the terms of the GNU General Public License, Version 3.0.

//! Generates the deterministic CBQ benchmark fixtures.
//!
//! ```bash
//! cargo run --release --example make_fixtures -- /tmp/bqc-fixtures
//! ```
//!
//! Every fixture is produced from a fixed seed, so a benchmark run is
//! reproducible across machines and over time. Reads are shaped like real
//! Illumina data in the way that matters for this tool: 3' adapter
//! contamination runs to the end of the read (read-through continues into the
//! index and flow-cell sequence), rather than being a fragment of adapter
//! surrounded by random bases.

use std::fs::File;
use std::path::Path;

use binseq::SequencingRecordBuilder;
use binseq::cbq::{ColumnarBlockWriter, FileHeaderBuilder};

/// `TruSeq` R1 adapter, followed by the index and flow-cell sequence a real read
/// runs into after the adapter.
const R1_READTHROUGH: &[u8] =
    b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAA";
/// The same for R2.
const R2_READTHROUGH: &[u8] =
    b"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTAAAAAAAAAAAAAAAA";
/// An adapter that is not in the built-in detection library.
const UNKNOWN_READTHROUGH: &[u8] =
    b"GTCGATCGTACGGCATCCGATCGTACGATCGGCATTAGCGCATTAGCGGATCGATCGTAAAAAAAAAAAAAAAA";

/// xorshift64*: deterministic and dependency free.
struct Rng(u64);

impl Rng {
    fn new(seed: u64) -> Self {
        Self(seed | 1)
    }

    fn next(&mut self) -> u64 {
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        self.0.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }

    fn below(&mut self, bound: usize) -> usize {
        (self.next() >> 33) as usize % bound
    }

    fn base(&mut self) -> u8 {
        b"ACGT"[(self.next() >> 62) as usize]
    }

    fn sequence(&mut self, length: usize) -> Vec<u8> {
        (0..length).map(|_| self.base()).collect()
    }
}

/// How a fixture's quality strings are shaped.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Quality {
    /// No quality column at all.
    Absent,
    /// High throughout, with the mild 3' decay of a healthy run.
    Good,
    /// A long degraded 3' tail plus a fraction of globally bad reads, so
    /// quality trimming and quality filters both do real work.
    Degraded,
}

/// One fixture's shape.
struct Spec {
    name: &'static str,
    records: usize,
    paired: bool,
    quality: Quality,
    headers: bool,
    /// Fraction of reads carrying adapter read-through, in percent.
    adapter_percent: usize,
    /// Read-through sequence appended after the insert.
    readthrough: (&'static [u8], &'static [u8]),
    /// Read length, or `None` for variable-length reads.
    read_length: Option<usize>,
    /// Give one in eight contaminated reads an indel inside the adapter.
    indels: bool,
    /// Build R2 as the reverse complement of the insert so the mates genuinely
    /// overlap, which is what paired-overlap inference needs.
    overlapping: bool,
    /// Plant miscalls at low-quality positions, giving correction work to do.
    miscalls: bool,
    seed: u64,
}

impl Spec {
    fn quality_string(&self, length: usize, index: usize, rng: &mut Rng) -> Option<Vec<u8>> {
        match self.quality {
            Quality::Absent => None,
            Quality::Good => Some(
                (0..length)
                    .map(|position| {
                        let decay = (position * 8 / length.max(1)) as u8;
                        b'!' + 38 - decay + (rng.next() % 3) as u8
                    })
                    .collect(),
            ),
            Quality::Degraded => {
                let bad_read = index.is_multiple_of(7);
                Some(
                    (0..length)
                        .map(|position| {
                            let tail = position * 3 >= length * 2;
                            let score = match (bad_read, tail) {
                                (true, _) => 8 + (rng.next() % 8) as u8,
                                (false, true) => 14 + (rng.next() % 14) as u8,
                                (false, false) => 34 + (rng.next() % 7) as u8,
                            };
                            b'!' + score
                        })
                        .collect(),
                )
            }
        }
    }
}

/// Builds one mate: a random insert, then adapter read-through, cut to length.
fn mate(spec: &Spec, rng: &mut Rng, length: usize, contaminated: bool, r2: bool) -> Vec<u8> {
    if !contaminated {
        return rng.sequence(length);
    }
    // Short inserts leave more adapter visible; keep a realistic spread.
    let insert = 20 + rng.below(length.saturating_sub(19).max(1));
    let mut read = rng.sequence(insert.min(length));
    let mut tail = if r2 {
        spec.readthrough.1.to_vec()
    } else {
        spec.readthrough.0.to_vec()
    };
    if spec.indels && rng.next().is_multiple_of(8) {
        if rng.next().is_multiple_of(2) {
            tail.insert(6 + rng.below(12), b'T');
        } else {
            tail.remove(6 + rng.below(12));
        }
    }
    read.extend_from_slice(&tail);
    read.truncate(length);
    while read.len() < length {
        read.push(rng.base()); // adapter shorter than the read: pad
    }
    read
}

/// Builds a genuinely overlapping pair from one insert.
///
/// R1 reads the insert's 5' end and R2 its 3' end, reverse complemented. An
/// insert shorter than twice the read length leaves the mates overlapping by
/// `2 * length - insert` bases; an insert shorter than the read length leaves
/// adapter read-through visible in both mates.
fn overlapping_pair(
    spec: &Spec,
    rng: &mut Rng,
    length: usize,
    contaminated: bool,
) -> (Vec<u8>, Option<Vec<u8>>) {
    // Uncontaminated inserts overlap by 40..=length bases; contaminated ones are
    // shorter than a read, so both mates run into the adapter.
    let insert_length = if contaminated {
        30 + rng.below(length - 40)
    } else {
        length + 40 + rng.below(length.saturating_sub(80).max(1))
    };
    let insert = rng.sequence(insert_length);
    let finish = |mut read: Vec<u8>, adapter: &[u8], rng: &mut Rng| {
        if read.len() < length {
            read.extend_from_slice(adapter);
            read.truncate(length);
            while read.len() < length {
                read.push(rng.base());
            }
        } else {
            read.truncate(length);
        }
        read
    };
    let r1 = finish(insert.clone(), spec.readthrough.0, rng);
    let reversed: Vec<u8> = insert.iter().rev().map(|&base| complement(base)).collect();
    let r2 = finish(reversed, spec.readthrough.1, rng);
    (r1, Some(r2))
}

/// Flips a base at a few of the read's lowest-quality positions.
fn plant_miscalls(read: &mut [u8], quality: &[u8], rng: &mut Rng) {
    const BASES: [u8; 4] = *b"ACGT";
    for position in 0..read.len() {
        // Miscall probability tracks the quality score, as on a real instrument.
        let score = quality[position].saturating_sub(b'!');
        let odds = match score {
            0..=9 => 4,
            10..=19 => 12,
            20..=29 => 60,
            _ => 400,
        };
        if rng.below(odds) == 0 {
            let current = read[position];
            let mut replacement = BASES[(rng.next() >> 62) as usize];
            if replacement == current {
                replacement = if current == b'A' { b'C' } else { b'A' };
            }
            read[position] = replacement;
        }
    }
}

fn complement(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        other => other,
    }
}

fn write_fixture(directory: &Path, spec: &Spec) -> binseq::Result<u64> {
    let path = directory.join(format!("{}.cbq", spec.name));
    let header = FileHeaderBuilder::default()
        .is_paired(spec.paired)
        .with_qualities(spec.quality != Quality::Absent)
        .with_headers(spec.headers)
        .with_flags(false)
        .with_block_size(1 << 20)
        .with_compression_level(3)
        .build();
    let mut writer = ColumnarBlockWriter::new(File::create(&path)?, header)?;
    let mut rng = Rng::new(spec.seed);

    for index in 0..spec.records {
        let length = spec
            .read_length
            .unwrap_or_else(|| 40 + rng.below(120).min(111));
        let contaminated = rng.below(100) < spec.adapter_percent;
        let (r1, r2) = if spec.paired && spec.overlapping {
            overlapping_pair(spec, &mut rng, length, contaminated)
        } else {
            (
                mate(spec, &mut rng, length, contaminated, false),
                spec.paired
                    .then(|| mate(spec, &mut rng, length, contaminated, true)),
            )
        };

        let mut s_qual = spec.quality_string(r1.len(), index, &mut rng);
        let mut x_qual = r2
            .as_ref()
            .and_then(|r2| spec.quality_string(r2.len(), index, &mut rng));
        // Overlapping fixtures also carry miscalls at low-quality positions,
        // which is what base correction exists to repair.
        let (mut r1, mut r2) = (r1, r2);
        if spec.miscalls
            && let (Some(r2), Some(q1), Some(q2)) = (r2.as_mut(), s_qual.as_mut(), x_qual.as_mut())
        {
            plant_miscalls(&mut r1, q1, &mut rng);
            plant_miscalls(r2, q2, &mut rng);
        }
        let s_header = spec
            .headers
            .then(|| format!("SIM:1:FC:1:1:{index}:{index} 1:N:0:ATCACG").into_bytes());
        let x_header = (spec.headers && spec.paired)
            .then(|| format!("SIM:1:FC:1:1:{index}:{index} 2:N:0:ATCACG").into_bytes());

        let record = SequencingRecordBuilder::default()
            .s_seq(&r1)
            .opt_s_qual(s_qual.as_deref())
            .opt_s_header(s_header.as_deref())
            .opt_x_seq(r2.as_deref())
            .opt_x_qual(x_qual.as_deref())
            .opt_x_header(x_header.as_deref())
            .build()?;
        writer.push(record)?;
    }
    writer.finish()?;
    Ok(std::fs::metadata(&path)?.len())
}

// One table of fixture shapes; splitting it would only scatter the data.
#[allow(clippy::too_many_lines)]
fn main() -> binseq::Result<()> {
    let directory = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "/tmp/bqc-fixtures".to_string());
    let directory = Path::new(&directory);
    std::fs::create_dir_all(directory)?;

    let truseq = (R1_READTHROUGH, R2_READTHROUGH);
    let base = Spec {
        name: "",
        records: 0,
        paired: false,
        quality: Quality::Good,
        headers: true,
        adapter_percent: 0,
        readthrough: truseq,
        read_length: Some(150),
        indels: false,
        overlapping: false,
        miscalls: false,
        seed: 1,
    };

    let specs = [
        // Single-end and paired-end 150 bp, moderate contamination.
        Spec {
            name: "se150",
            records: 500_000,
            adapter_percent: 20,
            seed: 0x5E15,
            ..base
        },
        Spec {
            name: "pe150",
            records: 250_000,
            paired: true,
            adapter_percent: 20,
            seed: 0x9E15,
            ..base
        },
        // Variable-length reads.
        Spec {
            name: "varlen",
            records: 400_000,
            read_length: None,
            adapter_percent: 20,
            seed: 0x7A21,
            ..base
        },
        // Adapter prevalence extremes.
        Spec {
            name: "adapter-heavy",
            records: 400_000,
            adapter_percent: 85,
            seed: 0xADEA,
            ..base
        },
        Spec {
            name: "adapter-light",
            records: 400_000,
            adapter_percent: 3,
            seed: 0x_AD_11,
            ..base
        },
        // Quality work: long degraded tails and some globally bad reads.
        Spec {
            name: "qual-heavy",
            records: 400_000,
            quality: Quality::Degraded,
            adapter_percent: 20,
            seed: 0x0BAD,
            ..base
        },
        // Sequence-only CBQ.
        Spec {
            name: "noqual",
            records: 400_000,
            quality: Quality::Absent,
            adapter_percent: 20,
            seed: 0x0F00,
            ..base
        },
        // Header-free CBQ, to exercise the schema-preserving writer path.
        Spec {
            name: "noheader",
            records: 400_000,
            headers: false,
            adapter_percent: 20,
            seed: 0x0DEA,
            ..base
        },
        // Genuinely overlapping mates for paired-overlap inference.
        Spec {
            name: "overlap",
            records: 200_000,
            paired: true,
            overlapping: true,
            miscalls: true,
            // Degraded tails plus quality-tracking miscalls are what make the
            // pair disagree, which is the input base correction acts on.
            quality: Quality::Degraded,
            adapter_percent: 40,
            seed: 0x0FAB,
            ..base
        },
        // Adapter indels. Smaller: the banded alignment is far more expensive.
        Spec {
            name: "indel",
            records: 100_000,
            adapter_percent: 40,
            indels: true,
            seed: 0x1DEA,
            ..base
        },
        // An adapter outside the detection library, for the consensus path.
        Spec {
            name: "unknown-adapter",
            records: 200_000,
            adapter_percent: 40,
            readthrough: (UNKNOWN_READTHROUGH, UNKNOWN_READTHROUGH),
            seed: 0x0BEE,
            ..base
        },
    ];

    for spec in &specs {
        let bytes = write_fixture(directory, spec)?;
        println!(
            "{:<16} {:>9} records  {:>4} MB  paired={} quality={} headers={}",
            spec.name,
            spec.records,
            bytes / (1024 * 1024),
            spec.paired,
            spec.quality != Quality::Absent,
            spec.headers,
        );
    }
    Ok(())
}
