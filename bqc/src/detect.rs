// Copyright (c) 2026 The Hiller Lab at the Senckenberg Gesellschaft für Naturforschung
// Distributed under the terms of the GNU General Public License, Version 3.0.

//! The known-adapter library, and the index that makes it cheap to search.
//!
//! This module holds the reference data and the primitives that read it. The
//! detector itself lives in [`crate::sniff::adapters`] — there is exactly one,
//! shared by `bqc sniff adapters`, which reports its findings, and by
//! `--auto-detect`, which consumes its recommendation. Both therefore always
//! agree about what is in a file.
//!
//! # The seed index
//!
//! The common exact case uses each adapter's 5' 12-mer: one rolling hash over a
//! read proposes plausible entries and coordinates for final verification.
//! Errors inside that seed, terminal partials and indels are covered by a
//! lossless set of shorter exact q-grams whose positions bound the possible
//! alignment. Matcher settings that cannot be indexed safely use the
//! exhaustive verifier instead, so indexing changes cost rather than results.
//!
//! # Voting
//!
//! `vote` and `majority_base` are the column arithmetic behind consensus
//! extension. They live here rather than beside the detector because the
//! thresholds they encode — half the voters must cover a column, and one base
//! must hold 60% of the votes — describe how much evidence is enough to call a
//! base, which is a property of the data rather than of any one caller.

use crate::adapter::{
    find_three_prime, indel_best_at, indel_improves, verify_at, Adapter, AdapterHit, AdapterParams,
};

/// Version of the built-in known-adapter library.
pub const KNOWN_ADAPTERS_VERSION: u32 = 1;

/// One entry of the known-adapter library.
pub struct KnownAdapter {
    pub name: &'static str,
    pub r1: &'static [u8],
    /// Mate-2 sequence; `None` means both mates use the R1 sequence.
    pub r2: Option<&'static [u8]>,
}

impl KnownAdapter {
    /// The sequence expected on `mate`.
    #[must_use]
    pub fn sequence(&self, mate: crate::process::Mate) -> &'static [u8] {
        match mate {
            crate::process::Mate::R1 => self.r1,
            crate::process::Mate::R2 => self.r2.unwrap_or(self.r1),
        }
    }

    /// Broad pipeline-facing category derived from the catalogue label.
    ///
    /// The bundled source distinguishes adapters, primers and PhiX/control
    /// sequences in its names but carries no separate taxonomy. Keeping this
    /// intentionally coarse avoids duplicating 234 records just to serialize the
    /// distinction already present in the catalogue.
    #[must_use]
    pub fn category(&self) -> &'static str {
        if self.name.contains("phix") {
            "control"
        } else if self.name.contains("primer") {
            "primer"
        } else {
            "adapter"
        }
    }
}

/// The built-in known-adapter library, version 1.
///
/// Includes the adapter and primer sequences catalogued by fastp, with
/// normalized lowercase kebab-case identifiers.
pub const KNOWN_ADAPTERS: &[KnownAdapter] = &[
    KnownAdapter {
        name: "illumina-truseq",
        r1: b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
        r2: Some(b"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"),
    },
    KnownAdapter {
        name: "nextera",
        r1: b"CTGTCTCTTATACACATCT",
        r2: None,
    },
    KnownAdapter {
        name: "small-rna",
        r1: b"TGGAATTCTCGGGTGCCAAGG",
        r2: Some(b"GATCGTCGGACTGTAGAACTCTGAACGTGTAGA"),
    },

    // General Illumina, TruSeq, PCR, and PhiX sequences.
    KnownAdapter {
        name: "illumina-expression-pcr-primer-2",
        r1: b"AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer",
        r1: b"AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGA",
        r2: None,
    },
    KnownAdapter {
        name: "truseq-universal-adapter",
        r1: b"AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-pcr-dimer",
        r1: b"AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-pcr-primer-pair",
        r1: b"AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTCAAGCAGAAGACGGCATACGAGCTCTTCCGATCT",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-paired-end-adapter-1",
        r1: b"ACACTCTTTCCCTACACGACGCTCTTCCGATCT",
        r2: None,
    },
    KnownAdapter {
        name: "truseq3-indexed-adapter",
        r1: b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
        r2: None,
    },
    KnownAdapter {
        name: "truseq-reverse-adapter",
        r1: b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "truseq2-paired-end-reverse",
        r1: b"AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG",
        r2: None,
    },
    KnownAdapter {
        name: "pcr-primer-2-reverse-complement",
        r1: b"AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "phix-read-1-adapter",
        r1: b"AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTGAAA",
        r2: None,
    },
    KnownAdapter {
        name: "truseq3-universal-adapter",
        r1: b"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA",
        r2: None,
    },
    KnownAdapter {
        name: "pcr-primer-1-reverse-complement",
        r1: b"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT",
        r2: None,
    },
    KnownAdapter {
        name: "phix-read-2-adapter",
        r1: b"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTAAAAAA",
        r2: None,
    },
    KnownAdapter {
        name: "truseq2-single-end",
        r1: b"AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },

    // Illumina RNA PCR primers and index primers.
    KnownAdapter {
        name: "illumina-rna-pcr-primer-index-35",
        r1: b"CAAGCAGAAGACGGCATACGAGATAAAATGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-pcr-primer-index-10",
        r1: b"CAAGCAGAAGACGGCATACGAGATAAGCTAGTGACTGGAGTTC",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-index-10",
        r1: b"CAAGCAGAAGACGGCATACGAGATAAGCTAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-pcr-primer-index-2",
        r1: b"CAAGCAGAAGACGGCATACGAGATACATCGGTGACTGGAGTTC",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-index-2",
        r1: b"CAAGCAGAAGACGGCATACGAGATACATCGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-index-38",
        r1: b"CAAGCAGAAGACGGCATACGAGATAGCTAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-index-27",
        r1: b"CAAGCAGAAGACGGCATACGAGATAGGAATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-index-25",
        r1: b"CAAGCAGAAGACGGCATACGAGATATCAGTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-index-31",
        r1: b"CAAGCAGAAGACGGCATACGAGATATCGTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-index-44",
        r1: b"CAAGCAGAAGACGGCATACGAGATATTATAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-index-37",
        r1: b"CAAGCAGAAGACGGCATACGAGATATTCCGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-pcr-primer-index-6",
        r1: b"CAAGCAGAAGACGGCATACGAGATATTGGCGTGACTGGAGTTC",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-index-6",
        r1: b"CAAGCAGAAGACGGCATACGAGATATTGGCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-pcr-primer-index-5",
        r1: b"CAAGCAGAAGACGGCATACGAGATCACTGTGTGACTGGAGTTC",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-index-5",
        r1: b"CAAGCAGAAGACGGCATACGAGATCACTGTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-index-23",
        r1: b"CAAGCAGAAGACGGCATACGAGATCCACTCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-index-30",
        r1: b"CAAGCAGAAGACGGCATACGAGATCCGGTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-index-21",
        r1: b"CAAGCAGAAGACGGCATACGAGATCGAAACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-index-42",
        r1: b"CAAGCAGAAGACGGCATACGAGATCGATTAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-index-33",
        r1: b"CAAGCAGAAGACGGCATACGAGATCGCCTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-paired-end-pcr-primer-2",
        r1: b"CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-index-22",
        r1: b"CAAGCAGAAGACGGCATACGAGATCGTACGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-pcr-primer-index-1",
        r1: b"CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTC",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-index-1",
        r1: b"CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-index-17",
        r1: b"CAAGCAGAAGACGGCATACGAGATCTCTACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-pcr-primer-index-9",
        r1: b"CAAGCAGAAGACGGCATACGAGATCTGATCGTGACTGGAGTTC",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-index-9",
        r1: b"CAAGCAGAAGACGGCATACGAGATCTGATCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-index-47",
        r1: b"CAAGCAGAAGACGGCATACGAGATCTTCGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-index-28",
        r1: b"CAAGCAGAAGACGGCATACGAGATCTTTTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-index-45",
        r1: b"CAAGCAGAAGACGGCATACGAGATGAATGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-pcr-primer-index-7",
        r1: b"CAAGCAGAAGACGGCATACGAGATGATCTGGTGACTGGAGTTC",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-index-7",
        r1: b"CAAGCAGAAGACGGCATACGAGATGATCTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-index-34",
        r1: b"CAAGCAGAAGACGGCATACGAGATGCCATGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-pcr-primer-index-3",
        r1: b"CAAGCAGAAGACGGCATACGAGATGCCTAAGTGACTGGAGTTC",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-index-3",
        r1: b"CAAGCAGAAGACGGCATACGAGATGCCTAAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-index-18",
        r1: b"CAAGCAGAAGACGGCATACGAGATGCGGACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-index-24",
        r1: b"CAAGCAGAAGACGGCATACGAGATGCTACCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-index-26",
        r1: b"CAAGCAGAAGACGGCATACGAGATGCTCATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-index-43",
        r1: b"CAAGCAGAAGACGGCATACGAGATGCTGTAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-index-14",
        r1: b"CAAGCAGAAGACGGCATACGAGATGGAACTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-index-16",
        r1: b"CAAGCAGAAGACGGCATACGAGATGGACGGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-index-20",
        r1: b"CAAGCAGAAGACGGCATACGAGATGGCCACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-pcr-primer-index-11",
        r1: b"CAAGCAGAAGACGGCATACGAGATGTAGCCGTGACTGGAGTTC",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-index-11",
        r1: b"CAAGCAGAAGACGGCATACGAGATGTAGCCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-index-39",
        r1: b"CAAGCAGAAGACGGCATACGAGATGTATAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-index-41",
        r1: b"CAAGCAGAAGACGGCATACGAGATGTCGTCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-pcr-primer-index-12",
        r1: b"CAAGCAGAAGACGGCATACGAGATTACAAGGTGACTGGAGTTC",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-index-12",
        r1: b"CAAGCAGAAGACGGCATACGAGATTACAAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-index-29",
        r1: b"CAAGCAGAAGACGGCATACGAGATTAGTTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-pcr-primer-index-8",
        r1: b"CAAGCAGAAGACGGCATACGAGATTCAAGTGTGACTGGAGTTC",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-index-8",
        r1: b"CAAGCAGAAGACGGCATACGAGATTCAAGTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-index-46",
        r1: b"CAAGCAGAAGACGGCATACGAGATTCGGGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-index-40",
        r1: b"CAAGCAGAAGACGGCATACGAGATTCTGAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-index-15",
        r1: b"CAAGCAGAAGACGGCATACGAGATTGACATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-index-32",
        r1: b"CAAGCAGAAGACGGCATACGAGATTGAGTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-index-48",
        r1: b"CAAGCAGAAGACGGCATACGAGATTGCCGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-pcr-primer-index-4",
        r1: b"CAAGCAGAAGACGGCATACGAGATTGGTCAGTGACTGGAGTTC",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-index-4",
        r1: b"CAAGCAGAAGACGGCATACGAGATTGGTCAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-index-36",
        r1: b"CAAGCAGAAGACGGCATACGAGATTGTTGGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-index-13",
        r1: b"CAAGCAGAAGACGGCATACGAGATTTGACTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-index-19",
        r1: b"CAAGCAGAAGACGGCATACGAGATTTTCACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
        r2: None,
    },

    // Single-end, SOLiD, Nextera i7, and general sequencing primers.
    KnownAdapter {
        name: "illumina-single-end-adapter-2",
        r1: b"CAAGCAGAAGACGGCATACGAGCTCTTCCGATCT",
        r2: None,
    },
    KnownAdapter {
        name: "abi-solid3-adapter-b",
        r1: b"CCACTACGCCTCCGCTTTCCTCTCTATGGGCAGTCGGTGAT",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-nlaiii-sequencing-primer",
        r1: b"CCGACAGGTTCAGAGTTCTACAGTCCGACATG",
        r2: None,
    },
    KnownAdapter {
        name: "nextera-i7-primer-n711",
        r1: b"CCGAGCCCACGAGACAAGAGGCAATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "nextera-i7-primer-n716",
        r1: b"CCGAGCCCACGAGACACTCGCTAATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "nextera-i7-primer-n724",
        r1: b"CCGAGCCCACGAGACACTGAGCGATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "nextera-i7-primer-n703",
        r1: b"CCGAGCCCACGAGACAGGCAGAAATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "nextera-i7-primer-n715",
        r1: b"CCGAGCCCACGAGACATCTCAGGATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "nextera-i7-primer-n722",
        r1: b"CCGAGCCCACGAGACATGCGCAGATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "nextera-i7-primer-n708",
        r1: b"CCGAGCCCACGAGACCAGAGAGGATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "nextera-i7-primer-n726",
        r1: b"CCGAGCCCACGAGACCCTAAGACATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "nextera-i7-primer-n710",
        r1: b"CCGAGCCCACGAGACCGAGGCTGATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "nextera-i7-primer-n727",
        r1: b"CCGAGCCCACGAGACCGATCAGTATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "nextera-i7-primer-n720",
        r1: b"CCGAGCCCACGAGACCGGAGCCTATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "nextera-i7-primer-n702",
        r1: b"CCGAGCCCACGAGACCGTACTAGATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "nextera-i7-primer-n707",
        r1: b"CCGAGCCCACGAGACCTCTCTACATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "nextera-i7-primer-n719",
        r1: b"CCGAGCCCACGAGACGCGTAGTAATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "nextera-i7-primer-n709",
        r1: b"CCGAGCCCACGAGACGCTACGCTATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "nextera-i7-primer-n714",
        r1: b"CCGAGCCCACGAGACGCTCATGAATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "nextera-i7-primer-n705",
        r1: b"CCGAGCCCACGAGACGGACTCCTATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "nextera-i7-primer-n718",
        r1: b"CCGAGCCCACGAGACGGAGCTACATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "nextera-i7-primer-n712",
        r1: b"CCGAGCCCACGAGACGTAGAGGAATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "nextera-i7-primer-n701",
        r1: b"CCGAGCCCACGAGACTAAGGCGAATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "nextera-i7-primer-n721",
        r1: b"CCGAGCCCACGAGACTACGCTGCATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "nextera-i7-primer-n723",
        r1: b"CCGAGCCCACGAGACTAGCGCTCATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "nextera-i7-primer-n706",
        r1: b"CCGAGCCCACGAGACTAGGCATGATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "nextera-i7-primer-n704",
        r1: b"CCGAGCCCACGAGACTCCTGAGCATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "nextera-i7-primer-n729",
        r1: b"CCGAGCCCACGAGACTCGACGTCATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "nextera-i7-primer-n728",
        r1: b"CCGAGCCCACGAGACTGCAGCTAATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-small-rna-sequencing-primer",
        r1: b"CGACAGGTTCAGAGTTCTACAGTCCGACGATC",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-paired-end-sequencing-primer-2",
        r1: b"CGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT",
        r2: None,
    },
    KnownAdapter {
        name: "clontech-universal-primer-mix-long",
        r1: b"CTAATACGACTCACTATAGGGCAAGCAGTGGTATCAACGCAGAGT",
        r2: None,
    },
    KnownAdapter {
        name: "nextera-i7-adapter-no-barcode",
        r1: b"CTGAGCGGGCTGGCAAGGCAGACCGATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "nextera-i5-adapter",
        r1: b"CTGATGGCGCGAGGGAGGCGTGTAGATCTCGGTGGTCGCCGTATCATT",
        r2: None,
    },
    KnownAdapter {
        name: "abi-solid3-adapter-a",
        r1: b"CTGCCCCGGGTTCCTCATTCTCTCAGCAGCATG",
        r2: None,
    },
    KnownAdapter {
        name: "nextera-i7-transposase-1",
        r1: b"CTGTCTCTTATACACATCTCCGAGCCCACGAGAC",
        r2: None,
    },
    KnownAdapter {
        name: "nextera-i7-transposase-2",
        r1: b"CTGTCTCTTATACACATCTCTGAGCGGGCTGGCAAGGC",
        r2: None,
    },
    KnownAdapter {
        name: "nextera-i5-transposase-2",
        r1: b"CTGTCTCTTATACACATCTCTGATGGCGCGAGGGAGGC",
        r2: None,
    },
    KnownAdapter {
        name: "nextera-i5-transposase-1",
        r1: b"CTGTCTCTTATACACATCTGACGCTGCCGACGA",
        r2: None,
    },

    // Nextera i5 index primers.
    KnownAdapter {
        name: "nextera-i5-primer-s516",
        r1: b"GACGCTGCCGACGAACTCTAGGGTGTAGATCTCGGTGGTCGCCGTATCATT",
        r2: None,
    },
    KnownAdapter {
        name: "nextera-i5-primer-s503",
        r1: b"GACGCTGCCGACGAAGAGGATAGTGTAGATCTCGGTGGTCGCCGTATCATT",
        r2: None,
    },
    KnownAdapter {
        name: "nextera-i5-primer-s515",
        r1: b"GACGCTGCCGACGAAGCTAGAAGTGTAGATCTCGGTGGTCGCCGTATCATT",
        r2: None,
    },
    KnownAdapter {
        name: "nextera-i5-primer-s508",
        r1: b"GACGCTGCCGACGAAGGCTTAGGTGTAGATCTCGGTGGTCGCCGTATCATT",
        r2: None,
    },
    KnownAdapter {
        name: "nextera-i5-primer-s502",
        r1: b"GACGCTGCCGACGAATAGAGAGGTGTAGATCTCGGTGGTCGCCGTATCATT",
        r2: None,
    },
    KnownAdapter {
        name: "nextera-i5-primer-s520",
        r1: b"GACGCTGCCGACGAATAGCCTTGTGTAGATCTCGGTGGTCGCCGTATCATT",
        r2: None,
    },
    KnownAdapter {
        name: "nextera-i5-primer-s510",
        r1: b"GACGCTGCCGACGAATTAGACGGTGTAGATCTCGGTGGTCGCCGTATCATT",
        r2: None,
    },
    KnownAdapter {
        name: "nextera-i5-primer-s511",
        r1: b"GACGCTGCCGACGACGGAGAGAGTGTAGATCTCGGTGGTCGCCGTATCATT",
        r2: None,
    },
    KnownAdapter {
        name: "nextera-i5-primer-s513",
        r1: b"GACGCTGCCGACGACTAGTCGAGTGTAGATCTCGGTGGTCGCCGTATCATT",
        r2: None,
    },
    KnownAdapter {
        name: "nextera-i5-primer-s505",
        r1: b"GACGCTGCCGACGACTCCTTACGTGTAGATCTCGGTGGTCGCCGTATCATT",
        r2: None,
    },
    KnownAdapter {
        name: "nextera-i5-primer-s518",
        r1: b"GACGCTGCCGACGACTTAATAGGTGTAGATCTCGGTGGTCGCCGTATCATT",
        r2: None,
    },
    KnownAdapter {
        name: "nextera-i5-primer-s501",
        r1: b"GACGCTGCCGACGAGCGATCTAGTGTAGATCTCGGTGGTCGCCGTATCATT",
        r2: None,
    },
    KnownAdapter {
        name: "nextera-i5-primer-s521",
        r1: b"GACGCTGCCGACGATAAGGCTCGTGTAGATCTCGGTGGTCGCCGTATCATT",
        r2: None,
    },
    KnownAdapter {
        name: "nextera-i5-primer-s507",
        r1: b"GACGCTGCCGACGATACTCCTTGTGTAGATCTCGGTGGTCGCCGTATCATT",
        r2: None,
    },
    KnownAdapter {
        name: "nextera-i5-primer-s506",
        r1: b"GACGCTGCCGACGATATGCAGTGTGTAGATCTCGGTGGTCGCCGTATCATT",
        r2: None,
    },
    KnownAdapter {
        name: "nextera-i5-primer-s522",
        r1: b"GACGCTGCCGACGATCGCATAAGTGTAGATCTCGGTGGTCGCCGTATCATT",
        r2: None,
    },
    KnownAdapter {
        name: "nextera-i5-primer-s504",
        r1: b"GACGCTGCCGACGATCTACTCTGTGTAGATCTCGGTGGTCGCCGTATCATT",
        r2: None,
    },
    KnownAdapter {
        name: "nextera-i5-primer-s517",
        r1: b"GACGCTGCCGACGATCTTACGCGTGTAGATCTCGGTGGTCGCCGTATCATT",
        r2: None,
    },

    // TruSeq indexed adapters and related external adapters.
    KnownAdapter {
        name: "nextera-lmp-read-1-external-adapter",
        r1: b"GATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
        r2: None,
    },
    KnownAdapter {
        name: "truseq-adapter-index-5",
        r1: b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "truseq-adapter-index-25-variant",
        r1: b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTGATATATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "truseq-adapter-index-25",
        r1: b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTGATATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "truseq-adapter-index-8",
        r1: b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTTGAATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "truseq-adapter-index-13-variant",
        r1: b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTCAACAATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "truseq-adapter-index-13",
        r1: b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTCAACTCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "truseq-adapter-index-14-variant",
        r1: b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTTCCGTATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "truseq-adapter-index-14",
        r1: b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTTCCGTCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "truseq-adapter-index-1",
        r1: b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "truseq-adapter-index-15-variant",
        r1: b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACATGTCAGAATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "truseq-adapter-index-15",
        r1: b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACATGTCAGTCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "truseq-adapter-index-27-variant",
        r1: b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACATTCCTTTATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "truseq-adapter-index-27",
        r1: b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACATTCCTTTCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "truseq-adapter-index-7",
        r1: b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATCATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "truseq-adapter-index-23",
        r1: b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACCCACTCTTCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "truseq-adapter-index-16-variant",
        r1: b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACCCGTCCCGATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "truseq-adapter-index-16",
        r1: b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACCCGTCCCTCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "truseq-adapter-index-2",
        r1: b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "truseq-adapter-index-22-variant",
        r1: b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGTACGTAATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "truseq-adapter-index-22",
        r1: b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGTACGTTCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "truseq-adapter-index-12",
        r1: b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACCTTGTAATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "truseq-adapter-index-23-variant",
        r1: b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACGAGTGGATATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "truseq-adapter-index-9",
        r1: b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACGATCAGATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "truseq-adapter-index-6",
        r1: b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAATATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "truseq-adapter-index-11",
        r1: b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGCTACATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "truseq-adapter-index-18-variant",
        r1: b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCCGCACATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "truseq-adapter-index-18",
        r1: b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCCGCATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "truseq-adapter-index-19-variant",
        r1: b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGAAACGATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "truseq-adapter-index-19",
        r1: b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGAAACTCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "truseq-adapter-index-20-variant",
        r1: b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGGCCTTATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "truseq-adapter-index-20",
        r1: b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGGCCTTCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "truseq-adapter-index-21-variant",
        r1: b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTTTCGGAATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "truseq-adapter-index-21",
        r1: b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTTTCGGTCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "truseq-adapter-index-10",
        r1: b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACTAGCTTATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "truseq-adapter-index-4",
        r1: b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACTGACCAATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "truseq-adapter-index-3",
        r1: b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACTTAGGCATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-paired-end-adapter-2",
        r1: b"GATCGGAAGAGCGGTTCAGCAGGAATGCCGAG",
        r2: None,
    },
    KnownAdapter {
        name: "nextera-lmp-read-2-external-adapter",
        r1: b"GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-single-end-adapter-1",
        r1: b"GATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "nextera-transposase-2",
        r1: b"GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-prefix-pe-2",
        r1: b"GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-prefix-pe-1",
        r1: b"TACACTCTTTCCCTACACGACGCTCTTCCGATCT",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-rp1",
        r1: b"TCGGACTGTAGAACTCTGAACGTGTAGATCTCGGTGGTCGCCGTATCATT",
        r2: None,
    },
    KnownAdapter {
        name: "nextera-transposase-1",
        r1: b"TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG",
        r2: None,
    },

    // Illumina small-RNA RPI index primers.
    KnownAdapter {
        name: "illumina-rna-pcr-primer-rpi-5",
        r1: b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-rpi-25",
        r1: b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACACTGATATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-rpi-8",
        r1: b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACACTTGAATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-rpi-13",
        r1: b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACAGTCAAATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-rpi-14",
        r1: b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACAGTTCCATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-rpi-1",
        r1: b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-rpi-26",
        r1: b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACATGAGCATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-rpi-15",
        r1: b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACATGTCAATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-rpi-27",
        r1: b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACATTCCTATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-rpi-28",
        r1: b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCAAAAGATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-rpi-29",
        r1: b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCAACTAATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-rpi-30",
        r1: b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCACCGGATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-rpi-31",
        r1: b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCACGATATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-rpi-32",
        r1: b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCACTCAATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-rpi-7",
        r1: b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCAGATCATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-rpi-33",
        r1: b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCAGGCGATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-rpi-34",
        r1: b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCATGGCATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-rpi-35",
        r1: b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCATTTTATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-rpi-36",
        r1: b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCCAACAATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-rpi-16",
        r1: b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCCGTCCATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-rpi-2",
        r1: b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-rpi-37",
        r1: b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCGGAATATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-rpi-22",
        r1: b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCGTACGATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-rpi-38",
        r1: b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCTAGCTATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-rpi-39",
        r1: b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCTATACATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-rpi-40",
        r1: b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCTCAGAATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-rpi-12",
        r1: b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACCTTGTAATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-rpi-41",
        r1: b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGACGACATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-rpi-23",
        r1: b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGAGTGGATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-rpi-9",
        r1: b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGATCAGATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-rpi-6",
        r1: b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGCCAATATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-rpi-11",
        r1: b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGGCTACATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-rpi-24",
        r1: b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGGTAGCATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-rpi-17",
        r1: b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGTAGAGATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-rpi-18",
        r1: b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGTCCGCATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-rpi-19",
        r1: b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGTGAAAATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-rpi-20",
        r1: b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGTGGCCATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-rpi-21",
        r1: b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGTTTCGATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-rpi-42",
        r1: b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTAATCGATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-rpi-43",
        r1: b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTACAGCATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-rpi-10",
        r1: b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTAGCTTATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-rpi-44",
        r1: b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTATAATATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-rpi-45",
        r1: b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTCATTCATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-rpi-46",
        r1: b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTCCCGAATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-rpi-47",
        r1: b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTCGAAGATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-rpi-48",
        r1: b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTCGGCAATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-rpi-4",
        r1: b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTGACCAATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-rna-pcr-primer-rpi-3",
        r1: b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTTAGGCATCTCGTATGCCGTCTTCTGCTTG",
        r2: None,
    },

    // Flow-cell, MGI/BGI, and QIAseq sequences.
    KnownAdapter {
        name: "illumina-flow-cell-1",
        r1: b"TTTTTTTTTTAATGATACGGCGACCACCGAGATCTACAC",
        r2: None,
    },
    KnownAdapter {
        name: "illumina-flow-cell-2",
        r1: b"TTTTTTTTTTCAAGCAGAAGACGGCATACGA",
        r2: None,
    },
    KnownAdapter {
        name: "mgi-bgi-forward",
        r1: b"AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA",
        r2: None,
    },
    KnownAdapter {
        name: "mgi-bgi-reverse",
        r1: b"AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG",
        r2: None,
    },
    KnownAdapter {
        name: "qiaseq-mirna",
        r1: b"AACTGTAGGCACCATCAAT",
        r2: None,
    },
];

/// Shortest exact seed worth indexing. More permissive matcher settings fall
/// back to the exhaustive matcher rather than filling the index with random
/// one- or two-base keys.
const MIN_INDEXED_SEED: usize = 4;
/// Longest seed representable by [`encode_kmer`].
const MAX_INDEXED_SEED: usize = 16;
/// Common exact matches take the original low-collision fast path before the
/// lossless q-gram fallback is needed.
const FAST_SEED: usize = 12;
/// Minimum column coverage (fraction of voting suffixes) to extend.
const MIN_COVERAGE: f64 = 0.5;
/// Minimum majority fraction to extend with a base.
const MIN_MAJORITY: f64 = 0.6;
/// Largest share of a seed one base may occupy before the seed is treated as a
/// homopolymer artefact rather than adapter evidence.
pub(crate) const MAX_SEED_BASE_SHARE: f64 = 0.75;

/// Encodes `bases` as a 2-bit k-mer; `None` when a non-ACGT base is present.
///
/// Any length up to 16 bases fits a `u32`. The library seed index and the de
/// novo k-mer counter use different lengths, so this is not fixed to one size.
pub(crate) fn encode_kmer(bases: &[u8]) -> Option<u32> {
    debug_assert!(bases.len() <= 16, "a 2-bit k-mer must fit a u32");
    let mut code = 0u32;
    for &base in bases {
        code = (code << 2) | base_slot(base)? as u32;
    }
    Some(code)
}

/// Decodes a 2-bit k-mer back into bases.
pub(crate) fn decode_kmer(mut code: u32, k: usize) -> Vec<u8> {
    let mut out = vec![0; k];
    for slot in out.iter_mut().rev() {
        *slot = b"ACGT"[(code & 3) as usize];
        code >>= 2;
    }
    out
}

/// One exact q-gram that can propose a library entry and its approximate start.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct SeedEntry {
    /// `(length << 32) | two-bit sequence`.
    key: u64,
    index: u32,
    offset: u16,
    /// Maximum coordinate drift caused by edits before this q-gram.
    slack: u16,
}

fn seed_key(bases: &[u8]) -> Option<u64> {
    let length = u64::try_from(bases.len()).ok()?;
    Some((length << 32) | u64::from(encode_kmer(bases)?))
}

/// One q-gram whose position is anchored relative to the read's 3' end.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct TerminalSeed {
    key: u64,
    distance: u16,
    seed: SeedEntry,
}

struct AdapterSeeds {
    full: Vec<SeedEntry>,
    terminal: Vec<TerminalSeed>,
}

/// Splits an adapter prefix into enough pieces that one must survive `errors`.
fn partition_seeds(
    adapter: &Adapter,
    index: u32,
    overlap: usize,
    errors: usize,
) -> Option<Vec<SeedEntry>> {
    if errors >= overlap {
        return None;
    }
    let mut seeds = Vec::new();
    let pieces = errors + 1;
    for piece in 0..pieces {
        let start = piece * overlap / pieces;
        let end = (piece + 1) * overlap / pieces;
        let length = (end - start).min(MAX_INDEXED_SEED);
        if length < MIN_INDEXED_SEED {
            return None;
        }
        let key = seed_key(&adapter.sequence[start..start + length])?;
        seeds.push(SeedEntry {
            key,
            index,
            offset: u16::try_from(start).ok()?,
            slack: u16::try_from(errors).ok()?,
        });
    }
    Some(seeds)
}

/// Lossless exact q-grams for one adapter under `params`.
///
/// Full adapters can occur before the read end, so their longest-prefix seeds
/// are searched everywhere. Partial adapters must end with the read; indexing
/// their seed's distance from that end preserves the same guarantee without
/// making ubiquitous five-base q-grams propose candidates throughout the read.
fn adapter_seeds(adapter: &Adapter, index: usize, params: AdapterParams) -> Option<AdapterSeeds> {
    if adapter.len() < params.min_overlap {
        return None;
    }
    let index = u32::try_from(index).ok()?;
    let mut full = partition_seeds(
        adapter,
        index,
        adapter.len(),
        params.error_limit(adapter.len()),
    )?;
    let mut terminal = Vec::new();
    for overlap in params.min_overlap..adapter.len() {
        let errors = params.error_limit(overlap);
        for seed in partition_seeds(adapter, index, overlap, errors)? {
            let nominal = overlap.checked_sub(seed.offset as usize)?;
            let drift = if params.allow_indels { errors } else { 0 };
            let seed_length = (seed.key >> 32) as usize;
            let first = nominal.saturating_sub(drift).max(seed_length);
            let last = nominal.checked_add(drift)?;
            for distance in first..=last {
                terminal.push(TerminalSeed {
                    key: seed.key,
                    distance: u16::try_from(distance).ok()?,
                    seed,
                });
            }
        }
    }
    full.sort_unstable();
    full.dedup();
    terminal.sort_unstable();
    terminal.dedup();
    Some(AdapterSeeds { full, terminal })
}

fn seed_lengths(keys: impl Iterator<Item = u64>) -> [bool; MAX_INDEXED_SEED + 1] {
    let mut lengths = [false; MAX_INDEXED_SEED + 1];
    for key in keys {
        lengths[(key >> 32) as usize] = true;
    }
    lengths
}

fn terminal_ranges(seeds: &[TerminalSeed]) -> Vec<std::ops::Range<usize>> {
    let Some(maximum) = seeds.last().map(|seed| seed.distance as usize) else {
        return Vec::new();
    };
    let mut ranges = vec![0..0; maximum + 1];
    let mut first = 0usize;
    for (distance, range) in ranges.iter_mut().enumerate() {
        while first < seeds.len() && (seeds[first].distance as usize) < distance {
            first += 1;
        }
        let mut end = first;
        while end < seeds.len() && seeds[end].distance as usize == distance {
            end += 1;
        }
        *range = first..end;
        first = end;
    }
    ranges
}

fn seed_mask(length: usize) -> u64 {
    // `length` never exceeds MAX_INDEXED_SEED, so the shift is at most 32.
    (1u64 << (2 * length)) - 1
}

fn lookup_by_key<T>(items: &[T], key: u64, key_of: impl Fn(&T) -> u64) -> &[T] {
    let first = items.partition_point(|item| key_of(item) < key);
    let rest = &items[first..];
    let count = rest.partition_point(|item| key_of(item) == key);
    &rest[..count]
}

fn lookup_seeds(seeds: &[SeedEntry], key: u64) -> &[SeedEntry] {
    lookup_by_key(seeds, key, |seed| seed.key)
}

fn record_best(best: &mut Option<(usize, AdapterHit)>, index: usize, mut hit: AdapterHit) {
    hit.adapter = index;
    let improves = best.as_ref().is_none_or(|(current_index, current)| {
        (hit.start, hit.errors, std::cmp::Reverse(hit.overlap), index)
            < (
                current.start,
                current.errors,
                std::cmp::Reverse(current.overlap),
                *current_index,
            )
    });
    if improves {
        *best = Some((index, hit));
    }
}

/// Mutable per-read state shared by the fast and lossless seed scans.
struct SeedSearch {
    params: AdapterParams,
    scratch: Option<Vec<u32>>,
    /// Coordinate-local winners used only by the indel matcher's shadow-window
    /// reduction. Substitution mode keeps its allocation-free global winner.
    by_start: Option<Vec<Option<AdapterHit>>>,
    tested: [(usize, usize); 256],
    tested_len: usize,
    best: Option<(usize, AdapterHit)>,
}

impl SeedSearch {
    fn new(params: AdapterParams, read_length: usize) -> Self {
        Self {
            params,
            scratch: params
                .allow_indels
                .then(|| vec![u32::MAX; 2 * (read_length + 1)]),
            by_start: params.allow_indels.then(|| vec![None; read_length + 1]),
            tested: [(usize::MAX, usize::MAX); 256],
            tested_len: 0,
            best: None,
        }
    }

    fn verify(&mut self, adapters: &[Adapter], index: usize, sequence: &[u8], start: usize) {
        let hit = if let Some(scratch) = self.scratch.as_deref_mut() {
            indel_best_at(
                std::slice::from_ref(&adapters[index]),
                self.params,
                sequence,
                start,
                scratch,
            )
        } else {
            verify_at(&adapters[index], self.params, sequence, start)
        };
        if let Some(mut hit) = hit {
            hit.adapter = index;
            if let Some(by_start) = self.by_start.as_mut() {
                let current = &mut by_start[start];
                let improves = current.as_ref().is_none_or(|current| {
                    (hit.errors, std::cmp::Reverse(hit.overlap), index)
                        < (
                            current.errors,
                            std::cmp::Reverse(current.overlap),
                            current.adapter,
                        )
                });
                if improves {
                    *current = Some(hit);
                }
            } else {
                record_best(&mut self.best, index, hit);
            }
        }
    }

    fn hit_at(&self, start: usize) -> Option<(usize, AdapterHit)> {
        if let Some(by_start) = &self.by_start {
            let hit = by_start.get(start).copied().flatten()?;
            Some((hit.adapter, hit))
        } else {
            self.best.filter(|(_, hit)| hit.start == start)
        }
    }

    /// Applies the same final ordering as `find_three_prime`.
    fn result(&self, last_start: usize) -> Option<(usize, AdapterHit)> {
        let Some(by_start) = &self.by_start else {
            return self.best;
        };
        let (first, first_hit) = by_start
            .iter()
            .enumerate()
            .find_map(|(start, hit)| hit.map(|hit| (start, hit)))?;
        let window_end = first.saturating_add(first_hit.errors).min(last_start);
        let mut best = first_hit;
        for &candidate in by_start[first..=window_end].iter().flatten() {
            if indel_improves(&candidate, &best) {
                best = candidate;
            }
        }
        Some((best.adapter, best))
    }

    fn test_seed(
        &mut self,
        adapters: &[Adapter],
        seed: SeedEntry,
        read_offset: usize,
        sequence: &[u8],
    ) {
        let offset = usize::from(seed.offset);
        let slack = if self.params.allow_indels {
            usize::from(seed.slack)
        } else {
            0
        };
        let first = read_offset.saturating_sub(offset.saturating_add(slack));
        let Some(end) = read_offset
            .checked_add(slack)
            .and_then(|position| position.checked_sub(offset))
        else {
            return;
        };
        let end = end.min(sequence.len() - self.params.min_overlap);
        if first > end {
            return;
        }
        for start in first..=end {
            let coordinate = (seed.index as usize, start);
            if self.tested[..self.tested_len].contains(&coordinate) {
                continue;
            }
            if self.tested_len < self.tested.len() {
                self.tested[self.tested_len] = coordinate;
                self.tested_len += 1;
            }
            self.verify(adapters, coordinate.0, sequence, start);
        }
    }
}

/// The known-adapter library, indexed by a lossless set of exact q-grams.
///
/// Fixed 5' seeding is fast but changes matcher behaviour: an error in that seed,
/// a terminal partial shorter than it, or an indel after it can prevent a valid
/// adapter from ever reaching final verification. The q-gram partition above is
/// still cheap to scan, while guaranteeing a proposal whenever the configured
/// matcher can accept the alignment. Pathologically permissive matcher settings
/// fall back to the exhaustive path rather than weakening correctness.
pub(crate) struct SeedIndex {
    adapters: Vec<Adapter>,
    /// Entries that cannot be indexed losslessly under the selected thresholds.
    unseeded: Vec<usize>,
    fast_seeds: Vec<SeedEntry>,
    full_seeds: Vec<SeedEntry>,
    terminal_seeds: Vec<TerminalSeed>,
    full_lengths: [bool; MAX_INDEXED_SEED + 1],
    terminal_lengths: [bool; MAX_INDEXED_SEED + 1],
    terminal_ranges: Vec<std::ops::Range<usize>>,
    /// Furthest a seed may lie beyond a fast hit while still proposing either
    /// an earlier match or a candidate inside its indel shadow window.
    max_scan_reach: usize,
}

impl SeedIndex {
    /// Builds the index for one mate's view of the library.
    pub(crate) fn new(mate: crate::process::Mate, params: AdapterParams) -> Self {
        let adapters: Vec<Adapter> = KNOWN_ADAPTERS
            .iter()
            .map(|entry| Adapter::new(entry.name, entry.sequence(mate)).expect("library is valid"))
            .collect();
        Self::build(adapters, params)
    }

    fn build(adapters: Vec<Adapter>, params: AdapterParams) -> Self {
        let mut full_seeds = Vec::new();
        let mut terminal_seeds = Vec::new();
        let mut fast_seeds = Vec::new();
        let mut unseeded = Vec::new();
        for (index, adapter) in adapters.iter().enumerate() {
            if let Some(key) = adapter.sequence.get(..FAST_SEED).and_then(seed_key) {
                fast_seeds.push(SeedEntry {
                    key,
                    index: index as u32,
                    offset: 0,
                    slack: 0,
                });
            }
            match adapter_seeds(adapter, index, params) {
                Some(seeds) => {
                    full_seeds.extend(seeds.full);
                    terminal_seeds.extend(seeds.terminal);
                }
                None => unseeded.push(index),
            }
        }
        fast_seeds.sort_unstable();
        full_seeds.sort_unstable();
        terminal_seeds.sort_unstable_by_key(|seed| (seed.distance, seed.key, seed.seed));
        let full_lengths = seed_lengths(full_seeds.iter().map(|seed| seed.key));
        let terminal_lengths = seed_lengths(terminal_seeds.iter().map(|seed| seed.key));
        let terminal_ranges = terminal_ranges(&terminal_seeds);
        let seed_reach = full_seeds
            .iter()
            .map(|seed| {
                usize::from(seed.offset)
                    + usize::from(seed.slack) * usize::from(params.allow_indels)
            })
            .chain(terminal_seeds.iter().map(|terminal| {
                usize::from(terminal.seed.offset)
                    + usize::from(terminal.seed.slack) * usize::from(params.allow_indels)
            }))
            .max()
            .unwrap_or(0);
        let indel_shadow = if params.allow_indels {
            adapters
                .iter()
                .map(|adapter| params.error_limit(adapter.len()))
                .max()
                .unwrap_or(0)
        } else {
            0
        };
        Self {
            adapters,
            unseeded,
            fast_seeds,
            full_seeds,
            terminal_seeds,
            full_lengths,
            terminal_lengths,
            terminal_ranges,
            max_scan_reach: seed_reach.saturating_add(indel_shadow),
        }
    }

    /// The library, in declaration order.
    pub(crate) fn adapters(&self) -> &[Adapter] {
        &self.adapters
    }

    /// The library entries worth verifying for one encoded q-gram.
    fn lookup_full(&self, key: u64) -> &[SeedEntry] {
        lookup_seeds(&self.full_seeds, key)
    }

    /// Partial-adapter entries carrying `key` at `distance` from the read end.
    fn lookup_terminal(&self, key: u64, distance: usize) -> &[TerminalSeed] {
        let Some(range) = self.terminal_ranges.get(distance) else {
            return &[];
        };
        let distance_seeds = &self.terminal_seeds[range.clone()];
        lookup_by_key(distance_seeds, key, |seed| seed.key)
    }

    fn scan_fast(&self, sequence: &[u8], search: &mut SeedSearch) -> Option<AdapterHit> {
        if sequence.len() < FAST_SEED {
            return None;
        }
        let mut code = 0u64;
        let mut valid = 0usize;
        for (position, &base) in sequence.iter().enumerate() {
            let Some(slot) = base_slot(base) else {
                valid = 0;
                continue;
            };
            code = ((code << 2) | slot as u64) & seed_mask(FAST_SEED);
            valid += 1;
            if valid < FAST_SEED {
                continue;
            }
            let start = position + 1 - FAST_SEED;
            let key = ((FAST_SEED as u64) << 32) | code;
            for seed in lookup_seeds(&self.fast_seeds, key) {
                search.verify(&self.adapters, seed.index as usize, sequence, start);
            }
            if let Some((_, hit)) = search.hit_at(start) {
                return Some(hit);
            }
        }
        None
    }

    fn scan_indexed(&self, sequence: &[u8], search: &mut SeedSearch, max_read_offset: usize) {
        for length in (MIN_INDEXED_SEED..=MAX_INDEXED_SEED).rev() {
            if (!self.full_lengths[length] && !self.terminal_lengths[length])
                || sequence.len() < length
            {
                continue;
            }
            let mask = seed_mask(length);
            let mut code = 0u64;
            let mut valid = 0usize;
            for (position, &base) in sequence.iter().enumerate() {
                let Some(slot) = base_slot(base) else {
                    valid = 0;
                    continue;
                };
                code = ((code << 2) | slot as u64) & mask;
                valid += 1;
                if valid < length {
                    continue;
                }
                let read_offset = position + 1 - length;
                if read_offset > max_read_offset {
                    break;
                }
                let key = ((length as u64) << 32) | code;
                if self.full_lengths[length] {
                    for &seed in self.lookup_full(key) {
                        search.test_seed(&self.adapters, seed, read_offset, sequence);
                    }
                }
                let distance = sequence.len() - read_offset;
                if self.terminal_lengths[length] && distance < self.terminal_ranges.len() {
                    for terminal in self.lookup_terminal(key, distance) {
                        search.test_seed(&self.adapters, terminal.seed, read_offset, sequence);
                    }
                }
            }
        }
    }

    /// The library entry that best explains `sequence`, and where it starts.
    ///
    /// The earliest qualifying coordinate wins, matching `find_three_prime`; a
    /// read is credited to exactly one entry so that support stays comparable
    /// between entries sharing a 5' prefix, which most of the library does.
    pub(crate) fn find(
        &self,
        sequence: &[u8],
        params: AdapterParams,
    ) -> Option<(usize, crate::adapter::AdapterHit)> {
        if sequence.len() < params.min_overlap {
            return None;
        }
        // A single unindexable entry invalidates candidate-only reduction: its
        // independently finalized hit cannot be merged losslessly with seeded
        // candidates under the indel shadow rule. This rare configuration uses
        // the authoritative matcher over the complete library instead.
        if !self.unseeded.is_empty() {
            return find_three_prime(&self.adapters, params, sequence)
                .map(|hit| (hit.adapter, hit));
        }
        let mut search = SeedSearch::new(params, sequence.len());
        // Exact canonical prefixes dominate real libraries. Preserve the old
        // low-collision fast path. Only a perfect hit at the first coordinate
        // is final; every other hit bounds a lossless indexed rescan.
        let fast_hit = self.scan_fast(sequence, &mut search);
        let last_start = sequence.len() - params.min_overlap;
        if fast_hit.is_some_and(|hit| hit.start == 0 && hit.errors == 0) {
            return search.result(last_start);
        }
        let max_read_offset = fast_hit.map_or(sequence.len(), |hit| {
            hit.start.saturating_add(self.max_scan_reach)
        });
        self.scan_indexed(sequence, &mut search, max_read_offset);
        search.result(last_start)
    }
}

/// Adds one vote for `base` in `columns[column]`; ambiguous bases abstain.
///
/// The column is created before the base is examined. An abstaining `N` must
/// still occupy its column: columns are positional, so skipping one would make
/// the next base land at an index past the end of the vector — and, worse, would
/// silently shift every later vote by one position.
pub(crate) fn vote(columns: &mut Vec<[u64; 4]>, column: usize, base: u8) {
    if columns.len() <= column {
        columns.resize(column + 1, [0; 4]);
    }
    if let Some(slot) = base_slot(base) {
        columns[column][slot] += 1;
    }
}

/// The `ACGT` slot index of a base, or `None` for an ambiguous base.
#[inline]
pub(crate) fn base_slot(base: u8) -> Option<usize> {
    match base {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

/// The majority base of a column, when coverage and majority both hold.
pub(crate) fn majority_base(column: &[u64; 4], voting: u64) -> Option<u8> {
    let votes: u64 = column.iter().sum();
    if voting == 0 || votes == 0 || (votes as f64) < MIN_COVERAGE * voting as f64 {
        return None;
    }
    // A, C, G, T order breaks ties deterministically.
    let (slot, &majority) = column
        .iter()
        .enumerate()
        .max_by_key(|&(slot, &count)| (count, std::cmp::Reverse(slot)))?;
    if (majority as f64) < MIN_MAJORITY * votes as f64 {
        return None;
    }
    Some(b"ACGT"[slot])
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::process::Mate;

    #[test]
    fn the_library_is_well_formed() {
        assert!(!KNOWN_ADAPTERS.is_empty());
        for entry in KNOWN_ADAPTERS {
            assert!(!entry.name.is_empty(), "an entry has no name");
            for mate in [Mate::R1, Mate::R2] {
                let sequence = entry.sequence(mate);
                assert!(!sequence.is_empty(), "{} has an empty mate", entry.name);
                assert!(
                    Adapter::new(entry.name, sequence).is_ok(),
                    "{} is not a usable adapter",
                    entry.name
                );
            }
        }
        // An entry with no distinct R2 uses its R1 sequence on both mates.
        let shared = KNOWN_ADAPTERS
            .iter()
            .find(|entry| entry.r2.is_none())
            .expect("the library has shared-sequence entries");
        assert_eq!(shared.sequence(Mate::R1), shared.sequence(Mate::R2));
    }

    #[test]
    fn kmer_encoding_round_trips() {
        for kmer in [b"ACGTACGTACGT".as_slice(), b"AAAAAAAAAAAA", b"TGCATGCATGCA"] {
            let code = encode_kmer(kmer).expect("encodes");
            assert_eq!(decode_kmer(code, kmer.len()), kmer);
        }
        // Any length up to sixteen bases fits a u32.
        let short = encode_kmer(b"ACGTACGTAC").expect("encodes");
        assert_eq!(decode_kmer(short, 10), b"ACGTACGTAC");
        // An ambiguous base has no encoding.
        assert!(encode_kmer(b"ACGTNCGTACGT").is_none());
    }

    #[test]
    fn the_seed_index_finds_an_adapter_wherever_read_through_starts() {
        let params = AdapterParams::default();
        let index = SeedIndex::new(Mate::R1, params);
        let truseq = b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
        for insert in [0usize, 1, 17, 60] {
            let mut read = vec![b'C'; insert];
            read.extend_from_slice(truseq);
            let (entry, hit) = index
                .find(&read, params)
                .unwrap_or_else(|| panic!("no hit with a {insert}-base insert"));
            assert_eq!(hit.start, insert);
            assert_eq!(index.adapters()[entry].sequence[..13], truseq[..13]);
        }
    }

    #[test]
    fn a_clean_read_matches_nothing() {
        let params = AdapterParams::default();
        let index = SeedIndex::new(Mate::R1, params);
        // Deterministic, non-adapter sequence.
        let read: Vec<u8> = (0..150u8)
            .map(|i| b"ACGT"[(i as usize * 7 + 3) % 4])
            .collect();
        assert!(index.find(&read, params).is_none());
    }

    #[test]
    fn a_seed_lookup_returns_every_entry_sharing_that_qgram() {
        let index = SeedIndex::new(Mate::R1, AdapterParams::default());
        let key = seed_key(b"AGATCGGA").expect("encodes");
        let entries = index.lookup_full(key);
        assert!(
            entries.len() > 1,
            "the library shares this prefix across entries"
        );
        for entry in entries {
            assert_eq!(entry.key, key);
            assert_eq!(entry.offset, 0);
            assert_eq!(
                index.adapters()[entry.index as usize].sequence[..8],
                *b"AGATCGGA"
            );
        }
        // A seed no entry carries returns nothing.
        assert!(index.lookup_full(u64::MAX).is_empty());
    }

    #[test]
    fn an_eight_base_terminal_partial_reaches_verification() {
        let params = AdapterParams::default();
        let index = SeedIndex::new(Mate::R1, params);
        let mut read = vec![b'C'; 37];
        read.extend_from_slice(b"AGATCGGA");
        let (entry, hit) = index.find(&read, params).expect("partial adapter found");
        assert_eq!(hit.start, 37);
        assert_eq!(hit.overlap, 8);
        assert_eq!(index.adapters()[entry].sequence[..8], *b"AGATCGGA");
    }

    #[test]
    fn an_indel_inside_the_old_prefix_seed_reaches_verification() {
        let params = AdapterParams {
            allow_indels: true,
            ..AdapterParams::default()
        };
        let index = SeedIndex::new(Mate::R1, params);
        let adapter = b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
        let mut read = vec![b'C'; 23];
        read.extend_from_slice(&adapter[..6]);
        read.push(b'T');
        read.extend_from_slice(&adapter[6..]);
        let (entry, hit) = index.find(&read, params).expect("indel adapter found");
        assert_eq!(index.adapters()[entry].sequence[..13], adapter[..13]);
        let verified = find_three_prime(
            std::slice::from_ref(&index.adapters()[entry]),
            params,
            &read,
        )
        .expect("proposed adapter verifies");
        assert_eq!(hit, verified);
        assert_eq!(verified.start, 23);
        assert_eq!(verified.errors, 1);
    }

    #[test]
    fn an_earlier_damaged_adapter_beats_a_later_fast_seed_hit() {
        let params = AdapterParams::default();
        let earlier = Adapter::new("earlier", b"GTCAGTACCGATGCTAGCTA").unwrap();
        let later = Adapter::new("later", b"TGCCTAACGGTACCAATGGC").unwrap();
        let adapters = vec![earlier.clone(), later.clone()];
        let index = SeedIndex::build(adapters.clone(), params);

        let mut damaged = earlier.sequence.clone();
        damaged[4] = if damaged[4] == b'A' { b'C' } else { b'A' };
        let read = [
            b"CCCCCCC".as_slice(),
            damaged.as_slice(),
            b"GGGG".as_slice(),
            later.sequence.as_slice(),
        ]
        .concat();

        let exhaustive = find_three_prime(&adapters, params, &read).expect("exhaustive hit");
        let indexed = index.find(&read, params).expect("indexed hit");
        assert_eq!(indexed.0, exhaustive.adapter);
        assert_eq!(indexed.1, exhaustive);
        assert_eq!(exhaustive.start, 7);
        assert_eq!(exhaustive.errors, 1);
    }

    #[test]
    fn an_indel_shadow_after_an_imperfect_origin_hit_is_still_scanned() {
        let params = AdapterParams {
            allow_indels: true,
            ..AdapterParams::default()
        };
        let exact = Adapter::new("exact", b"CGTACCTGATCGGACTTACG").unwrap();
        let read = [b"A".as_slice(), exact.sequence.as_slice()].concat();
        let mut shadow_sequence = read[..20].to_vec();
        shadow_sequence[15] = if shadow_sequence[15] == b'A' {
            b'C'
        } else {
            b'A'
        };
        let shadow = Adapter::new("shadow", &shadow_sequence).unwrap();
        let adapters = vec![shadow, exact];
        let index = SeedIndex::build(adapters.clone(), params);

        let exhaustive = find_three_prime(&adapters, params, &read).expect("exhaustive hit");
        let indexed = index.find(&read, params).expect("indexed hit");
        assert_eq!(indexed.0, exhaustive.adapter);
        assert_eq!(indexed.1, exhaustive);
        assert_eq!(exhaustive.start, 1);
        assert_eq!(exhaustive.errors, 0);
    }

    #[test]
    fn the_index_preserves_exhaustive_matches_across_overlap_boundaries() {
        let adapter = Adapter::new("test", b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA").unwrap();
        for allow_indels in [false, true] {
            let params = AdapterParams {
                allow_indels,
                ..AdapterParams::default()
            };
            let index = SeedIndex::build(vec![adapter.clone()], params);
            let check = |read: &[u8], label: &str| {
                let exhaustive = find_three_prime(std::slice::from_ref(&adapter), params, read);
                let indexed = index.find(read, params).map(|(entry, hit)| {
                    assert_eq!(entry, hit.adapter, "{label}");
                    hit
                });
                assert_eq!(indexed, exhaustive, "{label}");
            };

            for overlap in params.min_overlap..=adapter.len() {
                let prefix = vec![b'C'; 23];
                let exact = [prefix.as_slice(), &adapter.sequence[..overlap]].concat();
                check(&exact, &format!("exact overlap {overlap}"));

                if params.error_limit(overlap) > 0 {
                    let mut substituted = exact.clone();
                    let position = prefix.len() + overlap / 2;
                    substituted[position] = match substituted[position] {
                        b'A' => b'C',
                        _ => b'A',
                    };
                    check(&substituted, &format!("substitution overlap {overlap}"));

                    if allow_indels {
                        let mut inserted = exact.clone();
                        inserted.insert(position, b'T');
                        check(&inserted, &format!("insertion overlap {overlap}"));

                        let mut deleted = exact.clone();
                        deleted.remove(position);
                        check(&deleted, &format!("deletion overlap {overlap}"));
                    }
                }
            }

            // A complete adapter may be followed by more read sequence, so its
            // q-grams cannot rely on terminal anchoring.
            let internal = [
                vec![b'C'; 23],
                adapter.sequence.clone(),
                b"TGCATGC".to_vec(),
            ]
            .concat();
            check(&internal, "complete adapter before read end");
        }
    }

    #[test]
    fn a_column_extends_only_with_coverage_and_a_majority() {
        // Ten voters, nine of them agreeing: extends.
        let strong = [9, 1, 0, 0];
        assert_eq!(majority_base(&strong, 10), Some(b'A'));
        // Ten voters, evenly split: no majority.
        let split = [5, 5, 0, 0];
        assert_eq!(majority_base(&split, 10), None);
        // Four of ten voters cover the column: below the coverage floor.
        let thin = [4, 0, 0, 0];
        assert_eq!(majority_base(&thin, 10), None);
        // Nothing voted at all.
        assert_eq!(majority_base(&[0; 4], 10), None);
        assert_eq!(majority_base(&strong, 0), None);
    }

    #[test]
    fn ambiguous_bases_abstain_without_shifting_a_column() {
        let mut columns = Vec::new();
        vote(&mut columns, 0, b'A');
        vote(&mut columns, 1, b'N');
        vote(&mut columns, 2, b'C');
        // The `N` still occupies its column, so `C` lands at index two.
        assert_eq!(columns.len(), 3);
        assert_eq!(columns[0], [1, 0, 0, 0]);
        assert_eq!(columns[1], [0, 0, 0, 0]);
        assert_eq!(columns[2], [0, 1, 0, 0]);
    }

    #[test]
    fn base_slots_cover_the_canonical_alphabet() {
        assert_eq!(base_slot(b'A'), Some(0));
        assert_eq!(base_slot(b'C'), Some(1));
        assert_eq!(base_slot(b'G'), Some(2));
        assert_eq!(base_slot(b'T'), Some(3));
        assert_eq!(base_slot(b'N'), None);
        assert_eq!(base_slot(b'a'), None);
    }
}
