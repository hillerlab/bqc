// Copyright (c) 2026 The Hiller Lab at the Senckenberg Gesellschaft für Naturforschung
// Distributed under the terms of the GNU General Public License, Version 3.0.

//! `bqc sniff` — non-destructive inspection of a CBQ file.
//!
//! Sniffing never trims, filters, reorders or rewrites anything. It samples a
//! deterministic subset of the input, infers a property from it, and writes a
//! report. The input file is opened read-only and is byte-identical afterwards.
//!
//! Two subcommands share the sampling, execution and report conventions and
//! nothing else, because their evidence has nothing in common:
//!
//! ```text
//! sniff adapters   reference-free, sequence tails      -> candidate adapters
//! sniff strand     needs a transcriptome, mappings     -> library orientation
//! ```
//!
//! Both may return an inconclusive answer, and both say so rather than guessing.
//! `--require-confident` turns an inconclusive answer into a distinct exit code
//! so a pipeline can branch on it without parsing text.

pub mod adapters;
pub mod report;
pub mod sample;
#[cfg(feature = "sniff-strand")]
pub mod strand;

use serde::Serialize;

/// Version of the report envelope. Bumped only for incompatible changes.
pub const SCHEMA_VERSION: u32 = 1;

/// Output projection.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, clap::ValueEnum)]
#[serde(rename_all = "snake_case")]
pub enum Format {
    /// Human-readable summary.
    Text,
    /// Structured report; the stable pipeline interface.
    Json,
    /// One row per candidate (`adapters`) or one summary row (`strand`), for
    /// cohort aggregation.
    Tsv,
}

impl Format {
    #[must_use]
    pub fn name(self) -> &'static str {
        match self {
            Self::Text => "text",
            Self::Json => "json",
            Self::Tsv => "tsv",
        }
    }
}

/// How strongly the evidence supports one candidate.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum Confidence {
    /// Failed a gate that only strong evidence clears.
    Low,
    /// Genuine evidence, short of a recommendation.
    Medium,
    /// Clears every gate; recommendable.
    High,
}

impl Confidence {
    #[must_use]
    pub fn name(self) -> &'static str {
        match self {
            Self::Low => "low",
            Self::Medium => "medium",
            Self::High => "high",
        }
    }
}

/// The mate-level or file-level conclusion.
///
/// Declaration order is severity order: `max` over the mates' decisions is the
/// file-level decision, because a pipeline gating on the whole file cannot be
/// satisfied by one good mate.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum Decision {
    /// Exactly one unrelated candidate is high-confidence.
    Confident,
    /// Nothing clears the high-confidence gates.
    Inconclusive,
    /// Two or more unrelated candidates are high-confidence. Never resolved
    /// automatically: a mixed library is a fact about the data, not a tie.
    Mixed,
}

impl Decision {
    #[must_use]
    pub fn name(self) -> &'static str {
        match self {
            Self::Confident => "confident",
            Self::Mixed => "mixed",
            Self::Inconclusive => "inconclusive",
        }
    }

    /// Whether `--require-confident` is satisfied.
    #[must_use]
    pub fn is_confident(self) -> bool {
        matches!(self, Self::Confident)
    }
}

/// Where a candidate's evidence came from.
///
/// Provenance only. Support is measured once, by the verification pass, because
/// the same read can appear in several evidence sources and summing them would
/// count it more than once.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum EvidenceSource {
    KnownDatabase,
    KmerConsensus,
    PairedOverlap,
}

impl EvidenceSource {
    #[must_use]
    pub fn name(self) -> &'static str {
        match self {
            Self::KnownDatabase => "known_database",
            Self::KmerConsensus => "kmer_consensus",
            Self::PairedOverlap => "paired_overlap",
        }
    }
}

/// The gates a candidate must clear, and the constants behind them.
///
/// These are the algorithmic constants the plan requires to be centralized and
/// serialized with every report, so a result can be re-derived from its own
/// output. They are not exposed as command line options: a user should not have
/// to understand them to read an answer.
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub struct Gates {
    /// Sampled reads below which no candidate is anything but `low`.
    pub min_sample: u64,
    /// Reads that must support a candidate.
    pub min_supporting_reads: u64,
    /// Fraction of sampled reads that must support a candidate.
    pub min_support_fraction: f64,
    /// Shortest candidate that can be recommended.
    pub min_candidate_length: usize,
    /// Largest share of a candidate one base may occupy.
    pub max_base_share: f64,
    /// Minimum adjacent-base complexity.
    pub min_complexity: f64,
    /// Fraction of matches that must reach the read's 3' end.
    pub min_tail_connected_fraction: f64,
    /// Minimum tail-to-body enrichment, normalized by searchable positions.
    pub min_tail_enrichment: f64,
    /// Bases a de novo consensus must extend to before it counts as stable
    /// evidence that a candidate is adapter rather than biology.
    pub min_consensus_length: usize,
    /// Highest mean matcher error rate accepted.
    pub max_mean_error_rate: f64,
    /// A runner-up above this share of the leader's support makes the result
    /// `mixed` rather than `confident`.
    pub competitor_share: f64,
}

impl Default for Gates {
    fn default() -> Self {
        Self {
            // Small enough that a genuinely small file can still be
            // characterised — 200 reads with 40 carrying an adapter is not
            // ambiguous — and large enough that a handful of reads cannot
            // conclude anything.
            min_sample: 100,
            // The absolute floor only guards small samples; on a full 262 144
            // read sample the fraction gate below is the binding one, at 2 621
            // reads. Its job is to stop a couple of chance matches in a tiny
            // file from crossing any fraction.
            min_supporting_reads: 20,
            // The threshold the existing `--detect-min-support` default uses.
            min_support_fraction: 0.01,
            // Shorter than a seed cannot be verified against anything.
            min_candidate_length: 12,
            max_base_share: 0.75,
            min_complexity: 0.3,
            min_tail_connected_fraction: 0.80,
            min_tail_enrichment: 5.0,
            // Twice the seed, so the extension had to survive ten columns of
            // coverage and majority beyond what selected the seed.
            min_consensus_length: 20,
            // `AdapterParams::default().max_error_rate`.
            max_mean_error_rate: 0.10,
            competitor_share: 0.5,
        }
    }
}

/// A share of a whole, with the empty case defined as zero.
#[must_use]
pub fn fraction(part: u64, whole: u64) -> f64 {
    if whole == 0 {
        0.0
    } else {
        part as f64 / whole as f64
    }
}

/// Longest offset-tolerant agreement accepted as "the same adapter".
const RELATED_MAX_ERROR_RATE: f64 = 0.1;
/// Shortest agreement that means anything, regardless of sequence length.
const RELATED_MIN_OVERLAP: usize = 12;

/// Whether two candidate sequences describe the same adapter.
///
/// Only *unrelated* candidates make a library `mixed`, so this decides whether
/// a second high-confidence hit is a second library or another spelling of the
/// first. The bundled database carries both, and they are not all aligned: the
/// same `TruSeq` adapter appears as `AGATCGGAAGAGCACACG...` and as
/// `GATCGGAAGAGCACACG...C`, a one-base frame shift. Comparing from position zero
/// calls those two thirty-three mismatches apart and reports an ordinary
/// single-adapter library as pooled.
///
/// So the sequences are slid past each other and related when *some* offset puts
/// them in agreement over a long enough stretch. This covers prefixes, suffixes,
/// containment and frame shifts. Candidates of nearly equal length get one
/// additional bounded edit-distance comparison, so an internal insertion or
/// deletion accepted by the matcher does not turn one adapter into two apparent
/// libraries.
#[must_use]
pub fn related(a: &[u8], b: &[u8]) -> bool {
    if a.is_empty() || b.is_empty() {
        return true;
    }
    // Never demand more agreement than the shorter sequence can supply, or a
    // sequence stops being related to itself. A de novo consensus is only
    // required to reach `KMER_K` bases, so without the cap two identical
    // ten-base candidates would be proposed separately instead of merging.
    let shortest = a.len().min(b.len());
    let required = RELATED_MIN_OVERLAP.max(shortest / 2).min(shortest);
    let (len_a, len_b) = (a.len().cast_signed(), b.len().cast_signed());
    for offset in -(len_b - 1)..len_a {
        // `b` sits at `offset` relative to `a`; compare where they cover.
        let lo = offset.max(0);
        let hi = (offset + len_b).min(len_a);
        let overlap = (hi - lo).max(0) as usize;
        if overlap < required {
            continue;
        }
        let mismatches = (lo..hi)
            .filter(|&i| a[i as usize] != b[(i - offset) as usize])
            .count();
        if mismatches as f64 <= RELATED_MAX_ERROR_RATE * overlap as f64 {
            return true;
        }
    }
    let edit_limit = (shortest as f64 * RELATED_MAX_ERROR_RATE) as usize;
    a.len().abs_diff(b.len()) <= edit_limit && edit_distance(a, b) <= edit_limit
}

/// Global edit distance for the small candidate set used during discovery.
fn edit_distance(a: &[u8], b: &[u8]) -> usize {
    let mut previous: Vec<usize> = (0..=b.len()).collect();
    let mut current = vec![0; b.len() + 1];
    for (i, &left) in a.iter().enumerate() {
        current[0] = i + 1;
        for (j, &right) in b.iter().enumerate() {
            current[j + 1] = (previous[j + 1] + 1)
                .min(current[j] + 1)
                .min(previous[j] + usize::from(left != right));
        }
        std::mem::swap(&mut previous, &mut current);
    }
    previous[b.len()]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn prefix_variants_are_one_family() {
        assert!(related(
            b"AGATCGGAAGAGC",
            b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
        ));
        // One substitution in thirteen bases is inside the 10% budget the
        // adapter matcher itself uses.
        assert!(related(b"AGATCGGAAGAGC", b"AGATCGGAAGAGT"));
        // Two is not: at 15% divergence these are different sequences, not two
        // spellings of one.
        assert!(!related(b"AGATCGGAAGAGC", b"AGATCGGAAGTGT"));
        // Over a longer alignment the same rate tolerates more.
        assert!(related(
            b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
            b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCT"
        ));
    }

    #[test]
    fn a_frame_shifted_spelling_is_the_same_adapter() {
        // Both of these are in the bundled library, one base out of phase.
        // Comparing from position zero makes them look entirely different, and
        // an ordinary TruSeq library then reports as `mixed`.
        assert!(related(
            b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
            b"GATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
        ));
        // A suffix of one is still the same adapter.
        assert!(related(
            b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
            b"CACACGTCTGAACTCCAGTCA"
        ));
    }

    #[test]
    fn internal_indels_are_one_adapter_family() {
        let canonical = b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
        assert!(related(canonical, b"AGATCGTGAAGAGCACACGTCTGAACTCCAGTCA"));
        assert!(related(canonical, b"AGATCGAAGAGCACACGTCTGAACTCCAGTCA"));
        assert!(!related(canonical, b"AGATCTTTTAGAGCACACGTCTGAACTCCAGTCA"));
    }

    #[test]
    fn a_sequence_is_always_related_to_itself() {
        // Including sequences shorter than the usual agreement floor: a de novo
        // consensus may be only `KMER_K` bases, and two identical ones must
        // merge rather than compete.
        for sequence in [
            b"A".as_slice(),
            b"ACGTACGTAC",
            b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
        ] {
            assert!(related(sequence, sequence), "{sequence:?}");
        }
        assert!(!related(b"A", b"T"));
    }

    #[test]
    fn unrelated_sequences_are_competitors() {
        assert!(!related(b"AGATCGGAAGAGC", b"CTGTCTCTTATACACATCT"));
        assert!(!related(b"AAAAAAAAAAAA", b"TTTTTTTTTTTT"));
        // TruSeq against Nextera: the real "two libraries" case.
        assert!(!related(
            b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
            b"CTGTCTCTTATACACATCTCCGAGCCCACGAGAC"
        ));
    }

    #[test]
    fn confidence_orders_low_to_high() {
        assert!(Confidence::High > Confidence::Medium);
        assert!(Confidence::Medium > Confidence::Low);
    }
}
