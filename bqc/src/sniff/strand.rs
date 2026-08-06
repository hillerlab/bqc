// Copyright (c) 2026 The Hiller Lab at the Senckenberg Gesellschaft für Naturforschung
// Distributed under the terms of the GNU General Public License, Version 3.0.

//! RNA-seq library strandedness.
//!
//! Strandedness cannot be read off raw sequence composition: it is a fact about
//! how reads relate to *oriented transcripts*, so it needs a reference. This
//! maps a deterministic sample of the input through Salmon's selective-alignment
//! mapper and counts the library format each fragment was observed in.
//!
//! Salmon's Rust crates are used directly, against borrowed CBQ slices. There is
//! no FASTQ conversion, no temporary file, no subprocess and no second mapper.
//!
//! # Determinism
//!
//! Salmon's own detector stops after the first N confidently mapped fragments,
//! which makes the answer depend on which thread got there first. Here the
//! sample is split into fixed ordered batches, every record in a batch is mapped
//! before the stopping condition is examined, and the counts are merged by
//! integer addition. The result is identical at every thread count, and early
//! stopping happens only at a batch boundary.
//!
//! # Evidence before inference
//!
//! [`salmon_model::infer_format_from_counts`] is a pure function over counts,
//! and with no usable observations it returns *unstranded* — `U` for single-end,
//! `IU` for paired. That is a sensible default for a quantifier and a terrible
//! answer for a detector, because it is indistinguishable from a real
//! unstranded library. So the evidence gates run first: below them the answer is
//! `undetermined`, and no library type is reported as though it were measured.

use piscem_rs::mapping::hit_searcher::HitSearcher;
use salmon_core::{LibraryFormat, ReadType};
use salmon_core::{ReadOrientation, ReadStrandedness};
use salmon_index::SalmonIndex;
use salmon_map::{MapConfig, ScoredMapping, map_read_pair_into, map_single_read_into};
use serde::Serialize;

use crate::error::{Error, Result};
use crate::io::CbqInput;
use crate::sniff::fraction;
use crate::sniff::sample::{BlockWork, SamplePlan, block_work, for_each_sampled};

/// Number of distinct Salmon library formats.
///
/// `infer_format_from_counts` indexes `0..=MAX_FORMAT_ID` unconditionally, so a
/// shorter slice panics inside the crate. This is that length, named.
const FORMATS: usize = LibraryFormat::MAX_FORMAT_ID as usize + 1;

/// Sampled records mapped before the stopping condition is re-examined.
///
/// Large enough that the barrier costs little against mapping, small enough that
/// a library with a high mapping rate stops early.
const BATCH: u64 = 8_192;

/// Strand-sniffing thresholds.
#[derive(Debug, Clone, Copy, Serialize)]
pub struct Params {
    /// Records examined at most.
    pub sample_size: u64,
    /// Informative observations after which sampling may stop.
    pub target_informative: u64,
    /// Informative observations below which no answer is given.
    pub min_informative: u64,
    /// Informative share of examined records below which no answer is given.
    pub min_informative_fraction: f64,
    /// Forward or reverse share at which a library is called stranded.
    pub stranded_threshold: f64,
    /// Forward/reverse difference below which a library is called unstranded.
    pub unstranded_threshold: f64,
}

impl Default for Params {
    fn default() -> Self {
        Self {
            // nf-core/rnaseq samples a million reads.
            sample_size: 1_000_000,
            // Salmon's own `DEFAULT_SAMPLES_NEEDED`.
            target_informative: 50_000,
            min_informative: 5_000,
            min_informative_fraction: 0.05,
            // nf-core/rnaseq's confident-assignment defaults.
            stranded_threshold: 0.80,
            unstranded_threshold: 0.10,
        }
    }
}

impl Params {
    /// Validates strand thresholds.
    pub fn validate(self) -> Result<Self> {
        if self.sample_size == 0 {
            return Err(Error::config("--sample-size must be at least 1"));
        }
        if self.target_informative == 0 {
            return Err(Error::config("--target-informative must be at least 1"));
        }
        if self.min_informative == 0 {
            return Err(Error::config("--min-informative must be at least 1"));
        }
        if self.target_informative < self.min_informative {
            return Err(Error::config(format!(
                "--target-informative ({}) must be at least --min-informative ({})",
                self.target_informative, self.min_informative
            )));
        }
        if !self.min_informative_fraction.is_finite()
            || !(0.0..=1.0).contains(&self.min_informative_fraction)
        {
            return Err(Error::config(format!(
                "--min-informative-fraction must be within 0.0..=1.0 (got {})",
                self.min_informative_fraction
            )));
        }
        if !self.stranded_threshold.is_finite() || !(0.5..=1.0).contains(&self.stranded_threshold) {
            return Err(Error::config(format!(
                "--stranded-threshold must be within 0.5..=1.0 (got {})",
                self.stranded_threshold
            )));
        }
        if !self.unstranded_threshold.is_finite()
            || !(0.0..0.5).contains(&self.unstranded_threshold)
        {
            return Err(Error::config(format!(
                "--unstranded-threshold must be within 0.0..0.5 (got {})",
                self.unstranded_threshold
            )));
        }
        Ok(self)
    }
}

/// The pipeline-facing answer.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum Strandedness {
    Forward,
    Reverse,
    Unstranded,
    Undetermined,
}

impl Strandedness {
    #[must_use]
    pub fn name(self) -> &'static str {
        match self {
            Self::Forward => "forward",
            Self::Reverse => "reverse",
            Self::Unstranded => "unstranded",
            Self::Undetermined => "undetermined",
        }
    }

    /// `featureCounts -s`. `None` when undetermined: a downstream parameter must
    /// never be manufactured from an answer that was not established.
    #[must_use]
    pub fn featurecounts(self) -> Option<u8> {
        match self {
            Self::Unstranded => Some(0),
            Self::Forward => Some(1),
            Self::Reverse => Some(2),
            Self::Undetermined => None,
        }
    }

    /// `HTSeq` `--stranded`.
    #[must_use]
    pub fn htseq(self) -> Option<&'static str> {
        match self {
            Self::Unstranded => Some("no"),
            Self::Forward => Some("yes"),
            Self::Reverse => Some("reverse"),
            Self::Undetermined => None,
        }
    }
}

/// Dominant relative orientation of the mates.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum Orientation {
    Inward,
    Outward,
    Matching,
    Undetermined,
}

impl Orientation {
    #[must_use]
    pub fn name(self) -> &'static str {
        match self {
            Self::Inward => "inward",
            Self::Outward => "outward",
            Self::Matching => "matching",
            Self::Undetermined => "undetermined",
        }
    }
}

/// Provenance of the index a result was measured against.
#[derive(Debug, Clone, Serialize)]
pub struct IndexMetadata {
    pub path: String,
    pub salmon_version: String,
    pub index_version: u32,
    pub num_refs: usize,
    pub num_decoys: usize,
    pub has_decoys: bool,
    pub k: usize,
    /// Digest of the reference sequences, for matching a result to an index.
    pub seq_hash: String,
    pub name_hash: String,
}

impl IndexMetadata {
    fn new(index: &SalmonIndex, path: &std::path::Path) -> Self {
        let info = index.info();
        Self {
            path: path.display().to_string(),
            salmon_version: info.salmon_version.clone(),
            index_version: info.index_version,
            num_refs: info.num_refs,
            num_decoys: info.num_decoys,
            has_decoys: info.first_decoy_index.is_some(),
            k: info.k,
            seq_hash: info.seq_hash.clone(),
            name_hash: info.name_hash.clone(),
        }
    }
}

/// One observed library format and how often it was seen.
#[derive(Debug, Clone, Serialize)]
pub struct FormatCount {
    pub format: &'static str,
    pub count: u64,
}

/// The whole `sniff strand` result.
#[derive(Debug, Clone, Serialize)]
pub struct StrandSniff {
    pub decision: &'static str,
    /// Why no recommendation was made, as the JSON reason string.
    pub failure_reason: Option<&'static str>,
    /// Salmon's canonical library type, e.g. `ISR`. Absent when undetermined.
    pub salmon_library_type: Option<&'static str>,
    pub strandedness: Strandedness,
    pub pair_orientation: Orientation,
    pub forward_count: u64,
    pub reverse_count: u64,
    pub forward_fraction: f64,
    pub reverse_fraction: f64,
    pub format_counts: Vec<FormatCount>,
    pub records_examined: u64,
    pub mapped_records: u64,
    pub informative_records: u64,
    pub ambiguous_records: u64,
    pub decoy_only_records: u64,
    pub unmapped_records: u64,
    pub mapping_errors: u64,
    pub informative_fraction: f64,
    pub featurecounts_strand: Option<u8>,
    pub htseq_stranded: Option<&'static str>,
    pub index_metadata: IndexMetadata,
}

impl StrandSniff {
    /// Whether `--require-confident` is satisfied.
    #[must_use]
    pub fn is_confident(&self) -> bool {
        self.strandedness != Strandedness::Undetermined
    }
}

/// Per-worker mapping tallies. Merged by addition, so the merge is independent
/// of which worker saw which batch.
#[derive(Clone)]
struct Tally {
    formats: [u64; FORMATS],
    examined: u64,
    mapped: u64,
    informative: u64,
    ambiguous: u64,
    decoy_only: u64,
    errors: u64,
}

impl Tally {
    fn new() -> Self {
        Self {
            formats: [0; FORMATS],
            examined: 0,
            mapped: 0,
            informative: 0,
            ambiguous: 0,
            decoy_only: 0,
            errors: 0,
        }
    }

    fn merge(&mut self, other: &Self) {
        for (slot, add) in self.formats.iter_mut().zip(&other.formats) {
            *slot += add;
        }
        self.examined += other.examined;
        self.mapped += other.mapped;
        self.informative += other.informative;
        self.ambiguous += other.ambiguous;
        self.decoy_only += other.decoy_only;
        self.errors += other.errors;
    }

    /// Records one fragment's mappings.
    ///
    /// A fragment counts only when the mappings that survived agree on what they
    /// saw. Decoy filtering already happened inside the mapper, so an empty
    /// result means either nothing mapped or the best placement was a decoy —
    /// distinguished by the mapper's own per-fragment statistics.
    fn record(&mut self, mappings: &[ScoredMapping]) {
        self.examined += 1;
        let stats = salmon_map::take_last_map_stats();
        if mappings.is_empty() {
            if stats.decoy_dominated {
                self.decoy_only += 1;
            }
            return;
        }
        self.mapped += 1;

        // Orphans carry no observed format; only a fragment whose surviving
        // mappings all agree is evidence.
        let mut observed: Option<LibraryFormat> = None;
        for mapping in mappings {
            let Some(format) = mapping.format else {
                self.ambiguous += 1;
                return;
            };
            match observed {
                None => observed = Some(format),
                Some(seen) if seen == format => {}
                Some(_) => {
                    self.ambiguous += 1;
                    return;
                }
            }
        }
        if let Some(format) = observed {
            self.formats[format.format_id() as usize] += 1;
            self.informative += 1;
        } else {
            self.ambiguous += 1;
        }
    }
}

/// Loads the index, refusing anything the mapper cannot read with an actionable
/// message.
///
/// `SalmonIndex::load` rejects C++ salmon (pufferfish) indexes and Rust salmon
/// indexes below format version 1. Those are the overwhelming majority of
/// indexes in existence, so the failure must say what to do rather than surface
/// an internal error.
fn load_index(path: &std::path::Path) -> Result<SalmonIndex> {
    if !path.is_dir() {
        return Err(Error::config(format!(
            "--index {} is not a directory; it must be a Salmon index built by \
             `salmon index -t <transcripts> -i <index>`",
            path.display()
        )));
    }
    SalmonIndex::load(path).map_err(|error| {
        Error::config(format!(
            "cannot read the Salmon index at {}: {error}\n\
             `bqc sniff strand` uses the Salmon 2.x index format. An index \
             built by salmon 1.x cannot be read and must be rebuilt with \
             `salmon index -t <transcripts> -i <index>`",
            path.display()
        ))
    })
}

/// Builds a Salmon index from a transcriptome FASTA into a fresh temp dir and
/// returns that dir's path, ready for [`load_index`]. The index stays alive for
/// the process and OS temp cleanup removes it afterwards.
// : no persistent index cache — rebuilt every run; cache on
// transcriptome path+mtime when build time starts to matter.
pub fn build_temp_index(
    transcriptome: &std::path::Path,
    threads: usize,
) -> Result<std::path::PathBuf> {
    if !transcriptome.is_file() {
        return Err(Error::config(format!(
            "--transcriptome {} is not a file; give a transcriptome FASTA",
            transcriptome.display()
        )));
    }
    let output = std::env::temp_dir().join(format!(
        "bqc-strand-{}-{:x}",
        std::process::id(),
        now_unique()
    ));
    let mut options =
        salmon_index::IndexBuildOptions::new(vec![transcriptome.to_path_buf()], output.clone());
    options.threads = threads;
    salmon_index::build(&options).map_err(|error| {
        Error::config(format!(
            "cannot build a Salmon index from {}: {error}",
            transcriptome.display()
        ))
    })?;
    Ok(output)
}

fn now_unique() -> u64 {
    use std::time::{SystemTime, UNIX_EPOCH};
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_or(0, |d| d.as_nanos() as u64)
}

/// Maps a deterministic sample and infers the library type from what it saw.
pub fn sniff(
    input: &CbqInput,
    index_path: &std::path::Path,
    params: Params,
    span: Option<&std::ops::Range<u64>>,
    threads: usize,
) -> Result<(StrandSniff, SamplePlan)> {
    let index = load_index(index_path)?;
    let paired = input.schema().paired;
    let read_type = if paired {
        ReadType::PairedEnd
    } else {
        ReadType::SingleEnd
    };
    let plan = SamplePlan::for_input(input, span, params.sample_size);
    let work = block_work(&plan, input);
    let config = MapConfig::default();
    let minimum = index.k();

    // Fixed ordered batches. Every record of a batch is mapped before the
    // stopping condition is examined, so early stopping cannot depend on which
    // worker finished first.
    let mut total = Tally::new();
    let mut batch_start = 0u64;
    while batch_start < plan.selected {
        let batch_end = (batch_start + BATCH).min(plan.selected);
        let slice: Vec<BlockWork> = work
            .iter()
            .filter_map(|item| clamp(*item, batch_start, batch_end))
            .collect();
        let batch = map_batch(input, &slice, &plan, threads, &index, &config, minimum)?;
        total.merge(&batch);
        batch_start = batch_end;
        if total.informative >= params.target_informative {
            break;
        }
    }

    Ok((
        conclude(&total, read_type, params, &index, index_path),
        plan,
    ))
}

/// Restricts one block's work to a batch's ordinal range.
fn clamp(item: BlockWork, first: u64, end: u64) -> Option<BlockWork> {
    let clamped = BlockWork {
        block: item.block,
        first_ordinal: item.first_ordinal.max(first),
        end_ordinal: item.end_ordinal.min(end),
    };
    (clamped.first_ordinal < clamped.end_ordinal).then_some(clamped)
}

/// Maps one batch in parallel.
///
/// The index is shared by reference — it is `Sync` — while each worker owns its
/// `HitSearcher` and reuses one output vector, so mapping allocates nothing per
/// fragment.
fn map_batch(
    input: &CbqInput,
    work: &[BlockWork],
    plan: &SamplePlan,
    threads: usize,
    index: &SalmonIndex,
    config: &MapConfig,
    minimum: usize,
) -> Result<Tally> {
    let states = for_each_sampled(
        input,
        work,
        plan,
        threads,
        || (Tally::new(), HitSearcher::new(index.inner()), Vec::new()),
        |(tally, searcher, mappings): &mut (Tally, HitSearcher<'_>, Vec<ScoredMapping>), r1, r2| {
            map_one(tally, mappings, searcher, index, config, minimum, r1, r2);
        },
    )?;
    let mut merged = Tally::new();
    for (tally, _, _) in &states {
        merged.merge(tally);
    }
    Ok(merged)
}

/// Maps one fragment and records what it showed.
#[allow(clippy::too_many_arguments)]
fn map_one<'idx>(
    tally: &mut Tally,
    mappings: &mut Vec<ScoredMapping>,
    searcher: &mut HitSearcher<'idx>,
    index: &'idx SalmonIndex,
    config: &MapConfig,
    minimum: usize,
    r1: &[u8],
    r2: Option<&[u8]>,
) {
    // The mapper has no error channel — it returns mappings, not a result — so
    // "mapping errors" can only mean reads it was never given: a read shorter
    // than the index k-mer has no seed and cannot be mapped at all.
    let too_short = r1.len() < minimum || r2.is_some_and(|mate| mate.len() < minimum);
    if too_short {
        tally.examined += 1;
        tally.errors += 1;
        return;
    }
    match r2 {
        Some(r2) => map_read_pair_into(mappings, index.inner(), searcher, index, r1, r2, config),
        None => map_single_read_into(mappings, index.inner(), searcher, index, r1, config),
    }
    tally.record(mappings);
}

/// Applies the evidence gates and turns counts into an answer.
fn conclude(
    tally: &Tally,
    read_type: ReadType,
    params: Params,
    index: &SalmonIndex,
    index_path: &std::path::Path,
) -> StrandSniff {
    let informative_fraction = fraction(tally.informative, tally.examined);
    let enough = enough_evidence(tally, params);

    // Strand evidence, counted over the formats that carry a strand. Salmon's
    // `MSF`/`MSR` use `S`/`A` where the inward and outward formats use `SA`/`AS`,
    // so the strandedness field is the thing to read, not the format id.
    let (mut forward, mut reverse) = (0u64, 0u64);
    let (mut inward, mut outward, mut matching) = (0u64, 0u64, 0u64);
    let mut format_counts = Vec::new();
    for id in 0..=LibraryFormat::MAX_FORMAT_ID {
        let count = tally.formats[id as usize];
        let format = LibraryFormat::from_format_id(id);
        if count > 0 {
            format_counts.push(FormatCount {
                format: format.canonical(),
                count,
            });
        }
        match format.strandedness {
            ReadStrandedness::SA | ReadStrandedness::S => forward += count,
            ReadStrandedness::AS | ReadStrandedness::A => reverse += count,
            ReadStrandedness::U => {}
        }
        match format.orientation {
            ReadOrientation::Toward => inward += count,
            ReadOrientation::Away => outward += count,
            ReadOrientation::Same => matching += count,
            ReadOrientation::None => {}
        }
    }

    let stranded_total = forward + reverse;
    let forward_fraction = fraction(forward, stranded_total);
    let reverse_fraction = fraction(reverse, stranded_total);

    let strandedness = if !enough {
        Strandedness::Undetermined
    } else if forward_fraction >= params.stranded_threshold {
        Strandedness::Forward
    } else if reverse_fraction >= params.stranded_threshold {
        Strandedness::Reverse
    } else if (forward_fraction - reverse_fraction).abs() < params.unstranded_threshold {
        Strandedness::Unstranded
    } else {
        Strandedness::Undetermined
    };

    // Salmon's inference runs only once the gates have passed, so its
    // no-evidence fallback can never be mistaken for a measurement.
    let salmon_library_type = enough
        .then(|| salmon_model::infer_format_from_counts(&tally.formats, read_type).canonical());

    let pair_orientation = match read_type {
        ReadType::SingleEnd => Orientation::Undetermined,
        ReadType::PairedEnd if !enough => Orientation::Undetermined,
        ReadType::PairedEnd => dominant(inward, outward, matching),
    };

    StrandSniff {
        decision: if strandedness == Strandedness::Undetermined {
            "undetermined"
        } else {
            "confident"
        },
        failure_reason: (!enough).then_some("insufficient_mapping_evidence"),
        salmon_library_type,
        strandedness,
        pair_orientation,
        forward_count: forward,
        reverse_count: reverse,
        forward_fraction,
        reverse_fraction,
        format_counts,
        records_examined: tally.examined,
        mapped_records: tally.mapped,
        informative_records: tally.informative,
        ambiguous_records: tally.ambiguous,
        decoy_only_records: tally.decoy_only,
        unmapped_records: tally.examined - tally.mapped - tally.errors,
        mapping_errors: tally.errors,
        informative_fraction,
        featurecounts_strand: strandedness.featurecounts(),
        htseq_stranded: strandedness.htseq(),
        index_metadata: IndexMetadata::new(index, index_path),
    }
}

/// Evidence gates with an unconditional no-observation guard.
///
/// Validation prevents zero thresholds at the CLI, but keeping the invariant
/// here protects library callers and future configuration paths as well.
fn enough_evidence(tally: &Tally, params: Params) -> bool {
    tally.informative > 0
        && tally.informative >= params.min_informative
        && fraction(tally.informative, tally.examined) >= params.min_informative_fraction
}

/// The dominant orientation, or `Undetermined` when nothing was observed.
///
/// An unusual orientation is reported as it was measured. Collapsing outward or
/// matching into inward would hide exactly the finding worth surfacing.
fn dominant(inward: u64, outward: u64, matching: u64) -> Orientation {
    if inward == 0 && outward == 0 && matching == 0 {
        return Orientation::Undetermined;
    }
    if inward >= outward && inward >= matching {
        Orientation::Inward
    } else if outward >= inward && outward >= matching {
        Orientation::Outward
    } else {
        Orientation::Matching
    }
}

#[cfg(test)]
#[allow(clippy::float_cmp)]
mod tests {
    use super::*;

    fn tally(formats: &[(&str, u64)], examined: u64) -> Tally {
        let mut tally = Tally::new();
        tally.examined = examined;
        for (name, count) in formats {
            let format = LibraryFormat::parse(name).expect("a valid library format");
            tally.formats[format.format_id() as usize] = *count;
            tally.informative += count;
            tally.mapped += count;
        }
        tally
    }

    fn params() -> Params {
        Params {
            min_informative: 100,
            min_informative_fraction: 0.05,
            ..Params::default()
        }
    }

    #[test]
    fn the_counts_array_is_exactly_as_long_as_salmon_indexes_it() {
        // `infer_format_from_counts` reads `0..=MAX_FORMAT_ID` unconditionally.
        assert_eq!(FORMATS, 12);
        let counts = [0u64; FORMATS];
        // Would panic on a shorter slice.
        let format = salmon_model::infer_format_from_counts(&counts, ReadType::PairedEnd);
        assert_eq!(format.canonical(), "IU");
    }

    #[test]
    fn strand_thresholds_follow_the_nf_core_defaults() {
        let defaults = Params::default();
        assert_eq!(defaults.stranded_threshold, 0.80);
        assert_eq!(defaults.unstranded_threshold, 0.10);
        assert_eq!(defaults.target_informative, 50_000);
        assert!(defaults.validate().is_ok());
    }

    #[test]
    fn out_of_range_thresholds_are_rejected() {
        for stranded in [0.4, 1.1, f64::NAN, f64::INFINITY] {
            let params = Params {
                stranded_threshold: stranded,
                ..Params::default()
            };
            assert!(params.validate().is_err(), "{stranded} was accepted");
        }
        for unstranded in [-0.1, 0.5, 0.9, f64::NAN, f64::INFINITY] {
            let params = Params {
                unstranded_threshold: unstranded,
                ..Params::default()
            };
            assert!(params.validate().is_err(), "{unstranded} was accepted");
        }
        for fraction in [-0.1, 1.1, f64::NAN, f64::INFINITY] {
            let params = Params {
                min_informative_fraction: fraction,
                ..Params::default()
            };
            assert!(params.validate().is_err(), "{fraction} was accepted");
        }
    }

    #[test]
    fn evidence_count_gates_are_coherent() {
        assert!(
            Params {
                target_informative: 0,
                ..Params::default()
            }
            .validate()
            .is_err()
        );
        assert!(
            Params {
                min_informative: 0,
                ..Params::default()
            }
            .validate()
            .is_err()
        );
        assert!(
            Params {
                target_informative: 99,
                min_informative: 100,
                ..Params::default()
            }
            .validate()
            .is_err()
        );
    }

    #[test]
    fn a_forward_library_is_called_forward() {
        let counts = tally(&[("ISF", 900), ("ISR", 100)], 2000);
        let (forward, reverse) = strand_counts(&counts);
        assert_eq!((forward, reverse), (900, 100));
    }

    #[test]
    fn matching_orientation_formats_still_carry_a_strand() {
        // MSF/MSR use `S`/`A` rather than `SA`/`AS`; reading the format id
        // instead of the strandedness field would miss them.
        let counts = tally(&[("MSF", 700), ("MSR", 300)], 2000);
        let (forward, reverse) = strand_counts(&counts);
        assert_eq!((forward, reverse), (700, 300));
    }

    /// Strand totals, extracted the way `conclude` does.
    fn strand_counts(tally: &Tally) -> (u64, u64) {
        use salmon_core::ReadStrandedness;
        let (mut forward, mut reverse) = (0, 0);
        for id in 0..=LibraryFormat::MAX_FORMAT_ID {
            let count = tally.formats[id as usize];
            match LibraryFormat::from_format_id(id).strandedness {
                ReadStrandedness::SA | ReadStrandedness::S => forward += count,
                ReadStrandedness::AS | ReadStrandedness::A => reverse += count,
                ReadStrandedness::U => {}
            }
        }
        (forward, reverse)
    }

    #[test]
    fn no_evidence_is_undetermined_not_unstranded() {
        // Salmon's pure inference returns `IU` from an all-zero array. A
        // detector must not present that as a measurement.
        let empty = Tally::new();
        let informative_fraction = fraction(empty.informative, empty.examined);
        assert_eq!(informative_fraction, 0.0);
        assert!(empty.informative < params().min_informative);

        // The core invariant does not rely only on CLI validation: even an
        // unvalidated all-zero gate cannot turn no observations into evidence.
        let zero_gates = Params {
            target_informative: 0,
            min_informative: 0,
            min_informative_fraction: 0.0,
            ..Params::default()
        };
        assert!(!enough_evidence(&empty, zero_gates));
    }

    #[test]
    fn the_dominant_orientation_is_reported_as_measured() {
        assert_eq!(dominant(900, 50, 50), Orientation::Inward);
        assert_eq!(dominant(50, 900, 50), Orientation::Outward);
        assert_eq!(dominant(50, 50, 900), Orientation::Matching);
        assert_eq!(dominant(0, 0, 0), Orientation::Undetermined);
        // Ties resolve deterministically, inward first.
        assert_eq!(dominant(100, 100, 100), Orientation::Inward);
    }

    #[test]
    fn downstream_parameters_are_absent_when_undetermined() {
        assert_eq!(Strandedness::Unstranded.featurecounts(), Some(0));
        assert_eq!(Strandedness::Forward.featurecounts(), Some(1));
        assert_eq!(Strandedness::Reverse.featurecounts(), Some(2));
        assert_eq!(Strandedness::Undetermined.featurecounts(), None);
        assert_eq!(Strandedness::Reverse.htseq(), Some("reverse"));
        assert_eq!(Strandedness::Undetermined.htseq(), None);
    }

    #[test]
    fn tallies_merge_by_addition() {
        let mut a = tally(&[("ISF", 10)], 20);
        let b = tally(&[("ISF", 5), ("ISR", 3)], 12);
        a.merge(&b);
        let format = LibraryFormat::parse("ISF").unwrap();
        assert_eq!(a.formats[format.format_id() as usize], 15);
        assert_eq!(a.examined, 32);
        assert_eq!(a.informative, 18);
    }

    #[test]
    fn a_batch_clamps_to_the_blocks_it_covers() {
        let item = BlockWork {
            block: 3,
            first_ordinal: 100,
            end_ordinal: 200,
        };
        // Wholly inside.
        assert_eq!(clamp(item, 0, 1000), Some(item));
        // Partly overlapping.
        assert_eq!(
            clamp(item, 150, 1000),
            Some(BlockWork {
                block: 3,
                first_ordinal: 150,
                end_ordinal: 200
            })
        );
        // Disjoint.
        assert_eq!(clamp(item, 200, 300), None);
        assert_eq!(clamp(item, 0, 100), None);
    }
}
