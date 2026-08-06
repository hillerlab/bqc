// Copyright (c) 2026 The Hiller Lab at the Senckenberg Gesellschaft für Naturforschung
// Distributed under the terms of the GNU General Public License, Version 3.0.

//! Command line front end.
//!
//! Each command is converted into a [`Config`], which is the only
//! configuration type the library understands. The shared option enums derive
//! `clap::ValueEnum` directly, so the CLI needs no mirror types for them.

use std::collections::HashSet;
use std::ops::Range;
use std::path::{Path, PathBuf};

use clap::{Args, Parser, Subcommand, ValueEnum};
#[cfg(feature = "sniff-strand")]
use clap::ArgGroup;

use crate::config::{
    resolve_threads, AdapterOptions, Config, CorrectionOptions, FilterOptions, LinkedOptions,
    Merge, OutputOptions, Plan, PolyOptions, QualityCutOptions, ReportFormat, SegmentOptions, Step,
    TrimOptions,
};
use crate::correct::LogDetail;
use crate::engine::{Outputs, RunOptions};
use crate::error::{Error, Result};
use crate::io::{CbqInput, CbqOutput, TextOutput};
use crate::linked::{Require, Unmatched};
use crate::process::{FailedMode, PairPolicy};
use crate::report::{Report, RunReport};
use crate::segment::Terminal;
use crate::sniff::Format;

/// CBQ-native adapter removal, read trimming and per-read filtering.
#[derive(Debug, Parser)]
#[command(name = "bqc", version, about, long_about = None)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Command,
}

#[derive(Debug, Subcommand)]
pub enum Command {
    /// Remove adapter sequences from the 3' end of each read.
    Adapter(Box<AdapterCommand>),
    /// Shorten reads by position, quality, terminal Ns or homopolymer tails.
    Trim(Box<TrimCommand>),
    /// Accept or reject reads and pairs.
    Filter(Box<FilterCommand>),
    /// Correct low-quality bases from the other mate where the pair overlaps.
    Correct(Box<CorrectCommand>),
    /// Split reads at internal adapter occurrences into separate records.
    Segment(Box<SegmentCommand>),
    /// Run any combination of the adapter, trim and filter stages in one pass.
    Workflow(Box<WorkflowCommand>),
    /// Inspect a file without modifying it.
    Sniff(Box<SniffCommand>),
}

/// Non-destructive inspection.
///
/// These are subcommands rather than flags on one command because their inputs
/// differ: adapter discovery is reference-free and works on any library, while
/// strandedness needs a transcriptome (or its index) and only means anything
/// for RNA.
#[derive(Debug, Args)]
pub struct SniffCommand {
    #[command(subcommand)]
    pub kind: SniffKind,
}

#[derive(Debug, Subcommand)]
pub enum SniffKind {
    /// Infer which adapter sequences contaminate the reads.
    Adapters(Box<SniffAdaptersCommand>),
    /// Infer RNA-seq library strandedness against a Salmon transcriptome index.
    #[cfg(feature = "sniff-strand")]
    Strand(Box<SniffStrandCommand>),
}

#[cfg(feature = "sniff-strand")]
#[derive(Debug, Args)]
#[command(group(
    ArgGroup::new("index-source")
        .required(true)
        .multiple(false)
        .args(["index", "transcriptome"])
))]
pub struct SniffStrandCommand {
    #[command(flatten)]
    pub common: SniffArgs,

    /// Salmon transcriptome index directory; mutually exclusive with
    /// --transcriptome.
    #[arg(long, value_name = "PATH")]
    pub index: Option<PathBuf>,

    /// Transcriptome FASTA; a Salmon index is built from it on the fly.
    /// Mutually exclusive with --index.
    #[arg(long, value_name = "PATH")]
    pub transcriptome: Option<PathBuf>,

    /// Informative observations after which sampling may stop [default: 50000].
    #[arg(long, value_name = "INT")]
    pub target_informative: Option<u64>,

    /// Informative observations below which no answer is given [default: 5000].
    #[arg(long, value_name = "INT")]
    pub min_informative: Option<u64>,

    /// Informative share of examined records required [default: 0.05].
    #[arg(long, value_name = "FLOAT")]
    pub min_informative_fraction: Option<f64>,

    /// Forward or reverse share at which a library is called stranded
    /// [default: 0.80].
    #[arg(long, value_name = "FLOAT")]
    pub stranded_threshold: Option<f64>,

    /// Forward/reverse difference below which a library is called unstranded
    /// [default: 0.10].
    #[arg(long, value_name = "FLOAT")]
    pub unstranded_threshold: Option<f64>,
}

/// Options shared by every `sniff` subcommand.
#[derive(Debug, Args)]
pub struct SniffArgs {
    /// Input CBQ file. Never modified.
    #[arg(value_name = "INPUT.cbq")]
    pub input: PathBuf,

    /// Records sampled across the input, spread evenly rather than taken from
    /// the start.
    #[arg(long, value_name = "INT")]
    pub sample_size: Option<u64>,

    /// Restrict sampling to original record indices START..END.
    #[arg(long, value_name = "START..END")]
    pub span: Option<String>,

    /// Worker threads; 0 uses every available core.
    #[arg(short = 'T', long, value_name = "INT")]
    pub threads: Option<usize>,

    /// Write the report here instead of to stdout.
    #[arg(short = 'o', long, value_name = "PATH")]
    pub output: Option<PathBuf>,

    /// Report projection.
    #[arg(long, value_enum, value_name = "FORMAT", default_value = "text")]
    pub format: Format,

    /// Exit with status 3 when the result is not confident.
    #[arg(long)]
    pub require_confident: bool,

    /// Overwrite an existing output file.
    #[arg(long)]
    pub force: bool,
}

#[derive(Debug, Args)]
pub struct SniffAdaptersCommand {
    #[command(flatten)]
    pub common: SniffArgs,

    /// Candidates reported per mate.
    #[arg(long, value_name = "INT")]
    pub top: Option<usize>,

    /// Write a `bqc` TOML fragment configuring the detected adapters.
    /// Only written when the result is uniquely confident.
    #[arg(long, value_name = "PATH")]
    pub emit_config: Option<PathBuf>,
}

/// How a command finished.
///
/// Separate from [`Result`] because "the answer is not confident" is a
/// successful analysis, not a failure: it becomes a distinct exit status only
/// when the caller asked for one with `--require-confident`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Outcome {
    Success,
    NotConfident,
}

/// Segmentation turns one record into many, so it is its own command rather than
/// a step of `workflow`: no other stage can assume one record in, one record out.
#[derive(Debug, Args)]
pub struct SegmentCommand {
    #[command(flatten)]
    pub common: CommonArgs,
    #[command(flatten)]
    pub segment: SegmentArgs,
    #[command(flatten)]
    pub trim: TrimArgs,
    #[command(flatten)]
    pub filter: FilterArgs,
}

#[derive(Debug, Args)]
pub struct AdapterCommand {
    #[command(flatten)]
    pub common: CommonArgs,
    #[command(flatten)]
    pub adapter: AdapterArgs,
}

#[derive(Debug, Args)]
pub struct TrimCommand {
    #[command(flatten)]
    pub common: CommonArgs,
    #[command(flatten)]
    pub trim: TrimArgs,
}

#[derive(Debug, Args)]
pub struct FilterCommand {
    #[command(flatten)]
    pub common: CommonArgs,
    #[command(flatten)]
    pub filter: FilterArgs,
}

#[derive(Debug, Args)]
pub struct CorrectCommand {
    #[command(flatten)]
    pub common: CommonArgs,
    #[command(flatten)]
    pub correction: CorrectionArgs,
    /// Minimum R1/R2 overlap accepted as evidence [default: 30].
    #[arg(long, value_name = "INT")]
    pub paired_overlap_min_overlap: Option<usize>,
    /// Maximum mismatch fraction over the overlap [default: 0.10].
    #[arg(long, value_name = "FLOAT")]
    pub max_error_rate: Option<f64>,
}

#[derive(Debug, Args)]
pub struct WorkflowCommand {
    #[command(flatten)]
    pub common: CommonArgs,
    #[command(flatten)]
    pub adapter: AdapterArgs,
    #[command(flatten)]
    pub trim: TrimArgs,
    #[command(flatten)]
    pub filter: FilterArgs,
    #[command(flatten)]
    pub correction_args: CorrectionArgs,

    /// Stages to run (canonical form, e.g. `correct,adapter,trim,filter`).
    #[arg(long, value_delimiter = ',', value_name = "LIST")]
    pub steps: Option<Vec<StepArg>>,
    /// Skip the adapter stage.
    #[arg(long, conflicts_with = "steps")]
    pub no_adapter: bool,
    /// Skip the trim stage.
    #[arg(long, conflicts_with = "steps")]
    pub no_trim: bool,
    /// Skip the filter stage.
    #[arg(long, conflicts_with = "steps")]
    pub no_filter: bool,
    /// Enable paired-overlap base correction.
    #[arg(long)]
    pub correction: bool,
    /// TOML configuration file; command line arguments take precedence.
    #[arg(long, value_name = "PATH")]
    pub config: Option<PathBuf>,
}

/// Options shared by every command.
#[derive(Debug, Args)]
pub struct CommonArgs {
    /// Input CBQ file.
    #[arg(value_name = "INPUT.cbq")]
    pub input: PathBuf,

    /// Accepted output CBQ file.
    #[arg(short = 'o', long, value_name = "PATH")]
    pub output: PathBuf,

    /// Worker threads; 0 uses every available core.
    #[arg(short = 'T', long, value_name = "INT")]
    pub threads: Option<usize>,

    /// Restrict processing to original record indices START..END.
    #[arg(long, value_name = "START..END")]
    pub span: Option<String>,

    /// Output zstd compression level (default: inherited from the input).
    #[arg(long, value_name = "INT")]
    pub compression_level: Option<usize>,

    /// Output block size, e.g. 1M (default: inherited from the input).
    #[arg(long, value_name = "SIZE")]
    pub block_size: Option<String>,

    /// Write a structured report.
    #[arg(long, value_name = "PATH")]
    pub report: Option<PathBuf>,

    /// Report format.
    #[arg(long, value_enum, value_name = "FORMAT")]
    pub report_format: Option<ReportFormat>,

    /// Write rejected records to this CBQ file (requires filtering).
    #[arg(long, value_name = "PATH")]
    pub failed: Option<PathBuf>,

    /// Write per-mate failure reasons to this TSV file (requires filtering).
    #[arg(long, value_name = "PATH")]
    pub failed_reasons: Option<PathBuf>,

    /// Whether the failed output holds original or processed records.
    #[arg(long, value_enum, value_name = "MODE")]
    pub failed_mode: Option<FailedMode>,

    /// Pair retention policy: strict keeps only pairs where both mates pass;
    /// orphan also writes surviving mates to single-end files.
    #[arg(long, value_enum, value_name = "POLICY")]
    pub pair_policy: Option<PairPolicy>,

    /// Path prefix for orphan outputs (writes `PREFIX.R1.cbq` and
    /// `PREFIX.R2.cbq`; requires --pair-policy orphan).
    #[arg(long, value_name = "PREFIX")]
    pub orphan_prefix: Option<PathBuf>,

    /// Overwrite existing output files.
    #[arg(long)]
    pub force: bool,

    /// Suppress the summary written to stderr.
    #[arg(short = 'q', long)]
    pub quiet: bool,
}

#[derive(Debug, Args)]
pub struct AdapterArgs {
    /// Adapter sequence trimmed from R1.
    #[arg(long, value_name = "SEQ")]
    pub adapter_r1: Option<String>,

    /// Adapter sequence trimmed from R2.
    #[arg(long, value_name = "SEQ")]
    pub adapter_r2: Option<String>,

    /// FASTA file of adapters applied to both mates.
    #[arg(long, value_name = "PATH")]
    pub adapter_fasta: Option<PathBuf>,

    /// Minimum adapter/read overlap [default: 8].
    #[arg(long, value_name = "INT")]
    pub min_overlap: Option<usize>,

    /// Maximum mismatch fraction over the overlap [default: 0.10].
    #[arg(long, value_name = "FLOAT")]
    pub max_error_rate: Option<f64>,

    /// Absolute cap on mismatches, applied in addition to --max-error-rate.
    #[arg(long, value_name = "INT")]
    pub max_errors: Option<usize>,

    /// Count insertions and deletions (in addition to substitutions) when
    /// matching adapters, using a banded edit-distance alignment.
    #[arg(long)]
    pub allow_indels: bool,

    /// Infer the insert length from the R1/R2 overlap and trim adapter
    /// read-through at the inferred boundary (paired input only).
    #[arg(long)]
    pub paired_overlap: bool,

    /// Minimum R1/R2 overlap accepted as evidence [default: 30].
    #[arg(long, value_name = "INT")]
    pub paired_overlap_min_overlap: Option<usize>,

    /// Infer adapter sequences from the data before trimming: match a
    /// built-in known-adapter library, then fall back to suffix k-mer
    /// enrichment and a consensus. Aborts below the confidence threshold.
    #[arg(long)]
    pub auto_detect: bool,

    /// Leading records sampled per mate for auto-detection [default: 100000].
    #[arg(long, value_name = "INT")]
    pub detect_sample_size: Option<usize>,

    /// Minimum confidence (support x 3' enrichment) accepted as detection
    /// evidence [default: 0.01].
    #[arg(long, value_name = "FLOAT")]
    pub detect_min_support: Option<f64>,

    /// 5' flank of a linked adapter on R1; retains the insert between the flanks.
    #[arg(long, value_name = "SEQ")]
    pub linked_5p_r1: Option<String>,
    /// 3' flank of a linked adapter on R1.
    #[arg(long, value_name = "SEQ")]
    pub linked_3p_r1: Option<String>,
    /// 5' flank of a linked adapter on R2.
    #[arg(long, value_name = "SEQ")]
    pub linked_5p_r2: Option<String>,
    /// 3' flank of a linked adapter on R2.
    #[arg(long, value_name = "SEQ")]
    pub linked_3p_r2: Option<String>,
    /// Whether both flanks must match, or either alone [default: both].
    #[arg(long, value_enum, value_name = "REQUIRE")]
    pub linked_require: Option<Require>,
    /// Bases the 5' flank may start within [default: 3].
    #[arg(long, value_name = "INT")]
    pub linked_max_5p_offset: Option<usize>,
    /// Bases of 3' flank that may hang past the read end.
    #[arg(long, value_name = "INT")]
    pub linked_max_3p_overhang: Option<usize>,
    /// Shortest insert retained by a linked match [default: 1].
    #[arg(long, value_name = "INT")]
    pub linked_min_insert_length: Option<usize>,
    /// What happens to a read no linked definition matched [default: continue].
    #[arg(long, value_enum, value_name = "POLICY")]
    pub linked_unmatched: Option<Unmatched>,
}

#[derive(Debug, Args)]
pub struct TrimArgs {
    /// Bases removed from the 5' end of both mates.
    #[arg(long, value_name = "INT", conflicts_with_all = ["front_r1", "front_r2"])]
    pub front: Option<usize>,
    /// Bases removed from the 5' end of R1.
    #[arg(long, value_name = "INT")]
    pub front_r1: Option<usize>,
    /// Bases removed from the 5' end of R2.
    #[arg(long, value_name = "INT")]
    pub front_r2: Option<usize>,

    /// Bases removed from the 3' end of both mates.
    #[arg(long, value_name = "INT", conflicts_with_all = ["tail_r1", "tail_r2"])]
    pub tail: Option<usize>,
    /// Bases removed from the 3' end of R1.
    #[arg(long, value_name = "INT")]
    pub tail_r1: Option<usize>,
    /// Bases removed from the 3' end of R2.
    #[arg(long, value_name = "INT")]
    pub tail_r2: Option<usize>,

    /// Drop 5' bases until a leading window reaches this Phred score.
    #[arg(long, value_name = "PHRED")]
    pub quality_front: Option<u8>,
    /// Window width for --quality-front [default: 4].
    #[arg(long, value_name = "INT")]
    pub quality_front_window: Option<usize>,

    /// Drop 3' bases until a trailing window reaches this Phred score.
    #[arg(long, value_name = "PHRED")]
    pub quality_tail: Option<u8>,
    /// Window width for --quality-tail [default: 4].
    #[arg(long, value_name = "INT")]
    pub quality_tail_window: Option<usize>,

    /// Truncate the read at the first window below this Phred score.
    #[arg(long, value_name = "PHRED")]
    pub quality_right: Option<u8>,
    /// Window width for --quality-right [default: 4].
    #[arg(long, value_name = "INT")]
    pub quality_right_window: Option<usize>,

    /// Remove contiguous N bases from both ends.
    #[arg(long)]
    pub trim_terminal_n: bool,

    /// Trim G-rich 3' tails.
    #[arg(long)]
    pub poly_g: bool,
    /// Shortest poly-G tail that is trimmed [default: 10].
    #[arg(long, value_name = "INT")]
    pub poly_g_min_length: Option<usize>,
    /// Mismatch fraction tolerated inside a poly-G tail [default: 0.10].
    #[arg(long, value_name = "FLOAT")]
    pub poly_g_max_mismatch_rate: Option<f64>,

    /// Trim homopolymer 3' tails of any canonical base.
    #[arg(long)]
    pub poly_x: bool,
    /// Shortest poly-X tail that is trimmed [default: 10].
    #[arg(long, value_name = "INT")]
    pub poly_x_min_length: Option<usize>,
    /// Mismatch fraction tolerated inside a poly-X tail [default: 0.10].
    #[arg(long, value_name = "FLOAT")]
    pub poly_x_max_mismatch_rate: Option<f64>,

    /// Truncate both mates to at most INT bases (a transformation, not a
    /// filter; use --length-limit to reject long reads instead).
    #[arg(long, value_name = "INT", conflicts_with_all = ["max_length_r1", "max_length_r2"])]
    pub max_length: Option<usize>,
    /// Truncate R1 to at most INT bases.
    #[arg(long, value_name = "INT")]
    pub max_length_r1: Option<usize>,
    /// Truncate R2 to at most INT bases.
    #[arg(long, value_name = "INT")]
    pub max_length_r2: Option<usize>,
}

#[derive(Debug, Args)]
pub struct FilterArgs {
    /// Reject reads shorter than INT bases.
    #[arg(long, value_name = "INT")]
    pub min_length: Option<usize>,

    /// Reject reads longer than INT bases (use --max-length to truncate
    /// instead of rejecting).
    #[arg(long, value_name = "INT")]
    pub length_limit: Option<usize>,

    /// Reject reads with more than INT ambiguous bases.
    #[arg(long, value_name = "INT")]
    pub max_n: Option<usize>,
    /// Reject reads whose ambiguous-base fraction exceeds FLOAT.
    #[arg(long, value_name = "FLOAT")]
    pub max_n_fraction: Option<f64>,

    /// Phred score at or above which a base counts as qualified [default: 15].
    #[arg(long, value_name = "PHRED")]
    pub qualified_quality: Option<u8>,
    /// Reject reads with more than INT unqualified bases.
    #[arg(long, value_name = "INT")]
    pub max_unqualified_bases: Option<usize>,
    /// Reject reads whose unqualified fraction exceeds FLOAT.
    #[arg(long, value_name = "FLOAT")]
    pub max_unqualified_fraction: Option<f64>,

    /// Reject reads whose mean Phred score is below PHRED.
    #[arg(long, value_name = "PHRED")]
    pub min_mean_quality: Option<u8>,

    /// Reject reads whose adjacent-base complexity is below FLOAT.
    #[arg(long, value_name = "FLOAT")]
    pub min_complexity: Option<f64>,
}

#[derive(Debug, Args)]
pub struct CorrectionArgs {
    /// Lowest Phred score accepted as donor evidence, inclusive [default: 30].
    #[arg(long, value_name = "PHRED")]
    pub donor_quality: Option<u8>,

    /// Highest Phred score that may be overwritten, inclusive [default: 14].
    #[arg(long, value_name = "PHRED")]
    pub recipient_quality: Option<u8>,

    /// Write a correction log to this TSV file.
    #[arg(long, value_name = "PATH")]
    pub correction_log: Option<PathBuf>,

    /// Whether the log records one row per corrected pair or per corrected base.
    #[arg(long, value_enum, value_name = "DETAIL")]
    pub correction_log_detail: Option<LogDetail>,
}

impl CorrectionArgs {
    /// Raw options, with `enabled` supplied by the caller: the `correct` command
    /// implies it, `workflow` takes it from `--correction`.
    fn options(&self, enabled: bool) -> Option<CorrectionOptions> {
        let options = CorrectionOptions {
            enabled: enabled.then_some(true),
            donor_quality: self.donor_quality,
            recipient_quality: self.recipient_quality,
            log: self.correction_log.clone(),
            log_detail: self.correction_log_detail,
        };
        (options != CorrectionOptions::default()).then_some(options)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum StepArg {
    Correct,
    Adapter,
    Trim,
    Filter,
}

impl From<StepArg> for Step {
    fn from(value: StepArg) -> Self {
        match value {
            StepArg::Correct => Self::Correct,
            StepArg::Adapter => Self::Adapter,
            StepArg::Trim => Self::Trim,
            StepArg::Filter => Self::Filter,
        }
    }
}

impl CommonArgs {
    fn output_options(&self) -> Result<OutputOptions> {
        Ok(OutputOptions {
            failed: self.failed.clone(),
            failed_reasons: self.failed_reasons.clone(),
            failed_mode: self.failed_mode,
            pair_policy: self.pair_policy,
            orphan_prefix: self.orphan_prefix.clone(),
            report: self.report.clone(),
            report_format: self.report_format,
            compression_level: self.compression_level,
            block_size: self.block_size.as_deref().map(parse_size).transpose()?,
        })
    }

    fn span(&self) -> Result<Option<Range<u64>>> {
        parse_span_opt(self.span.as_deref())
    }
}

#[derive(Debug, Args)]
pub struct SegmentArgs {
    /// Adapter sequence used as a delimiter.
    #[arg(long, value_name = "SEQ")]
    pub adapter_r1: Option<String>,

    /// FASTA file of delimiter sequences.
    #[arg(long, value_name = "PATH")]
    pub adapter_fasta: Option<PathBuf>,

    /// Minimum adapter/read overlap [default: 8].
    #[arg(long, value_name = "INT")]
    pub min_overlap: Option<usize>,

    /// Maximum mismatch fraction over the overlap [default: 0.10].
    #[arg(long, value_name = "FLOAT")]
    pub max_error_rate: Option<f64>,

    /// Absolute cap on mismatches, applied in addition to --max-error-rate.
    #[arg(long, value_name = "INT")]
    pub max_errors: Option<usize>,

    /// Count insertions and deletions when matching delimiters.
    #[arg(long)]
    pub allow_indels: bool,

    /// Keep or discard the fragments before the first and after the last
    /// delimiter; `discard` keeps only fragments flanked on both sides, which
    /// also drops reads with no delimiter at all [default: keep].
    #[arg(long, value_enum, value_name = "MODE")]
    pub terminal_fragments: Option<Terminal>,

    /// Shortest fragment emitted [default: 1].
    #[arg(long, value_name = "INT")]
    pub min_segment_length: Option<usize>,

    /// Safety limit on fragments per source read [default: 64].
    #[arg(long, value_name = "INT")]
    pub max_segments_per_read: Option<usize>,

    /// Provenance sidecar: one row per emitted fragment. Required when the input
    /// has no headers.
    #[arg(long, value_name = "PATH")]
    pub segments: Option<PathBuf>,
}

impl SegmentArgs {
    fn adapter_options(&self) -> AdapterOptions {
        AdapterOptions {
            r1: self.adapter_r1.clone(),
            fasta: self.adapter_fasta.clone(),
            min_overlap: self.min_overlap,
            max_error_rate: self.max_error_rate,
            max_errors: self.max_errors,
            allow_indels: self.allow_indels.then_some(true),
            ..AdapterOptions::default()
        }
    }

    fn options(&self) -> SegmentOptions {
        SegmentOptions {
            terminal_fragments: self.terminal_fragments,
            min_segment_length: self.min_segment_length,
            max_segments_per_read: self.max_segments_per_read,
            segments: self.segments.clone(),
        }
    }
}

impl AdapterArgs {
    fn options(&self) -> Option<AdapterOptions> {
        if self.adapter_r1.is_none()
            && self.adapter_r2.is_none()
            && self.adapter_fasta.is_none()
            && self.linked_5p_r1.is_none()
            && self.linked_3p_r1.is_none()
            && self.linked_5p_r2.is_none()
            && self.linked_3p_r2.is_none()
        {
            // Threshold-only arguments still need to reach the merge step so
            // that a TOML file can supply the adapter sequences themselves.
            if self.min_overlap.is_none()
                && self.max_error_rate.is_none()
                && self.max_errors.is_none()
                && !self.allow_indels
                && !self.paired_overlap
                && self.paired_overlap_min_overlap.is_none()
                && !self.auto_detect
                && self.detect_sample_size.is_none()
                && self.detect_min_support.is_none()
            {
                return None;
            }
        }
        // One definition per mate from the command line; a configuration file can
        // declare several.
        let linked = |five: &Option<String>, three: &Option<String>| {
            (five.is_some() || three.is_some()).then(|| {
                vec![LinkedOptions {
                    name: None,
                    five_prime: five.clone(),
                    three_prime: three.clone(),
                    require: self.linked_require,
                    max_five_prime_offset: self.linked_max_5p_offset,
                    max_three_prime_overhang: self.linked_max_3p_overhang,
                    minimum_insert_length: self.linked_min_insert_length,
                }]
            })
        };
        Some(AdapterOptions {
            linked_r1: linked(&self.linked_5p_r1, &self.linked_3p_r1),
            linked_r2: linked(&self.linked_5p_r2, &self.linked_3p_r2),
            linked_unmatched: self.linked_unmatched,
            r1: self.adapter_r1.clone(),
            r2: self.adapter_r2.clone(),
            fasta: self.adapter_fasta.clone(),
            min_overlap: self.min_overlap,
            max_error_rate: self.max_error_rate,
            max_errors: self.max_errors,
            allow_indels: self.allow_indels.then_some(true),
            paired_overlap: self.paired_overlap.then_some(true),
            paired_overlap_min_overlap: self.paired_overlap_min_overlap,
            auto_detect: self.auto_detect.then_some(true),
            detect_sample_size: self.detect_sample_size,
            detect_min_support: self.detect_min_support,
            detected: None,
        })
    }
}

impl TrimArgs {
    fn options(&self) -> Result<Option<TrimOptions>> {
        let quality_cut = |threshold: Option<u8>, window: Option<usize>| {
            if threshold.is_none() && window.is_none() {
                None
            } else {
                Some(QualityCutOptions {
                    minimum_phred: threshold,
                    window,
                })
            }
        };
        let poly = |enabled: bool, min_length: Option<usize>, rate: Option<f64>, flag: &str| {
            if !enabled && (min_length.is_some() || rate.is_some()) {
                return Err(Error::config(format!(
                    "--{flag}-min-length and --{flag}-max-mismatch-rate require --{flag}"
                )));
            }
            Ok(if enabled {
                Some(PolyOptions {
                    min_length,
                    max_mismatch_rate: rate,
                })
            } else {
                None
            })
        };

        let options = TrimOptions {
            front: self.front,
            front_r1: self.front_r1,
            front_r2: self.front_r2,
            tail: self.tail,
            tail_r1: self.tail_r1,
            tail_r2: self.tail_r2,
            quality_front: quality_cut(self.quality_front, self.quality_front_window),
            quality_right: quality_cut(self.quality_right, self.quality_right_window),
            quality_tail: quality_cut(self.quality_tail, self.quality_tail_window),
            terminal_n: if self.trim_terminal_n {
                Some(true)
            } else {
                None
            },
            poly_g: poly(
                self.poly_g,
                self.poly_g_min_length,
                self.poly_g_max_mismatch_rate,
                "poly-g",
            )?,
            poly_x: poly(
                self.poly_x,
                self.poly_x_min_length,
                self.poly_x_max_mismatch_rate,
                "poly-x",
            )?,
            max_length: self.max_length,
            max_length_r1: self.max_length_r1,
            max_length_r2: self.max_length_r2,
        };
        // Every field is `Option`, so the default value *is* "nothing was given".
        Ok((options != TrimOptions::default()).then_some(options))
    }
}

impl FilterArgs {
    fn options(&self) -> Option<FilterOptions> {
        let options = FilterOptions {
            min_length: self.min_length,
            max_length: self.length_limit,
            max_n: self.max_n,
            max_n_fraction: self.max_n_fraction,
            qualified_quality: self.qualified_quality,
            max_unqualified_bases: self.max_unqualified_bases,
            max_unqualified_fraction: self.max_unqualified_fraction,
            min_mean_quality: self.min_mean_quality,
            min_complexity: self.min_complexity,
        };
        (options != FilterOptions::default()).then_some(options)
    }
}

/// A command reduced to the library's configuration types.
struct Invocation<'a> {
    name: &'static str,
    common: &'a CommonArgs,
    config: Config,
    steps: Vec<Step>,
    require_configured: bool,
    /// Compile through `compile_segmentation` instead of the step pipeline.
    segmenting: bool,
}

impl Cli {
    fn invocation(&self) -> Result<Invocation<'_>> {
        Ok(match &self.command {
            Command::Adapter(command) => Invocation {
                name: "adapter",
                common: &command.common,
                config: Config {
                    threads: command.common.threads,
                    adapter: command.adapter.options(),
                    output: Some(command.common.output_options()?),
                    ..Config::default()
                },
                steps: vec![Step::Adapter],
                require_configured: true,
                segmenting: false,
            },
            Command::Trim(command) => Invocation {
                name: "trim",
                common: &command.common,
                config: Config {
                    threads: command.common.threads,
                    trim: command.trim.options()?,
                    output: Some(command.common.output_options()?),
                    ..Config::default()
                },
                steps: vec![Step::Trim],
                require_configured: true,
                segmenting: false,
            },
            Command::Filter(command) => Invocation {
                name: "filter",
                common: &command.common,
                config: Config {
                    threads: command.common.threads,
                    filter: command.filter.options(),
                    output: Some(command.common.output_options()?),
                    ..Config::default()
                },
                steps: vec![Step::Filter],
                require_configured: true,
                segmenting: false,
            },
            Command::Correct(command) => Invocation {
                name: "correct",
                common: &command.common,
                config: Config {
                    threads: command.common.threads,
                    correction: command.correction.options(true),
                    // Correction shares the adapter stage's overlap thresholds;
                    // the `correct` command exposes them directly.
                    adapter: (command.paired_overlap_min_overlap.is_some()
                        || command.max_error_rate.is_some())
                    .then(|| AdapterOptions {
                        paired_overlap_min_overlap: command.paired_overlap_min_overlap,
                        max_error_rate: command.max_error_rate,
                        ..AdapterOptions::default()
                    }),
                    output: Some(command.common.output_options()?),
                    ..Config::default()
                },
                steps: vec![Step::Correct],
                require_configured: true,
                segmenting: false,
            },
            Command::Segment(command) => Invocation {
                name: "segment",
                common: &command.common,
                config: Config {
                    threads: command.common.threads,
                    adapter: Some(command.segment.adapter_options()),
                    segment: Some(command.segment.options()),
                    trim: command.trim.options()?,
                    filter: command.filter.options(),
                    output: Some(command.common.output_options()?),
                    ..Config::default()
                },
                steps: vec![Step::Segment],
                require_configured: true,
                segmenting: true,
            },
            Command::Workflow(command) => return workflow_invocation(command),
            // Sniffing writes no CBQ and compiles no processing plan, so it has
            // no `Invocation`; `run` dispatches it before reaching here.
            Command::Sniff(_) => unreachable!("sniff is dispatched before invocation()"),
        })
    }
}

/// The `workflow` arm: it merges a configuration file and resolves the stage
/// selection, which is more work than the single-stage commands do.
fn workflow_invocation(command: &WorkflowCommand) -> Result<Invocation<'_>> {
    let cli = Config {
        threads: command.common.threads,
        steps: command
            .steps
            .as_ref()
            .map(|steps| steps.iter().copied().map(Into::into).collect()),
        adapter: command.adapter.options(),
        trim: command.trim.options()?,
        filter: command.filter.options(),
        correction: command.correction_args.options(command.correction),
        segment: None,
        output: Some(command.common.output_options()?),
    };
    let config = match &command.config {
        Some(path) => cli.merge(Config::from_toml_file(path)?),
        None => cli,
    };
    // `--steps` is the canonical selection; `--no-*` refines it.
    let explicit = config.steps.is_some();
    let mut steps = config.steps.clone().unwrap_or_else(|| Step::ALL.to_vec());
    for (skip, step) in [
        (command.no_adapter, Step::Adapter),
        (command.no_trim, Step::Trim),
        (command.no_filter, Step::Filter),
    ] {
        if skip {
            steps.retain(|candidate| *candidate != step);
        }
    }
    // Correction stays in the requested set even when it is off, so the stage
    // compiler can still report threshold options given without `--correction`.
    // It does not, however, count as something to do.
    let correcting = config
        .correction
        .as_ref()
        .and_then(|options| options.enabled)
        == Some(true);
    if steps.iter().all(|step| *step == Step::Correct) && !correcting {
        return Err(Error::config("every stage was disabled; nothing to do"));
    }
    Ok(Invocation {
        name: "workflow",
        common: &command.common,
        config,
        steps,
        require_configured: explicit,
        segmenting: false,
    })
}

/// Parses an optional `START..END` span argument.
fn parse_span_opt(span: Option<&str>) -> Result<Option<Range<u64>>> {
    span.map(parse_span).transpose()
}

/// Parses `START..END`, `START..` or `..END` into a record range.
fn parse_span(text: &str) -> Result<Range<u64>> {
    let invalid = || {
        Error::config(format!(
            "invalid --span '{text}': expected START..END with 0-based record indices"
        ))
    };
    let (start, end) = text.split_once("..").ok_or_else(invalid)?;
    let start = if start.is_empty() {
        0
    } else {
        start.trim().parse::<u64>().map_err(|_| invalid())?
    };
    let end = if end.is_empty() {
        u64::MAX
    } else {
        end.trim().parse::<u64>().map_err(|_| invalid())?
    };
    if start > end {
        return Err(Error::config(format!(
            "invalid --span '{text}': start ({start}) is greater than end ({end})"
        )));
    }
    Ok(start..end)
}

/// Parses a byte size with an optional `K`, `M` or `G` suffix.
fn parse_size(text: &str) -> Result<usize> {
    let text = text.trim();
    let invalid = || {
        Error::config(format!(
            "invalid size '{text}': expected an integer with optional K/M/G"
        ))
    };
    let (digits, scale) = match text.as_bytes().last() {
        Some(b'K' | b'k') => (&text[..text.len() - 1], 1024),
        Some(b'M' | b'm') => (&text[..text.len() - 1], 1024 * 1024),
        Some(b'G' | b'g') => (&text[..text.len() - 1], 1024 * 1024 * 1024),
        _ => (text, 1),
    };
    let value: usize = digits.trim().parse().map_err(|_| invalid())?;
    value.checked_mul(scale).ok_or_else(invalid)
}

/// Rejects duplicate output destinations.
///
/// Paths are compared after resolving their parent directory, so `out.cbq` and
/// `./out.cbq` are recognised as the same destination. Two outputs sharing a
/// path would make one logical output overwrite the other at commit.
fn check_distinct_outputs(paths: &[&Path]) -> Result<()> {
    let mut seen = HashSet::new();
    for path in paths {
        if !seen.insert(resolve_destination(path)) {
            return Err(Error::config(format!(
                "{} is used for more than one output",
                path.display()
            )));
        }
    }
    Ok(())
}

/// Resolves a not-yet-existing output path as far as the filesystem allows.
fn resolve_destination(path: &Path) -> PathBuf {
    let parent = match path.parent() {
        Some(parent) if parent.as_os_str().is_empty() => Path::new("."),
        Some(parent) => parent,
        None => return path.to_path_buf(),
    };
    match (std::fs::canonicalize(parent), path.file_name()) {
        (Ok(parent), Some(name)) => parent.join(name),
        _ => path.to_path_buf(),
    }
}

/// Runs a parsed command line.
pub fn run(cli: &Cli) -> Result<Outcome> {
    // Sniffing shares no machinery with the processing commands: it writes no
    // CBQ, compiles no plan and never touches the input.
    if let Command::Sniff(command) = &cli.command {
        return run_sniff(command);
    }
    run_processing(cli)?;
    Ok(Outcome::Success)
}

fn run_sniff(command: &SniffCommand) -> Result<Outcome> {
    match &command.kind {
        SniffKind::Adapters(command) => run_sniff_adapters(command),
        #[cfg(feature = "sniff-strand")]
        SniffKind::Strand(command) => run_sniff_strand(command),
    }
}

#[cfg(feature = "sniff-strand")]
fn run_sniff_strand(command: &SniffStrandCommand) -> Result<Outcome> {
    use crate::sniff::strand;

    let common = &command.common;
    let input = CbqInput::open(&common.input)?;
    let span = common.span()?;
    let defaults = strand::Params::default();
    let params = strand::Params {
        sample_size: common.sample_size.unwrap_or(defaults.sample_size),
        target_informative: command
            .target_informative
            .unwrap_or(defaults.target_informative),
        min_informative: command.min_informative.unwrap_or(defaults.min_informative),
        min_informative_fraction: command
            .min_informative_fraction
            .unwrap_or(defaults.min_informative_fraction),
        stranded_threshold: command
            .stranded_threshold
            .unwrap_or(defaults.stranded_threshold),
        unstranded_threshold: command
            .unstranded_threshold
            .unwrap_or(defaults.unstranded_threshold),
    }
    .validate()?;

    let threads = resolve_threads(common.threads);
    let index_path = match (&command.index, &command.transcriptome) {
        (Some(index), None) => index.clone(),
        (None, Some(fasta)) => strand::build_temp_index(fasta, threads)?,
        _ => unreachable!("clap ArgGroup enforces exactly one of --index/--transcriptome"),
    };
    let (result, plan) = strand::sniff(
        &input,
        &index_path,
        params,
        span.as_ref(),
        threads,
    )?;
    let confident = result.is_confident();
    let report = crate::sniff::report::StrandReport::new(&input, &plan, params, result);
    emit(&report.render(common.format)?, common)?;

    Ok(if common.require_confident && !confident {
        Outcome::NotConfident
    } else {
        Outcome::Success
    })
}

fn run_sniff_adapters(command: &SniffAdaptersCommand) -> Result<Outcome> {
    let common = &command.common;
    if let (Some(report), Some(config)) = (common.output.as_deref(), command.emit_config.as_deref())
    {
        check_distinct_outputs(&[report, config])?;
    }
    let input = CbqInput::open(&common.input)?;
    let span = common.span()?;
    let defaults = crate::sniff::adapters::Params::default();
    let params = crate::sniff::adapters::Params {
        sample_size: common.sample_size.unwrap_or(defaults.sample_size),
        top: command.top.unwrap_or(defaults.top),
        ..defaults
    }
    .validate()?;

    let (result, plan) = crate::sniff::adapters::sniff(
        &input,
        params,
        span.as_ref(),
        resolve_threads(common.threads),
    )?;
    let confident = result.decision().is_confident();
    let report = crate::sniff::report::AdapterReport::new(&input, &plan, params, result);
    let rendered_report = report.render(common.format)?;
    if let Some(path) = command.emit_config.as_deref() {
        emit_to(path, &render_adapter_config(&report.result)?, common)?;
    }
    emit(&rendered_report, common)?;

    Ok(if common.require_confident && !confident {
        Outcome::NotConfident
    } else {
        Outcome::Success
    })
}

/// Writes the detected adapters as a `bqc` configuration fragment.
///
/// Only a uniquely confident result is written. A mixed or inconclusive one is
/// an error rather than an empty file: the caller asked for a configuration to
/// feed the next command, and silently producing one that trims nothing — or
/// one that picks arbitrarily between two libraries — would be worse than
/// stopping. Only the values needed to reproduce the recommendation go in;
/// detector internals stay in the report.
fn render_adapter_config(result: &crate::sniff::adapters::AdapterSniff) -> Result<String> {
    use std::fmt::Write as _;
    let Some((r1, r2)) = result.recommendation() else {
        return Err(Error::config(format!(
            "--emit-config needs a uniquely confident result, but the outcome was {}; \
             inspect the candidates in the report and choose an adapter explicitly",
            result.decision().name()
        )));
    };
    let mut fragment = String::from("# Written by `bqc sniff adapters`.\n[adapter]\n");
    let _ = writeln!(fragment, "r1 = \"{r1}\"");
    if let Some(r2) = r2 {
        let _ = writeln!(fragment, "r2 = \"{r2}\"");
    }
    Ok(fragment)
}

/// Writes a rendered report to a file.
///
/// The file is written atomically, through the same temporary-then-rename path
/// every other output uses, so a failed run leaves nothing behind.
fn emit_to(path: &Path, rendered: &str, common: &SniffArgs) -> Result<()> {
    let mut file = TextOutput::create(path, &common.input, common.force)?;
    file.write(rendered.as_bytes())?;
    file.finish()
}

/// Writes a rendered report to its destination.
fn emit(rendered: &str, common: &SniffArgs) -> Result<()> {
    use std::io::Write as _;
    if let Some(path) = common.output.as_deref() {
        return emit_to(path, rendered, common);
    }
    let mut stdout = std::io::stdout().lock();
    stdout
        .write_all(rendered.as_bytes())
        .and_then(|()| stdout.flush())
        .map_err(|e| Error::write(Path::new("<stdout>"), e))
}

impl SniffArgs {
    /// Parses `--span`.
    fn span(&self) -> Result<Option<Range<u64>>> {
        parse_span_opt(self.span.as_deref())
    }
}

fn run_processing(cli: &Cli) -> Result<()> {
    let started = std::time::Instant::now();
    let invocation = cli.invocation()?;
    let common = invocation.common;

    let input = CbqInput::open(&common.input)?;
    let span = common.span()?;

    // Auto-detection runs before compilation so the inferred adapters become
    // part of the resolved plan and report. It is the same detector
    // `sniff adapters` reports, invoked here for its recommendation instead.
    let mut config = invocation.config;
    let detection = match config.adapter.as_mut() {
        Some(adapter) if adapter.auto_detect == Some(true) => {
            let (detection, _) = crate::sniff::adapters::sniff(
                &input,
                crate::config::detection_params(adapter)?,
                span.as_ref(),
                resolve_threads(common.threads),
            )?;
            adapter.detected = Some(detection);
            adapter.detected.clone()
        }
        _ => None,
    };

    let plan = if invocation.segmenting {
        config.compile_segmentation(&input.header())?
    } else {
        config.compile(
            &invocation.steps,
            invocation.require_configured,
            &input.header(),
        )?
    };

    let mut paths = vec![common.output.as_path()];
    paths.extend(plan.output.failed.as_deref());
    paths.extend(plan.output.failed_reasons.as_deref());
    paths.extend(plan.output.orphan_r1.as_deref());
    paths.extend(plan.output.orphan_r2.as_deref());
    paths.extend(plan.output.correction_log.as_deref());
    paths.extend(plan.output.segments.as_deref());
    paths.extend(plan.output.report.as_deref());
    check_distinct_outputs(&paths)?;

    let mut destinations = Destinations::create(common, &plan, &input)?;
    let stats = {
        let mut outputs = Outputs {
            accepted: &mut destinations.accepted,
            failed: destinations.failed.as_mut(),
            reasons: destinations.reasons.as_mut(),
            orphan_r1: destinations.orphan_r1.as_mut(),
            orphan_r2: destinations.orphan_r2.as_mut(),
            corrections: destinations.corrections.as_mut(),
            segments: destinations.segments.as_mut(),
        };
        let options = RunOptions {
            threads: plan.threads,
            span: span.clone(),
            failed_mode: plan.output.failed_mode,
        };
        crate::engine::run(&input, &plan.workflow, &mut outputs, &options)?
    };
    destinations.finish()?;

    let report = Report::new(RunReport {
        command: invocation.name,
        input: &input,
        plan: &plan,
        span,
        accepted: &common.output,
        stats: &stats,
        elapsed: started.elapsed(),
        detection: detection.as_ref(),
    });
    write_report(&report, &plan, common, destinations.report.as_mut())
}

/// Every file a run writes to, created before processing starts.
struct Destinations {
    accepted: CbqOutput,
    failed: Option<CbqOutput>,
    orphan_r1: Option<CbqOutput>,
    orphan_r2: Option<CbqOutput>,
    reasons: Option<TextOutput>,
    corrections: Option<TextOutput>,
    segments: Option<TextOutput>,
    report: Option<TextOutput>,
}

impl Destinations {
    fn create(common: &CommonArgs, plan: &Plan, input: &CbqInput) -> Result<Self> {
        let schema = input.schema();
        let header = schema.to_file_header(plan.output.compression_level, plan.output.block_size);
        // Orphan outputs are single-end; the rest of the schema is preserved.
        let orphan_header = schema
            .unpaired()
            .to_file_header(plan.output.compression_level, plan.output.block_size);
        let cbq = |path: Option<&Path>, header| {
            path.map(|path| CbqOutput::create(path, &common.input, header, common.force))
                .transpose()
        };
        let text = |path: Option<&Path>| {
            path.map(|path| TextOutput::create(path, &common.input, common.force))
                .transpose()
        };
        Ok(Self {
            accepted: CbqOutput::create(&common.output, &common.input, header, common.force)?,
            failed: cbq(plan.output.failed.as_deref(), header)?,
            orphan_r1: cbq(plan.output.orphan_r1.as_deref(), orphan_header)?,
            orphan_r2: cbq(plan.output.orphan_r2.as_deref(), orphan_header)?,
            reasons: text(plan.output.failed_reasons.as_deref())?,
            corrections: text(plan.output.correction_log.as_deref())?,
            segments: text(plan.output.segments.as_deref())?,
            // The report is written last, but its destination is reserved up
            // front: a report path conflict must abort the run before any
            // processing happens, not after every output has been committed.
            report: text(plan.output.report.as_deref())?,
        })
    }

    /// Commits every output except the report, which is written afterwards.
    fn finish(&mut self) -> Result<()> {
        self.accepted.finish()?;
        for output in [
            self.failed.as_mut(),
            self.orphan_r1.as_mut(),
            self.orphan_r2.as_mut(),
        ]
        .into_iter()
        .flatten()
        {
            output.finish()?;
        }
        for output in [
            self.reasons.as_mut(),
            self.corrections.as_mut(),
            self.segments.as_mut(),
        ]
        .into_iter()
        .flatten()
        {
            output.finish()?;
        }
        Ok(())
    }
}

fn write_report(
    report: &Report,
    plan: &Plan,
    common: &CommonArgs,
    file: Option<&mut TextOutput>,
) -> Result<()> {
    if let Some(file) = file {
        file.write(report.render(plan.output.report_format)?.as_bytes())?;
        file.finish()?;
    }
    if !common.quiet {
        let mut stderr = std::io::stderr().lock();
        report
            .write_summary(&mut stderr)
            .map_err(|e| Error::write(Path::new("<stderr>"), e))?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;

    #[test]
    fn cli_definition_is_valid() {
        Cli::command().debug_assert();
    }

    fn parse(args: &[&str]) -> Cli {
        Cli::try_parse_from(args).unwrap()
    }

    #[test]
    fn span_parsing_accepts_open_ends() {
        assert_eq!(parse_span("0..100").unwrap(), 0..100);
        assert_eq!(parse_span("10..").unwrap(), 10..u64::MAX);
        assert_eq!(parse_span("..50").unwrap(), 0..50);
        assert_eq!(parse_span("..").unwrap(), 0..u64::MAX);
        assert!(parse_span("100..10").is_err());
        assert!(parse_span("abc").is_err());
        assert!(parse_span("1-2").is_err());
    }

    #[test]
    fn size_parsing_accepts_suffixes() {
        assert_eq!(parse_size("1024").unwrap(), 1024);
        assert_eq!(parse_size("1K").unwrap(), 1024);
        assert_eq!(parse_size("2m").unwrap(), 2 * 1024 * 1024);
        assert_eq!(parse_size("1G").unwrap(), 1024 * 1024 * 1024);
        assert!(parse_size("").is_err());
        assert!(parse_size("1X").is_err());
    }

    #[test]
    fn front_alias_stays_distinct_from_the_per_mate_options() {
        // The alias must not masquerade as an R2-specific option, which would
        // make `--front` unusable on single-end files.
        let cli = parse(&["bqc", "trim", "in.cbq", "-o", "out.cbq", "--front", "5"]);
        let Command::Trim(command) = &cli.command else {
            panic!("expected trim")
        };
        let options = command.trim.options().unwrap().unwrap();
        assert_eq!(options.front, Some(5));
        assert_eq!(options.front_r1, None);
        assert_eq!(options.front_r2, None);
    }

    #[test]
    fn mate_specific_arguments_conflict_with_the_alias() {
        assert!(Cli::try_parse_from([
            "bqc",
            "trim",
            "in.cbq",
            "-o",
            "out.cbq",
            "--front",
            "5",
            "--front-r1",
            "3",
        ])
        .is_err());
    }

    #[test]
    fn poly_parameters_require_their_flag() {
        let cli = parse(&[
            "bqc",
            "trim",
            "in.cbq",
            "-o",
            "out.cbq",
            "--poly-g-min-length",
            "12",
        ]);
        let Command::Trim(command) = &cli.command else {
            panic!("expected trim")
        };
        let err = command.trim.options().unwrap_err();
        assert!(format!("{err}").contains("require --poly-g"), "{err}");
    }

    #[test]
    fn filter_length_limit_maps_to_the_max_length_predicate() {
        let cli = parse(&[
            "bqc",
            "filter",
            "in.cbq",
            "-o",
            "out.cbq",
            "--min-length",
            "30",
            "--length-limit",
            "150",
        ]);
        let Command::Filter(command) = &cli.command else {
            panic!("expected filter")
        };
        let options = command.filter.options().unwrap();
        assert_eq!(options.min_length, Some(30));
        assert_eq!(options.max_length, Some(150));
    }

    #[test]
    fn empty_filter_and_correction_arguments_produce_no_options() {
        let cli = parse(&["bqc", "workflow", "in.cbq", "-o", "out.cbq"]);
        let Command::Workflow(command) = &cli.command else {
            panic!("expected workflow")
        };
        assert!(command.filter.options().is_none());
        assert!(command.correction_args.options(false).is_none());
    }

    #[test]
    fn one_correction_argument_produces_options() {
        let cli = parse(&[
            "bqc",
            "workflow",
            "in.cbq",
            "-o",
            "out.cbq",
            "--donor-quality",
            "35",
        ]);
        let Command::Workflow(command) = &cli.command else {
            panic!("expected workflow")
        };
        let options = command.correction_args.options(false).unwrap();
        assert_eq!(options.donor_quality, Some(35));
        assert_eq!(options.enabled, None);
    }

    #[test]
    fn workflow_steps_and_negations_resolve() {
        let cli = parse(&[
            "bqc",
            "workflow",
            "in.cbq",
            "-o",
            "out.cbq",
            "--adapter-r1",
            "ACGTACGTACGT",
            "--min-length",
            "30",
        ]);
        let invocation = cli.invocation().unwrap();
        assert_eq!(
            invocation.steps,
            Step::ALL.to_vec(),
            "every stage is requested; whether each runs is decided at compile time"
        );
        assert!(!invocation.require_configured);

        let cli = parse(&[
            "bqc",
            "workflow",
            "in.cbq",
            "-o",
            "out.cbq",
            "--steps",
            "adapter,filter",
            "--adapter-r1",
            "ACGTACGTACGT",
            "--min-length",
            "30",
        ]);
        let invocation = cli.invocation().unwrap();
        assert_eq!(invocation.steps, vec![Step::Adapter, Step::Filter]);
        assert!(
            invocation.require_configured,
            "explicit steps must be configured"
        );

        let cli = parse(&[
            "bqc",
            "workflow",
            "in.cbq",
            "-o",
            "out.cbq",
            "--no-adapter",
            "--min-length",
            "30",
        ]);
        let invocation = cli.invocation().unwrap();
        assert_eq!(
            invocation.steps,
            vec![Step::Correct, Step::Trim, Step::Filter],
            "--no-adapter removes only the adapter stage"
        );
    }

    #[test]
    fn workflow_rejects_disabling_every_stage() {
        let cli = parse(&[
            "bqc",
            "workflow",
            "in.cbq",
            "-o",
            "out.cbq",
            "--no-adapter",
            "--no-trim",
            "--no-filter",
        ]);
        assert!(cli.invocation().is_err());
    }

    #[test]
    fn steps_and_negations_cannot_be_combined() {
        assert!(Cli::try_parse_from([
            "bqc",
            "workflow",
            "in.cbq",
            "-o",
            "out.cbq",
            "--steps",
            "adapter",
            "--no-trim",
        ])
        .is_err());
    }

    #[test]
    fn duplicate_output_destinations_are_rejected() {
        let dir = tempfile::tempdir().unwrap();
        let a = dir.path().join("a.cbq");
        let b = dir.path().join("b.cbq");
        assert!(check_distinct_outputs(&[a.as_path(), b.as_path()]).is_ok());
        assert!(check_distinct_outputs(&[a.as_path(), a.as_path()]).is_err());

        // The same destination spelled two ways is still one destination.
        let spelled_differently = dir.path().join(".").join("a.cbq");
        assert!(
            check_distinct_outputs(&[a.as_path(), spelled_differently.as_path()]).is_err(),
            "two spellings of one path would target the same file"
        );
    }

    #[test]
    fn adapter_thresholds_alone_still_reach_the_merge_step() {
        let cli = parse(&[
            "bqc",
            "workflow",
            "in.cbq",
            "-o",
            "out.cbq",
            "--min-overlap",
            "10",
        ]);
        let Command::Workflow(command) = &cli.command else {
            panic!("expected workflow")
        };
        let options = command.adapter.options().unwrap();
        assert_eq!(options.min_overlap, Some(10));
        assert!(options.r1.is_none());
    }
}
