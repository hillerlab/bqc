// Copyright (c) 2026 The Hiller Lab at the Senckenberg Gesellschaft für Naturforschung
// Distributed under the terms of the GNU General Public License, Version 3.0.

//! Raw configuration and its compilation into a runnable [`Plan`].
//!
//! There are exactly two configuration layers:
//!
//! 1. [`Config`] — the raw, fully optional form. TOML files deserialize into
//!    it and the CLI converts its arguments into it, so both front ends share
//!    one validation path. CLI values override file values field by field.
//! 2. [`Plan`] — the compiled form: every value resolved, every adapter
//!    normalized, every threshold checked. The record loop only ever sees this.
//!
//! The compiled plan is what the JSON report serializes, so a report always
//! documents the complete resolved configuration rather than the flags the user
//! happened to type.

use std::path::{Path, PathBuf};

use binseq::cbq::FileHeader;
use serde::{Deserialize, Serialize};

use crate::adapter::{Adapter, AdapterParams, AdapterStage, read_adapter_fasta};
use crate::correct::{CorrectionStage, LogDetail};
use crate::error::{Error, Result};
use crate::filter::FilterStage;
use crate::io::Schema;
use crate::linked::{DEFAULT_FIVE_PRIME_OFFSET, LinkedAdapter, LinkedStage, Require, Unmatched};
use crate::process::{FailedMode, PairPolicy, Workflow};
use crate::segment::Terminal;
use crate::sniff::adapters::{AdapterSniff, Params as SniffParams};
use crate::trim::{MateTrim, PolyParams, QualityCut, TrimStage};
use crate::umi::{UmiLocation, UmiStage};

/// Default sliding-window width for quality trimming.
pub const DEFAULT_QUALITY_WINDOW: usize = 4;
/// Default Phred threshold separating qualified from unqualified bases.
pub const DEFAULT_QUALIFIED_QUALITY: u8 = 15;

/// A processing stage.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Step {
    Correct,
    Adapter,
    Trim,
    Filter,
    Umi,
    /// Internal adapter splitting. Reachable only through the `segment` command,
    /// so it is deliberately absent from [`Step::ALL`] and from `--steps`: it
    /// changes output cardinality and cannot be composed with the rest.
    Segment,
}

impl Step {
    pub const ALL: [Step; 5] = [
        Step::Correct,
        Step::Adapter,
        Step::Trim,
        Step::Filter,
        Step::Umi,
    ];
}

/// Structured report format.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize, clap::ValueEnum)]
#[serde(rename_all = "lowercase")]
pub enum ReportFormat {
    #[default]
    Json,
    Tsv,
}

/// Merges a lower-priority configuration into a higher-priority one.
pub(crate) trait Merge {
    /// Keeps every value set on `self`, filling gaps from `fallback`.
    #[must_use]
    fn merge(self, fallback: Self) -> Self;
}

fn merge_option<T: Merge>(primary: Option<T>, fallback: Option<T>) -> Option<T> {
    match (primary, fallback) {
        (Some(primary), Some(fallback)) => Some(primary.merge(fallback)),
        (primary, fallback) => primary.or(fallback),
    }
}

/// Implements field-wise precedence without serializing paths or runtime state.
///
/// `direct` fields use the primary value when present, `nested` fields recurse,
/// and `primary` fields never inherit from the fallback. The generated struct
/// literal is deliberately exhaustive so adding a field requires deciding its
/// merge behaviour at compile time.
macro_rules! impl_merge {
    (
        $type:ident {
            direct: [$($direct:ident),* $(,)?],
            nested: [$($nested:ident),* $(,)?],
            primary: [$($primary:ident),* $(,)?],
        }
    ) => {
        impl Merge for $type {
            fn merge(self, fallback: Self) -> Self {
                Self {
                    $($direct: self.$direct.or(fallback.$direct),)*
                    $($nested: merge_option(self.$nested, fallback.$nested),)*
                    $($primary: self.$primary,)*
                }
            }
        }
    };
}

/// The raw configuration shared by the CLI and TOML front ends.
#[derive(Debug, Clone, Default, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Config {
    pub threads: Option<usize>,
    pub steps: Option<Vec<Step>>,
    pub adapter: Option<AdapterOptions>,
    pub trim: Option<TrimOptions>,
    pub filter: Option<FilterOptions>,
    pub correction: Option<CorrectionOptions>,
    pub umi: Option<UmiOptions>,
    pub segment: Option<SegmentOptions>,
    pub output: Option<OutputOptions>,
}

impl_merge!(Config {
    direct: [threads, steps],
    nested: [adapter, trim, filter, correction, umi, segment, output],
    primary: [],
});

/// Internal adapter splitting options.
///
/// The delimiter sequences and matcher thresholds come from `[adapter]`: a
/// delimiter is matched exactly the way a 3' adapter is, so there is one set of
/// adapter options and one matcher, not two.
#[derive(Debug, Clone, Default, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SegmentOptions {
    pub terminal_fragments: Option<Terminal>,
    pub min_segment_length: Option<usize>,
    pub max_segments_per_read: Option<usize>,
    /// Provenance sidecar path.
    pub segments: Option<PathBuf>,
}

impl_merge!(SegmentOptions {
    direct: [
        terminal_fragments,
        min_segment_length,
        max_segments_per_read,
        segments,
    ],
    nested: [],
    primary: [],
});

/// Adapter removal options.
#[derive(Debug, Clone, Default, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AdapterOptions {
    pub r1: Option<String>,
    pub r2: Option<String>,
    pub fasta: Option<PathBuf>,
    pub min_overlap: Option<usize>,
    pub max_error_rate: Option<f64>,
    pub max_errors: Option<usize>,
    pub allow_indels: Option<bool>,
    pub paired_overlap: Option<bool>,
    pub paired_overlap_min_overlap: Option<usize>,
    pub auto_detect: Option<bool>,
    pub detect_sample_size: Option<usize>,
    pub detect_min_support: Option<f64>,
    /// Linked-adapter definitions for R1.
    pub linked_r1: Option<Vec<LinkedOptions>>,
    /// Linked-adapter definitions for R2.
    pub linked_r2: Option<Vec<LinkedOptions>>,
    /// What happens to a read no linked definition matched.
    pub linked_unmatched: Option<Unmatched>,
    /// Detection outcome, filled in by the front end before compilation.
    #[serde(skip)]
    pub detected: Option<AdapterSniff>,
}

impl AdapterOptions {
    fn is_empty(&self) -> bool {
        self.r1.is_none()
            && self.r2.is_none()
            && self.fasta.is_none()
            && self.paired_overlap != Some(true)
            && self.auto_detect != Some(true)
    }

    /// Whether any linked definition was supplied.
    fn has_linked(&self) -> bool {
        self.linked_r1.as_ref().is_some_and(|list| !list.is_empty())
            || self.linked_r2.as_ref().is_some_and(|list| !list.is_empty())
    }
}

impl_merge!(AdapterOptions {
    direct: [
        r1,
        r2,
        fasta,
        min_overlap,
        max_error_rate,
        max_errors,
        allow_indels,
        paired_overlap,
        paired_overlap_min_overlap,
        auto_detect,
        detect_sample_size,
        detect_min_support,
        linked_r1,
        linked_r2,
        linked_unmatched,
    ],
    nested: [],
    // Detection is runtime state, not fallback configuration.
    primary: [detected],
});

/// A windowed quality cut. Presence of the table enables the operation.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct QualityCutOptions {
    pub minimum_phred: Option<u8>,
    pub window: Option<usize>,
}

impl_merge!(QualityCutOptions {
    direct: [minimum_phred, window],
    nested: [],
    primary: [],
});

impl QualityCutOptions {
    fn compile(self, flag: &str) -> Result<QualityCut> {
        let Some(minimum_phred) = self.minimum_phred else {
            return Err(Error::config(format!(
                "--{flag}-window requires --{flag} <PHRED>"
            )));
        };
        QualityCut {
            minimum_phred,
            window: self.window.unwrap_or(DEFAULT_QUALITY_WINDOW),
        }
        .validate(flag)
    }
}

/// Homopolymer tail trimming. Presence of the table enables the operation.
#[derive(Debug, Clone, Copy, Default, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PolyOptions {
    pub min_length: Option<usize>,
    pub max_mismatch_rate: Option<f64>,
}

impl_merge!(PolyOptions {
    direct: [min_length, max_mismatch_rate],
    nested: [],
    primary: [],
});

impl PolyOptions {
    fn compile(self, flag: &str) -> Result<PolyParams> {
        let defaults = PolyParams::default();
        PolyParams {
            min_length: self.min_length.unwrap_or(defaults.min_length),
            max_mismatch_rate: self.max_mismatch_rate.unwrap_or(defaults.max_mismatch_rate),
        }
        .validate(flag)
    }
}

/// Trimming options.
#[derive(Debug, Clone, Copy, Default, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TrimOptions {
    /// Applies to both mates; the per-mate keys below take precedence.
    pub front: Option<usize>,
    pub front_r1: Option<usize>,
    pub front_r2: Option<usize>,
    /// Applies to both mates; the per-mate keys below take precedence.
    pub tail: Option<usize>,
    pub tail_r1: Option<usize>,
    pub tail_r2: Option<usize>,
    pub quality_front: Option<QualityCutOptions>,
    pub quality_right: Option<QualityCutOptions>,
    pub quality_tail: Option<QualityCutOptions>,
    pub terminal_n: Option<bool>,
    pub poly_g: Option<PolyOptions>,
    pub poly_x: Option<PolyOptions>,
    /// Applies to both mates; the per-mate keys below take precedence.
    pub max_length: Option<usize>,
    pub max_length_r1: Option<usize>,
    pub max_length_r2: Option<usize>,
}

impl_merge!(TrimOptions {
    direct: [
        front,
        front_r1,
        front_r2,
        tail,
        tail_r1,
        tail_r2,
        terminal_n,
        max_length,
        max_length_r1,
        max_length_r2,
    ],
    nested: [quality_front, quality_right, quality_tail, poly_g, poly_x],
    primary: [],
});

impl TrimOptions {
    /// Whether a mate-2 specific value is set, which requires a paired input.
    ///
    /// The both-mates aliases (`front`, `tail`, `max_length`) are deliberately
    /// excluded: they are meaningful on single-end files.
    fn uses_r2_only_options(&self) -> bool {
        self.front_r2.is_some() || self.tail_r2.is_some() || self.max_length_r2.is_some()
    }

    fn compile(self) -> Result<TrimStage> {
        if self.quality_right.is_some() && self.quality_tail.is_some() {
            return Err(Error::config(
                "--quality-right and --quality-tail are mutually exclusive: \
                 --quality-right already truncates the 3' end",
            ));
        }
        let shared = MateTrim {
            quality_front: self
                .quality_front
                .map(|cut| cut.compile("quality-front"))
                .transpose()?,
            quality_right: self
                .quality_right
                .map(|cut| cut.compile("quality-right"))
                .transpose()?,
            quality_tail: self
                .quality_tail
                .map(|cut| cut.compile("quality-tail"))
                .transpose()?,
            terminal_n: self.terminal_n.unwrap_or(false),
            poly_g: self.poly_g.map(|poly| poly.compile("poly-g")).transpose()?,
            poly_x: self.poly_x.map(|poly| poly.compile("poly-x")).transpose()?,
            ..MateTrim::default()
        };
        let stage = TrimStage {
            r1: MateTrim {
                front: self.front_r1.or(self.front).unwrap_or(0),
                tail: self.tail_r1.or(self.tail).unwrap_or(0),
                max_length: self.max_length_r1.or(self.max_length),
                ..shared
            },
            r2: MateTrim {
                front: self.front_r2.or(self.front).unwrap_or(0),
                tail: self.tail_r2.or(self.tail).unwrap_or(0),
                max_length: self.max_length_r2.or(self.max_length),
                ..shared
            },
        };
        if stage.r1.max_length == Some(0) || stage.r2.max_length == Some(0) {
            return Err(Error::config("--max-length must be at least 1"));
        }
        if stage.r1.is_noop() && stage.r2.is_noop() {
            return Err(Error::config(
                "trimming requires at least one operation \
                 (--front, --tail, --quality-front, --quality-right, --quality-tail, \
                 --trim-terminal-n, --poly-g, --poly-x or --max-length)",
            ));
        }
        Ok(stage)
    }
}

/// Filtering options.
#[derive(Debug, Clone, Copy, Default, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FilterOptions {
    pub min_length: Option<usize>,
    pub max_length: Option<usize>,
    pub max_n: Option<usize>,
    pub max_n_fraction: Option<f64>,
    pub qualified_quality: Option<u8>,
    pub max_unqualified_bases: Option<usize>,
    pub max_unqualified_fraction: Option<f64>,
    pub min_mean_quality: Option<u8>,
    pub min_complexity: Option<f64>,
}

impl_merge!(FilterOptions {
    direct: [
        min_length,
        max_length,
        max_n,
        max_n_fraction,
        qualified_quality,
        max_unqualified_bases,
        max_unqualified_fraction,
        min_mean_quality,
        min_complexity,
    ],
    nested: [],
    primary: [],
});

impl FilterOptions {
    fn compile(self) -> Result<FilterStage> {
        if self.qualified_quality.is_some()
            && self.max_unqualified_bases.is_none()
            && self.max_unqualified_fraction.is_none()
        {
            return Err(Error::config(
                "--qualified-quality requires --max-unqualified-bases \
                 or --max-unqualified-fraction",
            ));
        }
        FilterStage {
            min_length: self.min_length,
            max_length: self.max_length,
            max_n: self.max_n,
            max_n_fraction: self.max_n_fraction,
            qualified_quality: self.qualified_quality.unwrap_or(DEFAULT_QUALIFIED_QUALITY),
            max_unqualified_bases: self.max_unqualified_bases,
            max_unqualified_fraction: self.max_unqualified_fraction,
            min_mean_quality: self.min_mean_quality,
            min_complexity: self.min_complexity,
        }
        .validate()
    }
}

/// One named linked-adapter definition, as it appears in TOML.
#[derive(Debug, Clone, Default, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct LinkedOptions {
    pub name: Option<String>,
    pub five_prime: Option<String>,
    pub three_prime: Option<String>,
    pub require: Option<Require>,
    pub max_five_prime_offset: Option<usize>,
    pub max_three_prime_overhang: Option<usize>,
    pub minimum_insert_length: Option<usize>,
}

impl LinkedOptions {
    /// Compiles one definition, naming it after its position when unnamed.
    fn compile(&self, position: usize) -> Result<LinkedAdapter> {
        let name = self
            .name
            .clone()
            .unwrap_or_else(|| format!("linked{}", position + 1));
        let (Some(five), Some(three)) = (self.five_prime.as_ref(), self.three_prime.as_ref())
        else {
            return Err(Error::config(format!(
                "linked adapter '{name}' needs both a 5' and a 3' sequence"
            )));
        };
        Ok(LinkedAdapter {
            five_prime: Adapter::new(format!("{name}.5p"), five.as_bytes())?,
            three_prime: Adapter::new(format!("{name}.3p"), three.as_bytes())?,
            name,
            require: self.require.unwrap_or(Require::Both),
            max_five_prime_offset: self
                .max_five_prime_offset
                .unwrap_or(DEFAULT_FIVE_PRIME_OFFSET),
            max_three_prime_overhang: self.max_three_prime_overhang,
            minimum_insert_length: self.minimum_insert_length.unwrap_or(1),
        })
    }
}

/// Paired-overlap base correction options.
#[derive(Debug, Clone, Default, PartialEq, Eq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CorrectionOptions {
    pub enabled: Option<bool>,
    pub donor_quality: Option<u8>,
    pub recipient_quality: Option<u8>,
    pub log: Option<PathBuf>,
    pub log_detail: Option<LogDetail>,
}

impl_merge!(CorrectionOptions {
    direct: [enabled, donor_quality, recipient_quality, log, log_detail,],
    nested: [],
    primary: [],
});

/// UMI extraction options.
#[derive(Debug, Clone, Default, PartialEq, Eq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct UmiOptions {
    pub location: Option<UmiLocation>,
    pub length: Option<usize>,
    pub skip: Option<usize>,
    pub prefix: Option<String>,
    pub delimiter: Option<String>,
}

impl_merge!(UmiOptions {
    direct: [location, length, skip, prefix, delimiter],
    nested: [],
    primary: [],
});

/// Output options.
#[derive(Debug, Clone, Default, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct OutputOptions {
    pub failed: Option<PathBuf>,
    pub failed_reasons: Option<PathBuf>,
    pub failed_mode: Option<FailedMode>,
    pub pair_policy: Option<PairPolicy>,
    pub orphan_prefix: Option<PathBuf>,
    pub report: Option<PathBuf>,
    pub report_format: Option<ReportFormat>,
    pub compression_level: Option<usize>,
    pub block_size: Option<usize>,
}

impl_merge!(OutputOptions {
    direct: [
        failed,
        failed_reasons,
        failed_mode,
        pair_policy,
        orphan_prefix,
        report,
        report_format,
        compression_level,
        block_size,
    ],
    nested: [],
    primary: [],
});

/// Resolves a thread count: `None` or `0` means every available core.
pub(crate) fn resolve_threads(threads: Option<usize>) -> usize {
    match threads {
        None | Some(0) => std::thread::available_parallelism().map_or(1, Into::into),
        Some(threads) => threads,
    }
}

/// Derives the orphan output paths for a prefix: `<prefix>.R1.cbq` and
/// `<prefix>.R2.cbq`.
#[must_use]
fn orphan_paths(prefix: &Path) -> [PathBuf; 2] {
    let base = prefix.as_os_str();
    [
        PathBuf::from(format!("{}.R1.cbq", base.to_string_lossy())),
        PathBuf::from(format!("{}.R2.cbq", base.to_string_lossy())),
    ]
}

/// Resolved output configuration.
#[derive(Debug, Clone, Serialize)]
pub struct ResolvedOutput {
    pub failed: Option<PathBuf>,
    pub failed_reasons: Option<PathBuf>,
    /// Correction log path, when one was requested.
    pub correction_log: Option<PathBuf>,
    /// Segmentation provenance sidecar path, when one was requested.
    pub segments: Option<PathBuf>,
    pub failed_mode: FailedMode,
    pub pair_policy: PairPolicy,
    pub orphan_r1: Option<PathBuf>,
    pub orphan_r2: Option<PathBuf>,
    pub report: Option<PathBuf>,
    pub report_format: ReportFormat,
    pub compression_level: usize,
    pub block_size: usize,
}

/// A fully compiled run configuration.
#[derive(Debug, Clone, Serialize)]
pub struct Plan {
    pub steps: Vec<Step>,
    pub threads: usize,
    pub workflow: Workflow,
    pub output: ResolvedOutput,
}

impl Config {
    /// Loads a TOML configuration file.
    pub fn from_toml_file(path: &Path) -> Result<Self> {
        let text = std::fs::read_to_string(path).map_err(|e| Error::read(path, e))?;
        toml::from_str(&text).map_err(|e| {
            Error::config(format!(
                "invalid configuration file {}: {e}",
                path.display()
            ))
        })
    }

    /// Compiles the configuration for the stages in `requested`.
    ///
    /// When `require_configured` is set, every requested stage must carry
    /// options; that is the behaviour of the standalone subcommands and of
    /// `workflow --steps`. Otherwise unconfigured stages are simply dropped.
    pub fn compile(
        self,
        requested: &[Step],
        require_configured: bool,
        input: &FileHeader,
    ) -> Result<Plan> {
        let schema = Schema::from_header(input);
        if self.segment.is_some() {
            // Silence here would be a trap: the table would parse and do nothing.
            return Err(Error::config(
                "[segment] requires the segment command: splitting a read into \
                 fragments changes how many records the output has, so it cannot \
                 run as a workflow step",
            ));
        }
        let Self {
            threads,
            adapter,
            trim,
            filter,
            correction,
            umi,
            output,
            ..
        } = self;

        let mut steps = Vec::new();
        // Resolved before the adapter options are consumed; both stages share it.
        let overlap = overlap_params(adapter.as_ref());
        let umi_stage = compile_umi_stage(umi.as_ref(), requested, require_configured, schema)?;
        if umi_stage.is_some() {
            steps.push(Step::Umi);
        }
        let correction_stage = if requested.contains(&Step::Correct) {
            compile_correction(correction.as_ref(), require_configured, schema)?
        } else {
            None
        };
        if correction_stage.is_some() {
            steps.push(Step::Correct);
        }
        let matcher = adapter.as_ref().map(adapter_params).unwrap_or_default();
        let linked_stage = if requested.contains(&Step::Adapter) {
            compile_linked(adapter.as_ref(), schema, matcher)?
        } else {
            None
        };
        let adapter_stage = if requested.contains(&Step::Adapter) {
            // Linked definitions on their own are enough to configure the stage,
            // so the ordinary matcher may legitimately be absent.
            compile_adapter(
                adapter,
                require_configured && linked_stage.is_none(),
                schema,
            )?
        } else {
            None
        };
        if adapter_stage.is_some() || linked_stage.is_some() {
            steps.push(Step::Adapter);
        }
        let trim_stage = if requested.contains(&Step::Trim) {
            compile_trim(trim, require_configured, schema)?
        } else {
            None
        };
        if trim_stage.is_some() {
            steps.push(Step::Trim);
        }
        let filter_stage = if requested.contains(&Step::Filter) {
            compile_filter(filter, require_configured)?
        } else {
            None
        };
        if filter_stage.is_some() {
            steps.push(Step::Filter);
        }

        let workflow = Workflow::new(
            adapter_stage,
            linked_stage,
            trim_stage,
            filter_stage,
            correction_stage,
            umi_stage,
        )?
        .with_overlap_params(overlap);
        if workflow.needs_quality() && !schema.quality {
            return Err(Error::MissingQuality("quality trimming or filtering"));
        }

        let (pair_policy, mut resolved) = resolve_output(
            output.unwrap_or_default(),
            workflow.can_reject(),
            schema,
            input,
        )?;
        let workflow = workflow.with_pair_policy(pair_policy);
        // The log lives in `[correction]`, not `[output]`, so it is attached here.
        resolved.correction_log = correction.as_ref().and_then(|options| options.log.clone());
        if resolved.correction_log.is_some() && workflow.correction.is_none() {
            return Err(Error::config(
                "--correction-log requires the correction stage (--correction)",
            ));
        }

        let threads = resolve_threads(threads);

        Ok(Plan {
            steps,
            threads,
            workflow,
            output: resolved,
        })
    }
}

impl Config {
    /// Compiles a segmentation run.
    ///
    /// Separate from [`Config::compile`] because segmentation is not a step that
    /// composes with the others: it turns one record into many, so it owns the
    /// whole run. Trimming and filtering still apply, per fragment.
    pub fn compile_segmentation(self, input: &FileHeader) -> Result<Plan> {
        let schema = Schema::from_header(input);
        if schema.paired {
            return Err(Error::config(
                "the segment command requires a single-end input file: splitting a \
                 read into fragments has no defined effect on its mate",
            ));
        }
        let Self {
            threads,
            adapter,
            trim,
            filter,
            segment,
            output,
            ..
        } = self;
        let options = segment.unwrap_or_default();
        let adapter_options = adapter.unwrap_or_default();
        let params = adapter_params(&adapter_options);
        let mut adapters = Vec::new();
        if let Some(sequence) = &adapter_options.r1 {
            adapters.push(Adapter::new("r1", sequence.as_bytes())?);
        }
        if let Some(path) = &adapter_options.fasta {
            adapters.extend(read_adapter_fasta(path)?);
        }
        let segment_stage = crate::segment::SegmentStage::new(
            adapters,
            params,
            options.terminal_fragments.unwrap_or_default(),
            options.min_segment_length.unwrap_or(1),
            options
                .max_segments_per_read
                .unwrap_or(crate::segment::DEFAULT_MAX_SEGMENTS),
        )?;

        let trim_stage = compile_trim(trim, false, schema)?;
        let filter_stage = compile_filter(filter, false)?;
        let mut steps = vec![Step::Segment];
        if trim_stage.is_some() {
            steps.push(Step::Trim);
        }
        if filter_stage.is_some() {
            steps.push(Step::Filter);
        }
        let workflow = Workflow::segmenting(segment_stage, trim_stage, filter_stage);
        if workflow.needs_quality() && !schema.quality {
            return Err(Error::MissingQuality("quality trimming or filtering"));
        }

        let (_, mut resolved) = resolve_output(
            output.unwrap_or_default(),
            workflow.can_reject(),
            schema,
            input,
        )?;
        resolved.segments = options.segments;
        // A header-free input cannot carry the `|segment=` suffix, so the sidecar
        // is the only surviving provenance and is therefore required.
        if resolved.segments.is_none() && !schema.headers {
            return Err(Error::config(
                "--segments is required for an input without headers: the sidecar \
                 is then the only record of where each fragment came from",
            ));
        }
        let threads = resolve_threads(threads);
        Ok(Plan {
            steps,
            threads,
            workflow,
            output: resolved,
        })
    }
}

/// Validates output options and resolves them against the input schema.
fn resolve_output(
    output: OutputOptions,
    can_reject: bool,
    schema: Schema,
    input: &FileHeader,
) -> Result<(PairPolicy, ResolvedOutput)> {
    let pair_policy = output.pair_policy.unwrap_or(PairPolicy::Strict);
    if !can_reject {
        for (path, flag) in [
            (output.failed.as_ref(), "--failed"),
            (output.failed_reasons.as_ref(), "--failed-reasons"),
        ] {
            if path.is_some() {
                return Err(Error::config(format!(
                    "{flag} requires filtering: no record can be rejected \
                     without a filter stage"
                )));
            }
        }
        if pair_policy == PairPolicy::Orphan {
            return Err(Error::config(
                "--pair-policy orphan requires filtering: no record can fail \
                 without a filter stage",
            ));
        }
    }
    if pair_policy == PairPolicy::Orphan {
        if !schema.paired {
            return Err(Error::config(
                "--pair-policy orphan requires a paired input file",
            ));
        }
        if output.orphan_prefix.is_none() {
            return Err(Error::config(
                "--pair-policy orphan requires --orphan-prefix",
            ));
        }
    } else if output.orphan_prefix.is_some() {
        return Err(Error::config(
            "--orphan-prefix requires --pair-policy orphan",
        ));
    }
    let [orphan_r1, orphan_r2] = match &output.orphan_prefix {
        Some(prefix) => {
            let [r1, r2] = orphan_paths(prefix);
            [Some(r1), Some(r2)]
        }
        None => [None, None],
    };
    if output.failed_mode.is_some() && output.failed.is_none() {
        return Err(Error::config("--failed-mode requires --failed"));
    }
    if let Some(level) = output.compression_level
        && level > 22
    {
        return Err(Error::config(format!(
            "--compression-level must be within 0..=22 (got {level})"
        )));
    }
    if let Some(size) = output.block_size
        && size == 0
    {
        return Err(Error::config("--block-size must be greater than zero"));
    }

    Ok((
        pair_policy,
        ResolvedOutput {
            segments: None,
            failed: output.failed,
            correction_log: None,
            failed_reasons: output.failed_reasons,
            failed_mode: output.failed_mode.unwrap_or(FailedMode::Original),
            pair_policy,
            orphan_r1,
            orphan_r2,
            report: output.report,
            report_format: output.report_format.unwrap_or_default(),
            compression_level: output
                .compression_level
                .unwrap_or(input.compression_level as usize),
            block_size: output.block_size.unwrap_or(input.block_size as usize),
        },
    ))
}

/// Resolves the overlap parameters both the adapter and correction stages share.
///
/// The plan forbids correction-specific overlap thresholds, so these come from
/// the adapter options whether or not adapter overlap trimming is enabled.
fn overlap_params(options: Option<&AdapterOptions>) -> crate::overlap::OverlapParams {
    let defaults = crate::overlap::OverlapParams::default();
    let Some(options) = options else {
        return defaults;
    };
    crate::overlap::OverlapParams {
        min_overlap: options
            .paired_overlap_min_overlap
            .unwrap_or(defaults.min_overlap),
        max_error_rate: options.max_error_rate.unwrap_or(defaults.max_error_rate),
    }
}

/// Builds the correction stage, rejecting inputs it cannot be applied to.
fn compile_correction(
    options: Option<&CorrectionOptions>,
    require_configured: bool,
    schema: Schema,
) -> Result<Option<CorrectionStage>> {
    let enabled = options.is_some_and(|options| options.enabled == Some(true));
    if !enabled {
        if require_configured {
            return Err(Error::config(
                "base correction requires --correction (or `[correction] enabled = true`)",
            ));
        }
        // Threshold-only options without the switch are a mistake worth naming.
        if let Some(options) = options
            && (options.donor_quality.is_some()
                || options.recipient_quality.is_some()
                || options.log.is_some()
                || options.log_detail.is_some())
        {
            return Err(Error::config(
                "--donor-quality, --recipient-quality, --correction-log and                  --correction-log-detail require --correction",
            ));
        }
        return Ok(None);
    }
    if !schema.paired {
        return Err(Error::config(
            "base correction requires a paired input file: it needs both mates of              an overlapping pair",
        ));
    }
    if !schema.quality {
        return Err(Error::MissingQuality("base correction"));
    }
    let options = options.expect("enabled implies present");
    let defaults = CorrectionStage::default();
    CorrectionStage {
        donor_quality: options.donor_quality.unwrap_or(defaults.donor_quality),
        recipient_quality: options
            .recipient_quality
            .unwrap_or(defaults.recipient_quality),
        log_detail: options.log_detail.unwrap_or(defaults.log_detail),
    }
    .validate()
    .map(Some)
}

/// Compiles the UMI stage when it was requested.
fn compile_umi_stage(
    umi: Option<&UmiOptions>,
    requested: &[Step],
    require_configured: bool,
    schema: Schema,
) -> Result<Option<UmiStage>> {
    if requested.contains(&Step::Umi) {
        compile_umi(umi, require_configured, schema)
    } else {
        Ok(None)
    }
}

/// Builds the UMI stage, rejecting inputs it cannot be applied to.
fn compile_umi(
    options: Option<&UmiOptions>,
    require_configured: bool,
    schema: Schema,
) -> Result<Option<UmiStage>> {
    let Some(options) = options else {
        if require_configured {
            return Err(Error::config(
                "UMI extraction requires --umi-location (or `[umi] location = ...`)",
            ));
        }
        return Ok(None);
    };
    let Some(location) = options.location else {
        if require_configured {
            return Err(Error::config(
                "UMI extraction requires --umi-location (or `[umi] location = ...`)",
            ));
        }
        return Ok(None);
    };
    if location.needs_paired() && !schema.paired {
        return Err(Error::config(format!(
            "{} UMI requires a paired input file",
            location.name()
        )));
    }
    if !schema.headers {
        return Err(Error::config("UMI processing requires stored read headers"));
    }
    let length = options.length.unwrap_or(0);
    if location.removes_sequence() && length == 0 {
        return Err(Error::config(format!(
            "{} UMI requires --umi-length > 0",
            location.name()
        )));
    }
    Ok(Some(UmiStage {
        location,
        length,
        skip: options.skip.unwrap_or(0),
        prefix: options.prefix.as_deref().unwrap_or("").as_bytes().to_vec(),
        delimiter: options
            .delimiter
            .as_deref()
            .unwrap_or(":")
            .as_bytes()
            .to_vec(),
    }))
}

fn compile_adapter(
    options: Option<AdapterOptions>,
    require_configured: bool,
    schema: Schema,
) -> Result<Option<AdapterStage>> {
    let options = options.unwrap_or_default();
    if options.paired_overlap_min_overlap.is_some() && options.paired_overlap != Some(true) {
        return Err(Error::config(
            "--paired-overlap-min-overlap requires --paired-overlap",
        ));
    }
    if (options.detect_sample_size.is_some() || options.detect_min_support.is_some())
        && options.auto_detect != Some(true)
    {
        return Err(Error::config(
            "--detect-sample-size and --detect-min-support require --auto-detect",
        ));
    }
    if options.is_empty() {
        if require_configured {
            return Err(Error::config(
                "adapter removal requires --adapter-r1, --adapter-r2, --adapter-fasta, \
                 --paired-overlap or --auto-detect",
            ));
        }
        return Ok(None);
    }
    if options.r2.is_some() && !schema.paired {
        return Err(Error::config("--adapter-r2 requires a paired input file"));
    }
    let params = adapter_params(&options);

    let overlap_defaults = crate::overlap::OverlapParams::default();
    let paired_overlap = if options.paired_overlap == Some(true) {
        if !schema.paired {
            return Err(Error::config(
                "--paired-overlap requires a paired input file",
            ));
        }
        Some(
            crate::overlap::OverlapParams {
                min_overlap: options
                    .paired_overlap_min_overlap
                    .unwrap_or(overlap_defaults.min_overlap),
                // The mismatch budget is shared with explicit adapter matching.
                max_error_rate: params.max_error_rate,
            }
            .validate()?,
        )
    } else {
        None
    };

    let mut r1 = Vec::new();
    let mut r2 = Vec::new();
    if let Some(sequence) = &options.r1 {
        r1.push(Adapter::new("r1", sequence.as_bytes())?);
    }
    if let Some(sequence) = &options.r2 {
        r2.push(Adapter::new("r2", sequence.as_bytes())?);
    }
    if let Some(path) = &options.fasta {
        let from_fasta = read_adapter_fasta(path)?;
        r1.extend(from_fasta.iter().cloned());
        if schema.paired {
            r2.extend(from_fasta);
        }
    }
    if options.auto_detect == Some(true) {
        add_detected(&options, schema, &mut r1, &mut r2)?;
    }
    Ok(Some(AdapterStage::new(r1, r2, params, paired_overlap)?))
}

/// Adds the adapters auto-detection recommended. A mate with no confident
/// candidate contributes nothing, and a file where neither mate has one passes
/// through untrimmed: detection outcomes never abort the run.
fn add_detected(
    options: &AdapterOptions,
    schema: Schema,
    r1: &mut Vec<Adapter>,
    r2: &mut Vec<Adapter>,
) -> Result<()> {
    let detection = options
        .detected
        .as_ref()
        .ok_or_else(|| Error::config("--auto-detect was requested, but detection did not run"))?;
    if let Some(sequence) = detection.r1.recommended_sequence.as_deref() {
        r1.push(Adapter::new(
            detection
                .r1
                .recommended_name
                .as_deref()
                .unwrap_or("detected"),
            sequence.as_bytes(),
        )?);
    }
    if schema.paired
        && let Some(mate) = detection.r2.as_ref()
        && let Some(sequence) = mate.recommended_sequence.as_deref()
    {
        r2.push(Adapter::new(
            mate.recommended_name.as_deref().unwrap_or("detected"),
            sequence.as_bytes(),
        )?);
    }
    Ok(())
}

/// Builds the linked-segmentation stage, when definitions were supplied.
fn compile_linked(
    options: Option<&AdapterOptions>,
    schema: Schema,
    params: AdapterParams,
) -> Result<Option<LinkedStage>> {
    let Some(options) = options.filter(|options| options.has_linked()) else {
        return Ok(None);
    };
    let compile = |list: Option<&Vec<LinkedOptions>>| -> Result<Vec<LinkedAdapter>> {
        list.map_or_else(
            || Ok(Vec::new()),
            |definitions| {
                definitions
                    .iter()
                    .enumerate()
                    .map(|(position, definition)| definition.compile(position))
                    .collect()
            },
        )
    };
    let r2 = compile(options.linked_r2.as_ref())?;
    if !r2.is_empty() && !schema.paired {
        return Err(Error::config(
            "--linked-5p-r2 and --linked-3p-r2 require a paired input file",
        ));
    }
    Ok(Some(LinkedStage::new(
        compile(options.linked_r1.as_ref())?,
        r2,
        options.linked_unmatched.unwrap_or(Unmatched::Continue),
        params,
    )?))
}

/// Resolves the matcher thresholds from raw adapter options.
pub(crate) fn adapter_params(options: &AdapterOptions) -> AdapterParams {
    let defaults = AdapterParams::default();
    AdapterParams {
        min_overlap: options.min_overlap.unwrap_or(defaults.min_overlap),
        max_error_rate: options.max_error_rate.unwrap_or(defaults.max_error_rate),
        max_errors: options.max_errors,
        allow_indels: options.allow_indels == Some(true),
    }
}

/// Resolves detection thresholds from raw adapter options.
pub(crate) fn detection_params(options: &AdapterOptions) -> Result<SniffParams> {
    let defaults = SniffParams::default();
    let mut gates = defaults.gates;
    if let Some(support) = options.detect_min_support {
        if !(0.0..=1.0).contains(&support) {
            return Err(Error::config(format!(
                "--detect-min-support must be within 0.0..=1.0 (got {support})"
            )));
        }
        gates.min_support_fraction = support;
    }
    SniffParams {
        sample_size: options
            .detect_sample_size
            .map_or(defaults.sample_size, |size| size as u64),
        matcher: adapter_params(options),
        gates,
        ..defaults
    }
    .validate()
}

// `TrimOptions` is a large `Copy` struct, but compilation happens once per run.
#[allow(clippy::large_types_passed_by_value)]
fn compile_trim(
    options: Option<TrimOptions>,
    require_configured: bool,
    schema: Schema,
) -> Result<Option<TrimStage>> {
    let options = options.unwrap_or_default();
    if options.uses_r2_only_options() && !schema.paired {
        return Err(Error::config(
            "--front-r2, --tail-r2 and --max-length-r2 require a paired input file",
        ));
    }
    match options.compile() {
        Ok(stage) => Ok(Some(stage)),
        Err(error) => {
            // A completely unconfigured stage is only an error when the user
            // asked for it explicitly. Any other problem is always reported.
            let unconfigured = matches!(&error, Error::InvalidConfiguration(message)
                if message.starts_with("trimming requires"));
            if unconfigured && !require_configured {
                Ok(None)
            } else {
                Err(error)
            }
        }
    }
}

fn compile_filter(
    options: Option<FilterOptions>,
    require_configured: bool,
) -> Result<Option<FilterStage>> {
    let options = options.unwrap_or_default();
    match options.compile() {
        Ok(stage) => Ok(Some(stage)),
        Err(error) => {
            let unconfigured = matches!(&error, Error::InvalidConfiguration(message)
                if message.starts_with("filtering requires"));
            if unconfigured && !require_configured {
                Ok(None)
            } else {
                Err(error)
            }
        }
    }
}

#[cfg(test)]
#[allow(clippy::float_cmp)]
mod tests {
    use super::*;
    use crate::trim::TrimOp;

    fn header(paired: bool, quality: bool) -> FileHeader {
        Schema {
            paired,
            quality,
            headers: true,
            flags: false,
        }
        .to_file_header(4, 1 << 20)
    }

    fn adapter_config() -> Config {
        Config {
            adapter: Some(AdapterOptions {
                r1: Some("AGATCGGAAGAGCACACGTCTGAACTCCAGTCA".to_string()),
                ..AdapterOptions::default()
            }),
            ..Config::default()
        }
    }

    #[test]
    fn cli_values_override_file_values_field_by_field() {
        let file: Config = toml::from_str(
            r#"
            threads = 4
            [adapter]
            r1 = "ACGTACGTACGT"
            min_overlap = 6
            [trim.quality_tail]
            minimum_phred = 20
            window = 4
            "#,
        )
        .unwrap();
        let cli = Config {
            adapter: Some(AdapterOptions {
                min_overlap: Some(10),
                ..AdapterOptions::default()
            }),
            trim: Some(TrimOptions {
                quality_tail: Some(QualityCutOptions {
                    minimum_phred: Some(30),
                    window: None,
                }),
                ..TrimOptions::default()
            }),
            ..Config::default()
        };
        let merged = cli.merge(file);
        let adapter = merged.adapter.as_ref().unwrap();
        assert_eq!(adapter.min_overlap, Some(10), "CLI wins");
        assert_eq!(
            adapter.r1.as_deref(),
            Some("ACGTACGTACGT"),
            "file value survives"
        );
        let quality_tail = merged.trim.unwrap().quality_tail.unwrap();
        assert_eq!(quality_tail.minimum_phred, Some(30));
        assert_eq!(
            quality_tail.window,
            Some(4),
            "unset CLI field keeps the file value"
        );
        assert_eq!(merged.threads, Some(4));
    }

    #[cfg(unix)]
    #[test]
    fn typed_merge_preserves_non_utf8_paths() {
        use std::ffi::OsString;
        use std::os::unix::ffi::OsStringExt as _;

        let report = PathBuf::from(OsString::from_vec(b"report-\xff.json".to_vec()));
        let cli = Config {
            output: Some(OutputOptions {
                report: Some(report.clone()),
                ..OutputOptions::default()
            }),
            ..Config::default()
        };
        let file = Config {
            output: Some(OutputOptions {
                compression_level: Some(7),
                ..OutputOptions::default()
            }),
            ..Config::default()
        };

        let output = cli.merge(file).output.unwrap();
        assert_eq!(output.report, Some(report));
        assert_eq!(output.compression_level, Some(7));
    }

    #[test]
    fn the_documented_toml_schema_parses() {
        let config: Config = toml::from_str(
            r#"
            threads = 8

            [adapter]
            r1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
            r2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
            min_overlap = 8
            max_error_rate = 0.10

            [trim.quality_tail]
            minimum_phred = 20
            window = 4

            [trim]
            terminal_n = true

            [filter]
            min_length = 30
            max_n = 5
            qualified_quality = 15
            max_unqualified_fraction = 0.40

            [output]
            failed = "rejected.cbq"
            report = "clean.bqc.json"
            "#,
        )
        .unwrap();
        let plan = config
            .compile(&Step::ALL, false, &header(true, true))
            .unwrap();
        assert_eq!(plan.threads, 8);
        assert_eq!(plan.steps, vec![Step::Adapter, Step::Trim, Step::Filter]);
        let trim = plan.workflow.trim.unwrap();
        assert!(trim.r1.terminal_n && trim.r2.terminal_n);
        assert_eq!(trim.r1.quality_tail.unwrap().window, 4);
        assert_eq!(plan.workflow.filter.unwrap().min_length, Some(30));
        assert_eq!(plan.output.failed.unwrap(), PathBuf::from("rejected.cbq"));
    }

    #[test]
    fn the_documented_optional_feature_keys_parse() {
        // Mirrors the optional-feature TOML block in the README, so the
        // documentation cannot drift from the accepted schema.
        let config: Config = toml::from_str(
            r#"
            [adapter]
            r1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
            allow_indels = true
            paired_overlap = true
            paired_overlap_min_overlap = 30
            auto_detect = true
            detect_sample_size = 100000
            detect_min_support = 0.01

            [filter]
            min_length = 30

            [output]
            pair_policy = "orphan"
            orphan_prefix = "surviving"
            "#,
        )
        .unwrap();
        let adapter = config.adapter.as_ref().unwrap();
        assert_eq!(adapter.allow_indels, Some(true));
        assert_eq!(adapter.paired_overlap, Some(true));
        assert_eq!(adapter.paired_overlap_min_overlap, Some(30));
        assert_eq!(adapter.auto_detect, Some(true));
        assert_eq!(adapter.detect_sample_size, Some(100_000));
        assert_eq!(adapter.detect_min_support, Some(0.01));
        let output = config.output.as_ref().unwrap();
        assert_eq!(output.pair_policy, Some(PairPolicy::Orphan));
        assert_eq!(output.orphan_prefix, Some(PathBuf::from("surviving")));
        // A detection result is runtime state and must not be injectable.
        assert!(
            toml::from_str::<Config>("[adapter]\ndetected = {}").is_err(),
            "detection results are not configuration"
        );
    }

    #[test]
    fn unknown_configuration_keys_are_rejected() {
        assert!(toml::from_str::<Config>("threadz = 4").is_err());
        assert!(toml::from_str::<Config>("[adapter]\nmode = \"three-prime\"").is_err());
    }

    #[test]
    fn output_defaults_are_inherited_from_the_input_header() {
        let plan = adapter_config()
            .compile(&[Step::Adapter], true, &header(false, true))
            .unwrap();
        assert_eq!(plan.output.compression_level, 4);
        assert_eq!(plan.output.block_size, 1 << 20);
        assert_eq!(plan.output.failed_mode, FailedMode::Original);
        assert_eq!(plan.output.report_format, ReportFormat::Json);
    }

    #[test]
    fn standalone_commands_require_their_own_configuration() {
        let empty = Config::default();
        let err = empty
            .clone()
            .compile(&[Step::Adapter], true, &header(false, true))
            .unwrap_err();
        assert!(
            format!("{err}").contains("adapter removal requires"),
            "{err}"
        );
        let err = empty
            .clone()
            .compile(&[Step::Trim], true, &header(false, true))
            .unwrap_err();
        assert!(format!("{err}").contains("trimming requires"), "{err}");
        let err = empty
            .compile(&[Step::Filter], true, &header(false, true))
            .unwrap_err();
        assert!(format!("{err}").contains("filtering requires"), "{err}");
    }

    #[test]
    fn workflow_drops_unconfigured_stages_but_rejects_an_empty_workflow() {
        let plan = adapter_config()
            .compile(&Step::ALL, false, &header(false, true))
            .unwrap();
        assert_eq!(plan.steps, vec![Step::Adapter]);
        assert!(plan.workflow.trim.is_none() && plan.workflow.filter.is_none());

        let err = Config::default()
            .compile(&Step::ALL, false, &header(false, true))
            .unwrap_err();
        assert!(
            format!("{err}").contains("no operation configured"),
            "{err}"
        );
    }

    #[test]
    fn requested_steps_ignore_configuration_for_other_stages() {
        let config = Config {
            filter: Some(FilterOptions {
                min_length: Some(30),
                ..FilterOptions::default()
            }),
            ..adapter_config()
        };
        let plan = config
            .compile(&[Step::Adapter], false, &header(false, true))
            .unwrap();
        assert_eq!(plan.steps, vec![Step::Adapter]);
        assert!(plan.workflow.filter.is_none());
    }

    #[test]
    fn quality_operations_require_a_quality_column() {
        let config = Config {
            filter: Some(FilterOptions {
                min_mean_quality: Some(20),
                ..FilterOptions::default()
            }),
            ..Config::default()
        };
        let err = config
            .compile(&[Step::Filter], true, &header(false, false))
            .unwrap_err();
        assert!(matches!(err, Error::MissingQuality(_)), "{err}");
    }

    #[test]
    fn mate_two_options_require_a_paired_input() {
        let config = Config {
            trim: Some(TrimOptions {
                front_r2: Some(4),
                ..TrimOptions::default()
            }),
            ..Config::default()
        };
        let err = config
            .compile(&[Step::Trim], true, &header(false, true))
            .unwrap_err();
        assert!(format!("{err}").contains("require a paired input"), "{err}");

        let config = Config {
            adapter: Some(AdapterOptions {
                r2: Some("ACGTACGTACGT".to_string()),
                ..AdapterOptions::default()
            }),
            ..Config::default()
        };
        let err = config
            .compile(&[Step::Adapter], true, &header(false, true))
            .unwrap_err();
        assert!(
            format!("{err}").contains("--adapter-r2 requires a paired input"),
            "{err}"
        );
    }

    #[test]
    fn failed_outputs_require_a_filter_stage() {
        let config = Config {
            output: Some(OutputOptions {
                failed: Some(PathBuf::from("rejected.cbq")),
                ..OutputOptions::default()
            }),
            ..adapter_config()
        };
        let err = config
            .compile(&[Step::Adapter], true, &header(false, true))
            .unwrap_err();
        assert!(
            format!("{err}").contains("--failed requires filtering"),
            "{err}"
        );
    }

    #[test]
    fn mutually_exclusive_quality_cuts_are_rejected() {
        let config = Config {
            trim: Some(TrimOptions {
                quality_right: Some(QualityCutOptions {
                    minimum_phred: Some(20),
                    window: None,
                }),
                quality_tail: Some(QualityCutOptions {
                    minimum_phred: Some(20),
                    window: None,
                }),
                ..TrimOptions::default()
            }),
            ..Config::default()
        };
        let err = config
            .compile(&[Step::Trim], true, &header(false, true))
            .unwrap_err();
        assert!(format!("{err}").contains("mutually exclusive"), "{err}");
    }

    #[test]
    fn dependent_options_require_their_primary_flag() {
        let config = Config {
            trim: Some(TrimOptions {
                quality_tail: Some(QualityCutOptions {
                    minimum_phred: None,
                    window: Some(4),
                }),
                ..TrimOptions::default()
            }),
            ..Config::default()
        };
        let err = config
            .compile(&[Step::Trim], true, &header(false, true))
            .unwrap_err();
        assert!(
            format!("{err}").contains("--quality-tail-window requires"),
            "{err}"
        );

        let config = Config {
            filter: Some(FilterOptions {
                qualified_quality: Some(20),
                ..FilterOptions::default()
            }),
            ..Config::default()
        };
        let err = config
            .compile(&[Step::Filter], true, &header(false, true))
            .unwrap_err();
        assert!(
            format!("{err}").contains("--qualified-quality requires"),
            "{err}"
        );
    }

    #[test]
    fn trim_defaults_are_applied_when_an_operation_is_enabled() {
        let config = Config {
            trim: Some(TrimOptions {
                quality_tail: Some(QualityCutOptions {
                    minimum_phred: Some(20),
                    window: None,
                }),
                poly_g: Some(PolyOptions::default()),
                ..TrimOptions::default()
            }),
            ..Config::default()
        };
        let plan = config
            .compile(&[Step::Trim], true, &header(false, true))
            .unwrap();
        let trim = plan.workflow.trim.unwrap();
        assert_eq!(trim.r1.quality_tail.unwrap().window, DEFAULT_QUALITY_WINDOW);
        assert_eq!(trim.r1.poly_g.unwrap().min_length, 10);
        assert!((trim.r1.poly_g.unwrap().max_mismatch_rate - 0.10).abs() < 1e-12);
        assert!(!trim.r1.is_noop());
        // Shared operations apply to both mates.
        assert_eq!(trim.r2.poly_g.unwrap().min_length, 10);
        assert_eq!(trim.r2.quality_tail.unwrap().minimum_phred, 20);
    }

    #[test]
    fn thread_count_zero_resolves_to_available_parallelism() {
        let config = Config {
            threads: Some(0),
            ..adapter_config()
        };
        let plan = config
            .compile(&[Step::Adapter], true, &header(false, true))
            .unwrap();
        assert!(plan.threads >= 1);
    }

    #[test]
    fn compiled_plans_serialize_the_resolved_configuration() {
        let config = Config {
            trim: Some(TrimOptions {
                terminal_n: Some(true),
                poly_g: Some(PolyOptions::default()),
                ..TrimOptions::default()
            }),
            ..adapter_config()
        };
        let plan = config
            .compile(&Step::ALL, false, &header(false, true))
            .unwrap();
        let json = serde_json::to_value(&plan).unwrap();
        assert!(json["threads"].is_number());
        assert_eq!(json["workflow"]["pair_policy"], "strict");
        assert_eq!(json["workflow"]["adapter"]["params"]["min_overlap"], 8);
        assert_eq!(json["workflow"]["adapter"]["r1"][0]["name"], "r1");
        assert_eq!(
            json["workflow"]["adapter"]["r1"][0]["sequence"],
            "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
        );
        assert_eq!(json["workflow"]["trim"]["r1"]["terminal_n"], true);
        assert_eq!(json["workflow"]["trim"]["r1"]["poly_g"]["min_length"], 10);
        assert_eq!(json["output"]["failed_mode"], "original");
        // A disabled stage is explicit rather than absent.
        assert!(json["workflow"]["filter"].is_null());
        // Trim operations that were not requested are visible as null/zero.
        assert_eq!(json["workflow"]["trim"]["r1"]["front"], 0);
        assert!(json["workflow"]["trim"]["r1"]["quality_tail"].is_null());
        assert_eq!(TrimOp::PolyG.name(), "poly_g");
    }
}
