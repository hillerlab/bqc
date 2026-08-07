// Copyright (c) 2026 The Hiller Lab at the Senckenberg Gesellschaft für Naturforschung
// Distributed under the terms of the GNU General Public License, Version 3.0.

//! Structured reports and the console summary.
//!
//! `bqc` reports operational statistics — what it did to the data — not
//! dataset-level quality distributions. Per-position quality tables, GC
//! distributions and duplication estimates belong to `bqtools qc`.

use std::path::{Path, PathBuf};
use std::time::Duration;

use serde::{Serialize, Serializer};

use crate::config::{Plan, ReportFormat};
use crate::error::Result;
use crate::filter::{FilterReason, REASONS};
use crate::io::Schema;
use crate::process::Mate;
use crate::stats::Stats;
use crate::trim::{TrimOp, TrimStage};

/// Version of the `binseq` crate this build links against.
///
/// Kept in sync with `Cargo.toml` by `binseq_version_matches_manifest`.
pub const BINSEQ_VERSION: &str = "0.9.4";

/// Serializes a byte string as UTF-8 text (adapters are ASCII by construction).
pub fn serialize_bytes<S: Serializer>(
    bytes: &[u8],
    serializer: S,
) -> std::result::Result<S::Ok, S::Error> {
    serializer.serialize_str(&String::from_utf8_lossy(bytes))
}

#[derive(Debug, Serialize)]
pub struct InputReport {
    pub path: PathBuf,
    pub bytes: u64,
    pub records: u64,
    pub blocks: usize,
    pub schema: Schema,
    pub cbq_version: u8,
    pub compression_level: u64,
    pub block_size: u64,
}

#[derive(Debug, Serialize)]
pub struct AdapterEntry {
    pub name: String,
    pub mate: &'static str,
    pub hits: u64,
}

#[derive(Debug, Serialize)]
pub struct AdapterReport {
    pub r1_reads_trimmed: u64,
    pub r2_reads_trimmed: u64,
    pub r1_bases_removed: u64,
    pub r2_bases_removed: u64,
    pub r1_overlap_reads_trimmed: u64,
    pub r2_overlap_reads_trimmed: u64,
    pub r1_overlap_bases_removed: u64,
    pub r2_overlap_bases_removed: u64,
    pub per_adapter: Vec<AdapterEntry>,
    /// Auto-detection evidence, when `--auto-detect` was enabled.
    pub detection: Option<crate::sniff::adapters::AdapterSniff>,
}

#[derive(Debug, Serialize)]
pub struct TrimEntry {
    pub operation: &'static str,
    pub reads_trimmed: u64,
    pub bases_removed: u64,
}

#[derive(Debug, Serialize)]
pub struct ReasonEntry {
    pub reasons: String,
    pub records: u64,
}

#[derive(Debug, Serialize)]
pub struct FilterReport {
    pub accepted_records: u64,
    pub rejected_records: u64,
    pub orphan_r1_records: u64,
    pub orphan_r2_records: u64,
    pub r1_only_failed: u64,
    pub r2_only_failed: u64,
    pub both_failed: u64,
    pub per_reason: Vec<ReasonEntry>,
    pub per_combination: Vec<ReasonEntry>,
}

#[derive(Debug, Serialize)]
pub struct CountsReport {
    pub records_in: u64,
    pub records_out: u64,
    pub records_rejected: u64,
    pub r1_bases_in: u64,
    pub r2_bases_in: u64,
    pub r1_bases_out: u64,
    pub r2_bases_out: u64,
    /// Bases written to the orphan outputs, which are separate files.
    pub r1_bases_orphaned: u64,
    pub r2_bases_orphaned: u64,
}

#[derive(Debug, Serialize)]
pub struct PerformanceReport {
    pub elapsed_seconds: f64,
    pub records_per_second: f64,
    pub bases_per_second: f64,
}

#[derive(Debug, Serialize)]
pub struct OutputReport {
    pub accepted: PathBuf,
    pub failed: Option<PathBuf>,
    pub failed_reasons: Option<PathBuf>,
    pub orphan_r1: Option<PathBuf>,
    pub orphan_r2: Option<PathBuf>,
    pub report: Option<PathBuf>,
}

/// Correction statistics. Present exactly when the correction stage ran.
///
/// Counts describe corrections applied to the read, before trimming: a corrected
/// base that a later stage trims away is still counted, and
/// `corrected_pairs_by_disposition` says where those pairs ended up.
#[derive(Debug, Serialize)]
pub struct CorrectionReport {
    pub pairs_examined: u64,
    pub pairs_with_overlap: u64,
    pub pairs_with_mismatches: u64,
    pub overlap_mismatches: u64,
    pub corrected_pairs: u64,
    pub corrected_pairs_both_mates: u64,
    pub corrected_reads: u64,
    pub corrected_r1_reads: u64,
    pub corrected_r2_reads: u64,
    pub corrected_bases: u64,
    /// Bases corrected in R1, which is also the R2-to-R1 direction count.
    pub corrected_r1_bases: u64,
    /// Bases corrected in R2, which is also the R1-to-R2 direction count.
    pub corrected_r2_bases: u64,
    pub unresolved_mismatches: u64,
    pub noncanonical_donors_skipped: u64,
    /// Corrections per pair, over pairs with at least one correction.
    pub corrections_per_pair: Vec<CountEntry>,
    /// Original-to-corrected substitutions, `A C G T N` in both directions.
    pub substitutions: Vec<SubstitutionEntry>,
    pub corrected_pairs_by_disposition: Vec<DispositionEntry>,
}

#[derive(Debug, Serialize)]
pub struct CountEntry {
    pub corrections: u32,
    pub pairs: u64,
}

#[derive(Debug, Serialize)]
pub struct SubstitutionEntry {
    pub original: char,
    pub corrected: char,
    pub bases: u64,
}

#[derive(Debug, Serialize)]
pub struct DispositionEntry {
    pub disposition: String,
    pub pairs: u64,
}

/// Linked-segmentation statistics. Present exactly when the stage ran.
#[derive(Debug, Serialize)]
pub struct LinkedReport {
    pub both_adapters: u64,
    pub five_prime_only: u64,
    pub three_prime_only: u64,
    pub invalid_order: u64,
    pub insert_too_short: u64,
    pub neither_adapter: u64,
    pub leading_bases_removed: u64,
    pub trailing_bases_removed: u64,
    /// Per mate, then per definition.
    pub per_definition: Vec<LinkedEntry>,
    pub unmatched_policy: crate::linked::Unmatched,
}

#[derive(Debug, Serialize)]
pub struct LinkedEntry {
    pub name: String,
    pub mate: &'static str,
    pub matches: u64,
    pub both_adapters_required: bool,
}

/// Segmentation statistics. Present exactly when the `segment` command ran.
///
/// Source records and fragments are counted separately throughout: one source
/// record can become any number of fragments, so a single "records" number would
/// be ambiguous. The top-level `counts` are output records — fragments.
#[derive(Debug, Serialize)]
pub struct SegmentReport {
    pub source_records: u64,
    /// Source records with at least one accepted delimiter.
    pub source_records_split: u64,
    pub delimiters_accepted: u64,
    /// Candidate matches dropped for overlapping an accepted delimiter.
    pub candidates_suppressed: u64,
    pub fragments_emitted: u64,
    /// Fragments with an adapter on both sides.
    pub internal_fragments: u64,
    /// Fragments touching a source-read end.
    pub terminal_fragments: u64,
    pub discarded_empty: u64,
    pub discarded_too_short: u64,
    pub discarded_terminal: u64,
    pub discarded_over_limit: u64,
    pub fragments_accepted: u64,
    pub fragments_rejected: u64,
    /// Most fragments any one source record produced.
    pub max_fragments_per_source: u64,
    /// Fragments emitted per source record.
    pub fragments_per_source: Vec<FragmentCountEntry>,
    /// Length summary over emitted fragments, as cut, before trimming.
    pub fragment_length: FragmentLengthSummary,
    /// Delimiters matched per declared adapter.
    pub per_adapter: Vec<AdapterEntry>,
    pub terminal_policy: crate::segment::Terminal,
    pub min_segment_length: u64,
    pub max_segments_per_read: u64,
}

/// Fragment lengths as segmentation produced them, before any trimming.
#[derive(Debug, Serialize)]
pub struct FragmentLengthSummary {
    pub shortest: u64,
    pub longest: u64,
    pub mean: f64,
}

#[derive(Debug, Serialize)]
pub struct FragmentCountEntry {
    pub fragments: u32,
    pub source_records: u64,
}

/// The complete structured report.
#[derive(Debug, Serialize)]
pub struct Report {
    pub tool: &'static str,
    pub version: &'static str,
    pub binseq_version: &'static str,
    pub command: String,
    pub input: InputReport,
    pub configuration: ReportConfiguration,
    pub outputs: OutputReport,
    pub counts: CountsReport,
    pub adapter: Option<AdapterReport>,
    pub trim: Option<Vec<TrimEntry>>,
    pub filter: Option<FilterReport>,
    pub linked: Option<LinkedReport>,
    pub segment: Option<SegmentReport>,
    pub correction: Option<CorrectionReport>,
    pub performance: PerformanceReport,
}

/// The resolved configuration, including the stage order actually executed.
#[derive(Debug, Serialize)]
pub struct ReportConfiguration {
    pub stage_order: Vec<&'static str>,
    pub span: Option<[u64; 2]>,
    #[serde(flatten)]
    pub plan: serde_json::Value,
}

/// Everything a finished run contributes to its report.
pub struct RunReport<'a> {
    /// Subcommand name.
    pub command: &'a str,
    pub input: &'a crate::io::CbqInput,
    pub plan: &'a Plan,
    /// Record range processed, when restricted.
    pub span: Option<std::ops::Range<u64>>,
    /// Path of the accepted output.
    pub accepted: &'a Path,
    pub stats: &'a Stats,
    pub elapsed: Duration,
    pub detection: Option<&'a crate::sniff::adapters::AdapterSniff>,
}

impl Report {
    /// Assembles a report from a finished run.
    #[must_use]
    pub fn new(run: RunReport<'_>) -> Self {
        let RunReport {
            command,
            input,
            plan,
            span,
            accepted,
            stats,
            elapsed,
            detection,
        } = run;
        let header = input.header();
        let seconds = elapsed.as_secs_f64();
        let per_second = |count: u64| {
            if seconds > 0.0 {
                count as f64 / seconds
            } else {
                0.0
            }
        };

        Self {
            tool: "bqc",
            version: env!("CARGO_PKG_VERSION"),
            binseq_version: BINSEQ_VERSION,
            command: command.to_string(),
            input: InputReport {
                path: input.path().to_path_buf(),
                bytes: std::fs::metadata(input.path()).map_or(0, |m| m.len()),
                records: input.num_records(),
                blocks: input.blocks().len(),
                schema: input.schema(),
                cbq_version: header.version,
                compression_level: header.compression_level,
                block_size: header.block_size,
            },
            configuration: ReportConfiguration {
                stage_order: plan.workflow.stage_order(),
                span: span.map(|span| [span.start, span.end]),
                plan: serde_json::to_value(plan).unwrap_or(serde_json::Value::Null),
            },
            outputs: OutputReport {
                accepted: accepted.to_path_buf(),
                failed: plan.output.failed.clone(),
                failed_reasons: plan.output.failed_reasons.clone(),
                orphan_r1: plan.output.orphan_r1.clone(),
                orphan_r2: plan.output.orphan_r2.clone(),
                report: plan.output.report.clone(),
            },
            counts: CountsReport {
                records_in: stats.records_in,
                records_out: stats.records_out,
                records_rejected: stats.records_rejected,
                r1_bases_in: stats.bases_in[0],
                r2_bases_in: stats.bases_in[1],
                r1_bases_out: stats.bases_out[0],
                r2_bases_out: stats.bases_out[1],
                r1_bases_orphaned: stats.bases_orphaned[0],
                r2_bases_orphaned: stats.bases_orphaned[1],
            },
            adapter: plan
                .workflow
                .adapter
                .as_ref()
                .map(|stage| adapter_report(stage, stats, detection)),
            trim: plan
                .workflow
                .trim
                .as_ref()
                .map(|stage| trim_entries(stage, stats)),
            filter: plan.workflow.filter.as_ref().map(|_| filter_report(stats)),
            linked: plan
                .workflow
                .linked
                .as_ref()
                .map(|stage| linked_report(stage, stats)),
            segment: plan
                .workflow
                .segment
                .as_ref()
                .map(|stage| segment_report(stage, stats)),
            correction: plan
                .workflow
                .correction
                .as_ref()
                .map(|_| correction_report(stats)),
            performance: PerformanceReport {
                elapsed_seconds: seconds,
                records_per_second: per_second(stats.records_in),
                bases_per_second: per_second(stats.total_bases_in()),
            },
        }
    }

    /// Renders the report in the requested format.
    pub fn render(&self, format: ReportFormat) -> Result<String> {
        let value = serde_json::to_value(self).unwrap_or(serde_json::Value::Null);
        Ok(match format {
            ReportFormat::Json => {
                let mut text = serde_json::to_string_pretty(&value).unwrap_or_default();
                text.push('\n');
                text
            }
            ReportFormat::Tsv => {
                let mut rows = Vec::new();
                flatten(&value, String::new(), &mut rows);
                let mut text = String::from("key\tvalue\n");
                for (key, value) in rows {
                    text.push_str(&key);
                    text.push('\t');
                    text.push_str(&value);
                    text.push('\n');
                }
                text
            }
        })
    }

    /// Writes the human-readable summary to `out`.
    pub fn write_summary<W: std::io::Write>(&self, out: &mut W) -> std::io::Result<()> {
        let unit = if self.input.schema.paired {
            "pairs"
        } else {
            "reads"
        };
        self.write_input(out, unit)?;
        self.write_stages(out, unit)?;
        self.write_linked(out)?;
        self.write_segment(out)?;
        self.write_correction(out)?;
        self.write_totals(out)
    }

    fn write_input<W: std::io::Write>(&self, out: &mut W, unit: &str) -> std::io::Result<()> {
        writeln!(out, "Input")?;
        writeln!(
            out,
            "  records:                {:>16} {unit}",
            commas(self.counts.records_in)
        )?;
        writeln!(
            out,
            "  bases R1:               {:>16}",
            commas(self.counts.r1_bases_in)
        )?;
        if self.input.schema.paired {
            writeln!(
                out,
                "  bases R2:               {:>16}",
                commas(self.counts.r2_bases_in)
            )?;
        }
        Ok(())
    }

    /// Summarises what auto-detection concluded, and on what evidence.
    fn write_detection<W: std::io::Write>(
        out: &mut W,
        detection: &crate::sniff::adapters::AdapterSniff,
    ) -> std::io::Result<()> {
        for mate in std::iter::once(&detection.r1).chain(detection.r2.as_ref()) {
            match (&mate.recommended_sequence, mate.candidates.first()) {
                (Some(sequence), Some(leader)) => writeln!(
                    out,
                    "  detected {}: {sequence} ({}, {:.2}% of reads)",
                    mate.mate,
                    leader
                        .evidence_sources
                        .iter()
                        .map(|source| source.name())
                        .collect::<Vec<_>>()
                        .join("+"),
                    leader.support_fraction * 100.0,
                )?,
                _ => writeln!(
                    out,
                    "  detected {}: none ({})",
                    mate.mate,
                    mate.decision.name()
                )?,
            }
        }
        Ok(())
    }

    fn write_stages<W: std::io::Write>(&self, out: &mut W, unit: &str) -> std::io::Result<()> {
        if let Some(adapter) = &self.adapter {
            writeln!(out, "\nAdapter")?;
            writeln!(
                out,
                "  R1 reads trimmed:       {:>16}",
                commas(adapter.r1_reads_trimmed)
            )?;
            if self.input.schema.paired {
                writeln!(
                    out,
                    "  R2 reads trimmed:       {:>16}",
                    commas(adapter.r2_reads_trimmed)
                )?;
            }
            writeln!(
                out,
                "  bases removed:          {:>16}",
                commas(adapter.r1_bases_removed + adapter.r2_bases_removed)
            )?;
            let overlap_reads = adapter.r1_overlap_reads_trimmed + adapter.r2_overlap_reads_trimmed;
            if overlap_reads > 0 {
                writeln!(
                    out,
                    "  overlap-trimmed reads:  {:>16}",
                    commas(overlap_reads)
                )?;
                writeln!(
                    out,
                    "  overlap bases removed:  {:>16}",
                    commas(adapter.r1_overlap_bases_removed + adapter.r2_overlap_bases_removed)
                )?;
            }
            for entry in &adapter.per_adapter {
                let label = format!("{} hits ({}):", entry.name, entry.mate);
                writeln!(out, "  {label:<23} {:>16}", commas(entry.hits))?;
            }
            if let Some(detection) = &adapter.detection {
                Self::write_detection(out, detection)?;
            }
        }

        if let Some(trim) = &self.trim {
            writeln!(out, "\nTrim")?;
            for entry in trim {
                writeln!(
                    out,
                    "  {:<22} {:>16} reads, {:>16} bases",
                    entry.operation,
                    commas(entry.reads_trimmed),
                    commas(entry.bases_removed)
                )?;
            }
        }

        if let Some(filter) = &self.filter {
            writeln!(out, "\nFilter")?;
            writeln!(
                out,
                "  accepted {unit}: {:>21}",
                commas(filter.accepted_records)
            )?;
            writeln!(
                out,
                "  rejected {unit}: {:>21}",
                commas(filter.rejected_records)
            )?;
            if filter.orphan_r1_records > 0 {
                writeln!(
                    out,
                    "  orphan R1 reads:      {:>16}",
                    commas(filter.orphan_r1_records)
                )?;
            }
            if filter.orphan_r2_records > 0 {
                writeln!(
                    out,
                    "  orphan R2 reads:      {:>16}",
                    commas(filter.orphan_r2_records)
                )?;
            }
            for entry in &filter.per_reason {
                if entry.records > 0 {
                    writeln!(out, "  {:<22} {:>16}", entry.reasons, commas(entry.records))?;
                }
            }
        }
        Ok(())
    }

    fn write_linked<W: std::io::Write>(&self, out: &mut W) -> std::io::Result<()> {
        let Some(linked) = &self.linked else {
            return Ok(());
        };
        writeln!(out, "\nLinked segmentation")?;
        for (label, value) in [
            ("both adapters:", linked.both_adapters),
            ("5' adapter only:", linked.five_prime_only),
            ("3' adapter only:", linked.three_prime_only),
            ("invalid order:", linked.invalid_order),
            ("insert too short:", linked.insert_too_short),
            ("neither adapter:", linked.neither_adapter),
        ] {
            writeln!(out, "  {label:<22} {:>16}", commas(value))?;
        }
        writeln!(
            out,
            "  bases removed:         {:>16} before, {} after",
            commas(linked.leading_bases_removed),
            commas(linked.trailing_bases_removed)
        )?;
        for entry in &linked.per_definition {
            let label = format!("{} matches ({}):", entry.name, entry.mate);
            writeln!(out, "  {label:<22} {:>16}", commas(entry.matches))?;
        }
        Ok(())
    }

    fn write_segment<W: std::io::Write>(&self, out: &mut W) -> std::io::Result<()> {
        let Some(segment) = &self.segment else {
            return Ok(());
        };
        writeln!(out, "\nSegmentation")?;
        for (label, value) in [
            ("source records:", segment.source_records),
            ("records split:", segment.source_records_split),
            ("delimiters:", segment.delimiters_accepted),
            ("suppressed candidates:", segment.candidates_suppressed),
            ("fragments emitted:", segment.fragments_emitted),
            ("  internal:", segment.internal_fragments),
            ("  terminal:", segment.terminal_fragments),
            ("discarded empty:", segment.discarded_empty),
            ("discarded too short:", segment.discarded_too_short),
            ("discarded terminal:", segment.discarded_terminal),
            ("over segment limit:", segment.discarded_over_limit),
            ("most from one read:", segment.max_fragments_per_source),
        ] {
            writeln!(out, "  {label:<22} {:>16}", commas(value))?;
        }
        writeln!(
            out,
            "  fragment length:       {:>16} shortest, {} longest, {:.1} mean",
            commas(segment.fragment_length.shortest),
            commas(segment.fragment_length.longest),
            segment.fragment_length.mean,
        )?;
        for entry in &segment.per_adapter {
            let label = format!("{} delimiters:", entry.name);
            writeln!(out, "  {label:<22} {:>16}", commas(entry.hits))?;
        }
        Ok(())
    }

    fn write_correction<W: std::io::Write>(&self, out: &mut W) -> std::io::Result<()> {
        let Some(correction) = &self.correction else {
            return Ok(());
        };
        writeln!(out, "\nCorrection")?;
        writeln!(
            out,
            "  pairs with overlap:     {:>16}",
            commas(correction.pairs_with_overlap)
        )?;
        writeln!(
            out,
            "  overlap mismatches:     {:>16}",
            commas(correction.overlap_mismatches)
        )?;
        writeln!(
            out,
            "  corrected pairs:        {:>16}",
            commas(correction.corrected_pairs)
        )?;
        writeln!(
            out,
            "  corrected bases R1/R2:  {:>16} / {}",
            commas(correction.corrected_r1_bases),
            commas(correction.corrected_r2_bases)
        )?;
        writeln!(
            out,
            "  unresolved mismatches:  {:>16}",
            commas(correction.unresolved_mismatches)
        )?;
        if correction.noncanonical_donors_skipped > 0 {
            writeln!(
                out,
                "  ambiguous donors:       {:>16}",
                commas(correction.noncanonical_donors_skipped)
            )?;
        }
        Ok(())
    }

    fn write_totals<W: std::io::Write>(&self, out: &mut W) -> std::io::Result<()> {
        writeln!(out, "\nOutput")?;
        writeln!(
            out,
            "  accepted records:       {:>16}",
            commas(self.counts.records_out)
        )?;
        writeln!(
            out,
            "  accepted bases:         {:>16}",
            commas(self.counts.r1_bases_out + self.counts.r2_bases_out)
        )?;
        let orphan_bases = self.counts.r1_bases_orphaned + self.counts.r2_bases_orphaned;
        if orphan_bases > 0 {
            writeln!(
                out,
                "  orphan bases:           {:>16}",
                commas(orphan_bases)
            )?;
        }

        writeln!(out, "\nPerformance")?;
        writeln!(
            out,
            "  elapsed:                {:>16.3} s",
            self.performance.elapsed_seconds
        )?;
        writeln!(
            out,
            "  records/second:         {:>16}",
            commas(self.performance.records_per_second as u64)
        )?;
        writeln!(
            out,
            "  bases/second:           {:>16}",
            commas(self.performance.bases_per_second as u64)
        )?;
        Ok(())
    }
}

fn adapter_report(
    stage: &crate::adapter::AdapterStage,
    stats: &Stats,
    detection: Option<&crate::sniff::adapters::AdapterSniff>,
) -> AdapterReport {
    let mut per_adapter = Vec::new();
    for (mate, adapters) in [(Mate::R1, &stage.r1), (Mate::R2, &stage.r2)] {
        for (index, adapter) in adapters.iter().enumerate() {
            per_adapter.push(AdapterEntry {
                name: adapter.name.clone(),
                mate: mate.name(),
                hits: stats.adapter_hit_count(mate, index),
            });
        }
    }
    AdapterReport {
        r1_reads_trimmed: stats.adapter_reads[0],
        r2_reads_trimmed: stats.adapter_reads[1],
        r1_bases_removed: stats.adapter_bases[0],
        r2_bases_removed: stats.adapter_bases[1],
        r1_overlap_reads_trimmed: stats.overlap_reads[0],
        r2_overlap_reads_trimmed: stats.overlap_reads[1],
        r1_overlap_bases_removed: stats.overlap_bases[0],
        r2_overlap_bases_removed: stats.overlap_bases[1],
        per_adapter,
        detection: detection.cloned(),
    }
}

fn filter_report(stats: &Stats) -> FilterReport {
    FilterReport {
        accepted_records: stats.records_out,
        rejected_records: stats.records_rejected,
        orphan_r1_records: stats.records_orphaned[0],
        orphan_r2_records: stats.records_orphaned[1],
        r1_only_failed: stats.r1_only_failed,
        r2_only_failed: stats.r2_only_failed,
        both_failed: stats.both_failed,
        per_reason: REASONS
            .iter()
            .enumerate()
            .map(|(bit, (_, name))| ReasonEntry {
                reasons: (*name).to_string(),
                records: stats.reason_counts[bit],
            })
            .collect(),
        per_combination: stats
            .reason_combinations
            .iter()
            .map(|(bits, count)| ReasonEntry {
                reasons: FilterReason::from_bits_truncate(*bits).label(),
                records: *count,
            })
            .collect(),
    }
}

/// Builds the correction section from accumulated statistics.
fn correction_report(stats: &Stats) -> CorrectionReport {
    const BASES: [char; 5] = ['A', 'C', 'G', 'T', 'N'];
    CorrectionReport {
        pairs_examined: stats.correction_examined,
        pairs_with_overlap: stats.correction_overlaps,
        pairs_with_mismatches: stats.correction_mismatched_pairs,
        overlap_mismatches: stats.correction_mismatches,
        corrected_pairs: stats.corrected_pairs,
        corrected_pairs_both_mates: stats.corrected_both_mates,
        corrected_reads: stats.corrected_reads[0] + stats.corrected_reads[1],
        corrected_r1_reads: stats.corrected_reads[0],
        corrected_r2_reads: stats.corrected_reads[1],
        corrected_bases: stats.corrected_bases[0] + stats.corrected_bases[1],
        corrected_r1_bases: stats.corrected_bases[0],
        corrected_r2_bases: stats.corrected_bases[1],
        unresolved_mismatches: stats.correction_unresolved,
        noncanonical_donors_skipped: stats.correction_noncanonical,
        corrections_per_pair: stats
            .correction_histogram
            .iter()
            .map(|(corrections, pairs)| CountEntry {
                corrections: *corrections,
                pairs: *pairs,
            })
            .collect(),
        substitutions: stats
            .correction_substitutions
            .iter()
            .enumerate()
            .flat_map(|(original, row)| {
                row.iter()
                    .enumerate()
                    .filter_map(move |(corrected, bases)| {
                        (*bases > 0).then_some(SubstitutionEntry {
                            original: BASES[original],
                            corrected: BASES[corrected],
                            bases: *bases,
                        })
                    })
            })
            .collect(),
        corrected_pairs_by_disposition: stats
            .corrected_by_disposition
            .iter()
            .map(|(disposition, pairs)| DispositionEntry {
                disposition: (*disposition).to_string(),
                pairs: *pairs,
            })
            .collect(),
    }
}

/// Builds the segmentation section from accumulated statistics.
fn segment_report(stage: &crate::segment::SegmentStage, stats: &Stats) -> SegmentReport {
    SegmentReport {
        source_records: stats.segment_sources,
        source_records_split: stats.segment_split_sources,
        delimiters_accepted: stats.segment_boundaries,
        candidates_suppressed: stats.segment_suppressed,
        fragments_emitted: stats.segment_internal + stats.segment_terminal,
        internal_fragments: stats.segment_internal,
        terminal_fragments: stats.segment_terminal,
        discarded_empty: stats.segment_empty,
        discarded_too_short: stats.segment_short,
        discarded_terminal: stats.segment_terminal_discarded,
        discarded_over_limit: stats.segment_over_limit,
        // For a segmenting run the output records are the fragments, so these are
        // the same counters the totals use, named for what they hold here.
        fragments_accepted: stats.records_out,
        fragments_rejected: stats.records_rejected,
        max_fragments_per_source: stats.segment_max_fragments,
        fragment_length: FragmentLengthSummary {
            shortest: stats.segment_length_min,
            longest: stats.segment_length_max,
            mean: {
                let emitted = stats.segment_internal + stats.segment_terminal;
                if emitted == 0 {
                    0.0
                } else {
                    stats.segment_length_total as f64 / emitted as f64
                }
            },
        },
        fragments_per_source: stats
            .segment_histogram
            .iter()
            .map(|(fragments, source_records)| FragmentCountEntry {
                fragments: *fragments,
                source_records: *source_records,
            })
            .collect(),
        per_adapter: stage
            .adapters
            .iter()
            .enumerate()
            .map(|(index, adapter)| AdapterEntry {
                name: adapter.name.clone(),
                mate: Mate::R1.name(),
                hits: stats.adapter_hits[0].get(index).copied().unwrap_or(0),
            })
            .collect(),
        terminal_policy: stage.terminal,
        min_segment_length: stage.min_segment_length as u64,
        max_segments_per_read: stage.max_segments as u64,
    }
}

/// Builds the linked-segmentation section from accumulated statistics.
fn linked_report(stage: &crate::linked::LinkedStage, stats: &Stats) -> LinkedReport {
    LinkedReport {
        both_adapters: stats.linked_both[0] + stats.linked_both[1],
        five_prime_only: stats.linked_five_only[0] + stats.linked_five_only[1],
        three_prime_only: stats.linked_three_only[0] + stats.linked_three_only[1],
        invalid_order: stats.linked_out_of_order[0] + stats.linked_out_of_order[1],
        insert_too_short: stats.linked_short_insert[0] + stats.linked_short_insert[1],
        neither_adapter: stats.linked_neither[0] + stats.linked_neither[1],
        leading_bases_removed: stats.linked_leading_bases[0] + stats.linked_leading_bases[1],
        trailing_bases_removed: stats.linked_trailing_bases[0] + stats.linked_trailing_bases[1],
        per_definition: [(Mate::R1, &stage.r1), (Mate::R2, &stage.r2)]
            .into_iter()
            .flat_map(|(mate, definitions)| {
                definitions
                    .iter()
                    .enumerate()
                    .map(move |(index, definition)| LinkedEntry {
                        name: definition.name.clone(),
                        mate: mate.name(),
                        matches: stats.linked_hits[mate.index()]
                            .get(index)
                            .copied()
                            .unwrap_or(0),
                        both_adapters_required: definition.require == crate::linked::Require::Both,
                    })
            })
            .collect(),
        unmatched_policy: stage.unmatched,
    }
}

fn trim_entries(stage: &TrimStage, stats: &Stats) -> Vec<TrimEntry> {
    TrimOp::ALL
        .iter()
        .filter(|op| operation_enabled(stage, **op))
        .map(|op| TrimEntry {
            operation: op.name(),
            reads_trimmed: stats.trim_reads[*op as usize],
            bases_removed: stats.trim_bases[*op as usize],
        })
        .collect()
}

fn operation_enabled(stage: &TrimStage, op: TrimOp) -> bool {
    [&stage.r1, &stage.r2].iter().any(|mate| match op {
        TrimOp::Fixed => mate.front > 0 || mate.tail > 0,
        TrimOp::Quality => mate.needs_quality(),
        TrimOp::TerminalN => mate.terminal_n,
        TrimOp::PolyG => mate.poly_g.is_some(),
        TrimOp::PolyX => mate.poly_x.is_some(),
        TrimOp::MaxLength => mate.max_length.is_some(),
    })
}

/// Flattens a JSON value into dotted key/value rows.
fn flatten(value: &serde_json::Value, prefix: String, rows: &mut Vec<(String, String)>) {
    match value {
        serde_json::Value::Object(map) => {
            for (key, value) in map {
                let key = if prefix.is_empty() {
                    key.clone()
                } else {
                    format!("{prefix}.{key}")
                };
                flatten(value, key, rows);
            }
        }
        serde_json::Value::Array(items) => {
            for (index, value) in items.iter().enumerate() {
                flatten(value, format!("{prefix}[{index}]"), rows);
            }
        }
        serde_json::Value::String(text) => rows.push((prefix, text.clone())),
        serde_json::Value::Null => rows.push((prefix, String::new())),
        other => rows.push((prefix, other.to_string())),
    }
}

/// Formats an integer with thousands separators.
#[must_use]
pub fn commas(value: u64) -> String {
    let digits = value.to_string();
    let mut out = String::with_capacity(digits.len() + digits.len() / 3);
    for (index, digit) in digits.chars().enumerate() {
        if index > 0 && (digits.len() - index).is_multiple_of(3) {
            out.push(',');
        }
        out.push(digit);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn binseq_version_matches_manifest() {
        let manifest = include_str!("../Cargo.toml");
        assert!(
            manifest.contains(&format!("binseq = {{ version = \"{BINSEQ_VERSION}\"")),
            "BINSEQ_VERSION must match the binseq dependency in Cargo.toml"
        );
    }

    #[test]
    fn thousands_separators() {
        assert_eq!(commas(0), "0");
        assert_eq!(commas(999), "999");
        assert_eq!(commas(1_000), "1,000");
        assert_eq!(commas(25_000_000), "25,000,000");
        assert_eq!(commas(1_234_567_890), "1,234,567,890");
    }

    #[test]
    fn json_flattening_produces_dotted_keys() {
        let value = serde_json::json!({
            "a": {"b": 1, "c": "x"},
            "d": [{"e": true}, {"e": false}],
            "f": null,
        });
        let mut rows = Vec::new();
        flatten(&value, String::new(), &mut rows);
        assert!(rows.contains(&("a.b".to_string(), "1".to_string())));
        assert!(rows.contains(&("a.c".to_string(), "x".to_string())));
        assert!(rows.contains(&("d[0].e".to_string(), "true".to_string())));
        assert!(rows.contains(&("d[1].e".to_string(), "false".to_string())));
        assert!(rows.contains(&("f".to_string(), String::new())));
    }
}
