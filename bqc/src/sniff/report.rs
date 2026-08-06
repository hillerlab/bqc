// Copyright (c) 2026 The Hiller Lab at the Senckenberg Gesellschaft für Naturforschung
// Distributed under the terms of the GNU General Public License, Version 3.0.

//! Report projection for both `sniff` subcommands.
//!
//! One common metadata structure and one report type per subcommand — no
//! generic report tree, no inheritance. JSON is the stable pipeline interface
//! and carries a `schema_version`; text is for a human reading a terminal; TSV
//! is for aggregating a cohort into one table.
//!
//! Every algorithmic constant that shaped a result is serialized with it, so a
//! report explains itself without the reader having to know this version's
//! defaults.

use std::fmt::Write as _;

use serde::Serialize;

use crate::error::Result;
use crate::io::CbqInput;
use crate::report::commas;
use crate::sniff::sample::SamplePlan;
use crate::sniff::{Confidence, Format, SCHEMA_VERSION, adapters};

/// Serializes a report as pretty JSON with a trailing newline.
fn render_json<T: Serialize>(value: &T) -> Result<String> {
    serde_json::to_string_pretty(value)
        .map(|json| json + "\n")
        .map_err(|e| crate::error::Error::config(format!("failed to serialize report: {e}")))
}

/// The tool that produced a report.
#[derive(Debug, Clone, Serialize)]
pub struct Tool {
    pub name: &'static str,
    pub version: &'static str,
}

impl Default for Tool {
    fn default() -> Self {
        Self {
            name: "bqc",
            version: env!("CARGO_PKG_VERSION"),
        }
    }
}

/// What was inspected.
#[derive(Debug, Clone, Serialize)]
pub struct Input {
    pub path: String,
    pub records: u64,
    pub paired: bool,
    pub quality: bool,
    pub headers: bool,
}

impl Input {
    #[must_use]
    pub fn new(input: &CbqInput) -> Self {
        let schema = input.schema();
        Self {
            path: input.path().display().to_string(),
            records: input.num_records(),
            paired: schema.paired,
            quality: schema.quality,
            headers: schema.headers,
        }
    }
}

/// Which records were inspected, and how they were chosen.
#[derive(Debug, Clone, Serialize)]
pub struct Sample {
    /// Always `deterministic-distributed`; recorded so a consumer can tell this
    /// apart from a leading-prefix sample without inspecting the version.
    pub method: &'static str,
    pub range_start: u64,
    pub range_end: u64,
    pub requested: u64,
    pub selected: u64,
}

impl Sample {
    #[must_use]
    pub fn new(plan: &SamplePlan) -> Self {
        Self {
            method: "deterministic-distributed",
            range_start: plan.range_start,
            range_end: plan.range_end,
            requested: plan.requested,
            selected: plan.selected,
        }
    }
}

/// Metadata shared by every `sniff` report.
#[derive(Debug, Clone, Serialize)]
pub struct Meta {
    pub schema_version: u32,
    pub command: &'static str,
    pub tool: Tool,
    pub input: Input,
    pub sample: Sample,
}

impl Meta {
    #[must_use]
    pub fn new(command: &'static str, input: &CbqInput, plan: &SamplePlan) -> Self {
        Self {
            schema_version: SCHEMA_VERSION,
            command,
            tool: Tool::default(),
            input: Input::new(input),
            sample: Sample::new(plan),
        }
    }
}

/// The `sniff adapters` report.
#[derive(Debug, Clone, Serialize)]
pub struct AdapterReport {
    #[serde(flatten)]
    pub meta: Meta,
    pub parameters: adapters::Params,
    pub result: adapters::AdapterSniff,
    pub warnings: Vec<String>,
}

impl AdapterReport {
    #[must_use]
    pub fn new(
        input: &CbqInput,
        plan: &SamplePlan,
        parameters: adapters::Params,
        result: adapters::AdapterSniff,
    ) -> Self {
        let warnings = adapter_warnings(&result, plan);
        Self {
            meta: Meta::new("sniff adapters", input, plan),
            parameters,
            result,
            warnings,
        }
    }

    /// Renders the report in the requested projection.
    pub fn render(&self, format: Format) -> Result<String> {
        match format {
            Format::Json => render_json(self),
            Format::Tsv => Ok(self.tsv()),
            Format::Text => Ok(self.text()),
        }
    }

    /// One row per candidate, for cohort aggregation.
    fn tsv(&self) -> String {
        let mut out = String::from(
            "input\tmate\tdecision\tsequence\tknown_name\tknown_category\tconfidence\tevidence_sources\t\
             supporting_reads\tsupport_fraction\ttail_connected_fraction\ttail_enrichment\t\
             median_start\tmedian_distance_to_end\texact_matches\tsubstitution_matches\t\
             indel_matches\tmean_error_rate\tpaired_overlap_support\n",
        );
        let mates = std::iter::once(&self.result.r1).chain(self.result.r2.as_ref());
        for mate in mates {
            let path = &self.meta.input.path;
            let decision = mate.decision.name();
            if mate.candidates.is_empty() {
                let _ = writeln!(
                    out,
                    "{path}\t{}\t{decision}\t.\t.\t.\t.\t.\t0\t0\t0\t.\t0\t0\t0\t0\t0\t0\t0",
                    mate.mate
                );
            }
            for candidate in &mate.candidates {
                let _ = writeln!(
                    out,
                    "{path}\t{}\t{decision}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.6}\t{}",
                    mate.mate,
                    String::from_utf8_lossy(&candidate.sequence),
                    candidate.known_name.as_deref().unwrap_or("."),
                    candidate.known_category.as_deref().unwrap_or("."),
                    candidate.confidence.name(),
                    sources(candidate, ","),
                    candidate.supporting_reads,
                    candidate.support_fraction,
                    candidate.tail_connected_fraction,
                    tsv_enrichment(candidate.tail_enrichment),
                    candidate.median_start,
                    candidate.median_distance_to_end,
                    candidate.exact_matches,
                    candidate.substitution_matches,
                    candidate.indel_matches,
                    candidate.mean_error_rate,
                    candidate.paired_overlap_support,
                );
            }
        }
        out
    }

    fn text(&self) -> String {
        let mut out = String::new();
        let _ = writeln!(out, "Input: {}", self.meta.input.path);
        let _ = writeln!(
            out,
            "Records available: {}",
            commas(self.meta.input.records)
        );
        let _ = writeln!(
            out,
            "Records sampled: {}",
            commas(self.meta.sample.selected)
        );
        let _ = writeln!(out, "Sampling: {}", self.meta.sample.method);

        let mates = std::iter::once(&self.result.r1).chain(self.result.r2.as_ref());
        for mate in mates {
            let _ = writeln!(
                out,
                "\n{} adapter result: {}",
                mate.mate,
                mate.decision.name()
            );
            if mate.candidates.is_empty() {
                out.push_str("  No candidate found.\n");
            }
            for candidate in &mate.candidates {
                let _ = writeln!(
                    out,
                    "  Sequence: {}",
                    String::from_utf8_lossy(&candidate.sequence)
                );
                if let Some(name) = &candidate.known_name {
                    let _ = writeln!(out, "  Known as: {name}");
                }
                if let Some(category) = &candidate.known_category {
                    let _ = writeln!(out, "  Category: {category}");
                }
                let _ = writeln!(out, "  Confidence: {}", candidate.confidence.name());
                let _ = writeln!(
                    out,
                    "  Supporting reads: {} ({:.2}%)",
                    commas(candidate.supporting_reads),
                    candidate.support_fraction * 100.0
                );
                let _ = writeln!(
                    out,
                    "  3' connected matches: {:.1}%",
                    candidate.tail_connected_fraction * 100.0
                );
                let _ = writeln!(
                    out,
                    "  Tail enrichment: {}",
                    enrichment(candidate.tail_enrichment)
                );
                let _ = writeln!(out, "  Evidence: {}", sources(candidate, ", "));
                if candidate.confidence != Confidence::High {
                    out.push('\n');
                }
            }
            if mate.poly_a_signal > 0 || mate.poly_g_signal > 0 {
                let _ = writeln!(
                    out,
                    "  Artifact signal: {} poly-A, {} poly-G reads",
                    commas(mate.poly_a_signal),
                    commas(mate.poly_g_signal)
                );
            }
        }
        for warning in &self.warnings {
            let _ = writeln!(out, "\nWarning: {warning}");
        }
        out
    }
}

/// A candidate's evidence sources, joined for display.
fn sources(candidate: &adapters::Candidate, separator: &str) -> String {
    candidate
        .evidence_sources
        .iter()
        .map(|source| source.name())
        .collect::<Vec<_>>()
        .join(separator)
}

fn enrichment(value: Option<f64>) -> String {
    value.map_or_else(|| "tail only".to_string(), |value| format!("{value:.1}x"))
}

fn tsv_enrichment(value: Option<f64>) -> String {
    value.map_or_else(|| "tail_only".to_string(), |value| format!("{value:.4}"))
}

/// Conditions worth surfacing that are not failures.
fn adapter_warnings(result: &adapters::AdapterSniff, plan: &SamplePlan) -> Vec<String> {
    let mut warnings = Vec::new();
    if plan.selected < plan.requested {
        warnings.push(format!(
            "requested {} records but the input holds {}; sampled all of them",
            plan.requested,
            plan.available()
        ));
    }
    let mates = std::iter::once(&result.r1).chain(result.r2.as_ref());
    for mate in mates {
        if mate.decision == crate::sniff::Decision::Mixed {
            warnings.push(format!(
                "{} has more than one unrelated high-confidence adapter; \
                 this usually means pooled libraries or concatenated runs, \
                 and no sequence is recommended automatically",
                mate.mate
            ));
        }
        if mate.poly_a_signal * 2 > mate.sampled_reads && mate.sampled_reads > 0 {
            warnings.push(format!(
                "{}: over half of sampled reads end in a poly-A run",
                mate.mate
            ));
        }
        if mate.poly_g_signal * 2 > mate.sampled_reads && mate.sampled_reads > 0 {
            warnings.push(format!(
                "{}: over half of sampled reads end in a poly-G run, \
                 which on a two-colour instrument indicates no signal",
                mate.mate
            ));
        }
    }
    warnings
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sniff::Decision;

    fn empty_mate(name: &'static str) -> adapters::MateResult {
        adapters::MateResult {
            mate: name,
            decision: Decision::Inconclusive,
            recommended_sequence: None,
            recommended_name: None,
            candidates: Vec::new(),
            sampled_reads: 100,
            informative_reads: 0,
            poly_a_signal: 0,
            poly_g_signal: 0,
        }
    }

    #[test]
    fn a_short_input_is_reported_as_a_warning_not_a_failure() {
        let plan = SamplePlan::new(0, 10, 1000);
        let result = adapters::AdapterSniff {
            database_version: 1,
            r1: empty_mate("R1"),
            r2: None,
        };
        let warnings = adapter_warnings(&result, &plan);
        assert_eq!(warnings.len(), 1);
        assert!(warnings[0].contains("sampled all of them"));
    }

    #[test]
    fn tsv_writes_a_row_even_when_nothing_was_found() {
        let plan = SamplePlan::new(0, 1000, 100);
        let result = adapters::AdapterSniff {
            database_version: 1,
            r1: empty_mate("R1"),
            r2: Some(empty_mate("R2")),
        };
        let report = AdapterReport {
            meta: Meta {
                schema_version: SCHEMA_VERSION,
                command: "sniff adapters",
                tool: Tool::default(),
                input: Input {
                    path: "x.cbq".to_string(),
                    records: 1000,
                    paired: true,
                    quality: true,
                    headers: true,
                },
                sample: Sample::new(&plan),
            },
            parameters: adapters::Params::default(),
            result,
            warnings: Vec::new(),
        };
        let tsv = report.render(Format::Tsv).unwrap();
        let lines: Vec<&str> = tsv.lines().collect();
        assert_eq!(lines.len(), 3, "header plus one row per mate");
        assert!(lines[1].starts_with("x.cbq\tR1\tinconclusive\t"));
        assert!(lines[2].starts_with("x.cbq\tR2\tinconclusive\t"));
        let columns = lines[0].split('\t').count();
        assert!(
            lines.iter().all(|line| line.split('\t').count() == columns),
            "every TSV row follows the header schema: {tsv}"
        );
    }

    #[test]
    fn tail_only_enrichment_has_explicit_projection_values() {
        assert_eq!(enrichment(None), "tail only");
        assert_eq!(tsv_enrichment(None), "tail_only");
        assert_eq!(enrichment(Some(5.25)), "5.2x");
        assert_eq!(tsv_enrichment(Some(5.25)), "5.2500");
    }

    #[test]
    fn json_carries_the_schema_version_and_the_gates() {
        let plan = SamplePlan::new(0, 1000, 100);
        let report = AdapterReport {
            meta: Meta {
                schema_version: SCHEMA_VERSION,
                command: "sniff adapters",
                tool: Tool::default(),
                input: Input {
                    path: "x.cbq".to_string(),
                    records: 1000,
                    paired: false,
                    quality: true,
                    headers: true,
                },
                sample: Sample::new(&plan),
            },
            parameters: adapters::Params::default(),
            result: adapters::AdapterSniff {
                database_version: 1,
                r1: empty_mate("R1"),
                r2: None,
            },
            warnings: Vec::new(),
        };
        let json: serde_json::Value =
            serde_json::from_str(&report.render(Format::Json).unwrap()).unwrap();
        assert_eq!(json["schema_version"], 1);
        assert_eq!(json["command"], "sniff adapters");
        assert_eq!(json["sample"]["method"], "deterministic-distributed");
        // The constants that shaped the result travel with it.
        assert_eq!(json["parameters"]["gates"]["min_supporting_reads"], 20);
        assert_eq!(json["result"]["r1"]["decision"], "inconclusive");
    }
}

/// The `sniff strand` report.
#[cfg(feature = "sniff-strand")]
#[derive(Debug, Clone, Serialize)]
pub struct StrandReport {
    #[serde(flatten)]
    pub meta: Meta,
    pub parameters: crate::sniff::strand::Params,
    pub result: crate::sniff::strand::StrandSniff,
    pub warnings: Vec<String>,
}

#[cfg(feature = "sniff-strand")]
impl StrandReport {
    #[must_use]
    pub fn new(
        input: &CbqInput,
        plan: &SamplePlan,
        parameters: crate::sniff::strand::Params,
        result: crate::sniff::strand::StrandSniff,
    ) -> Self {
        let warnings = strand_warnings(&result);
        Self {
            meta: Meta::new("sniff strand", input, plan),
            parameters,
            result,
            warnings,
        }
    }

    /// Renders the report in the requested projection.
    pub fn render(&self, format: Format) -> Result<String> {
        match format {
            Format::Json => render_json(self),
            Format::Tsv => Ok(self.tsv()),
            Format::Text => Ok(self.text()),
        }
    }

    /// One summary row, for aggregating a cohort into one table.
    fn tsv(&self) -> String {
        let result = &self.result;
        let mut out = String::from(
            "input\tdecision\tsalmon_library_type\tstrandedness\tpair_orientation\t\
             forward_fraction\treverse_fraction\tinformative_records\tinformative_fraction\t\
             records_examined\tfeaturecounts_strand\thtseq_stranded\n",
        );
        let _ = writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{}\t{:.6}\t{}\t{}\t{}",
            self.meta.input.path,
            result.decision,
            result.salmon_library_type.unwrap_or("."),
            result.strandedness.name(),
            result.pair_orientation.name(),
            result.forward_fraction,
            result.reverse_fraction,
            result.informative_records,
            result.informative_fraction,
            result.records_examined,
            result
                .featurecounts_strand
                .map_or_else(|| ".".to_string(), |value| value.to_string()),
            result.htseq_stranded.unwrap_or("."),
        );
        out
    }

    fn text(&self) -> String {
        let result = &self.result;
        let mut out = String::new();
        let _ = writeln!(out, "Input: {}", self.meta.input.path);
        let _ = writeln!(out, "Index: {}", result.index_metadata.path);
        let _ = writeln!(
            out,
            "Records available: {}",
            commas(self.meta.input.records)
        );
        let _ = writeln!(out, "Records examined: {}", commas(result.records_examined));
        let _ = writeln!(out, "Sampling: {}", self.meta.sample.method);
        let _ = writeln!(out);
        let _ = writeln!(
            out,
            "Library type: {}",
            result.salmon_library_type.unwrap_or("undetermined")
        );
        let _ = writeln!(out, "Strandedness: {}", result.strandedness.name());
        let _ = writeln!(out, "Pair orientation: {}", result.pair_orientation.name());
        let _ = writeln!(
            out,
            "Informative fragments: {} ({:.1}%)",
            commas(result.informative_records),
            result.informative_fraction * 100.0
        );
        let _ = writeln!(
            out,
            "Forward fraction: {:.1}%",
            result.forward_fraction * 100.0
        );
        let _ = writeln!(
            out,
            "Reverse fraction: {:.1}%",
            result.reverse_fraction * 100.0
        );
        let _ = writeln!(out, "Decision: {}", result.decision);
        if let Some(reason) = &result.failure_reason {
            let _ = writeln!(out, "Reason: {reason}");
        }
        for warning in &self.warnings {
            let _ = writeln!(out, "\nWarning: {warning}");
        }
        out
    }
}

/// Conditions worth surfacing that are not failures.
#[cfg(feature = "sniff-strand")]
fn strand_warnings(result: &crate::sniff::strand::StrandSniff) -> Vec<String> {
    use crate::sniff::strand::Orientation;
    let mut warnings = Vec::new();
    // An unusual orientation is a real finding, not noise to be smoothed away.
    match result.pair_orientation {
        Orientation::Outward => warnings.push(
            "the dominant pair orientation is outward, not inward; \
             check the library preparation before using these reads"
                .to_string(),
        ),
        Orientation::Matching => warnings.push(
            "the dominant pair orientation is matching, not inward; \
             check the library preparation before using these reads"
                .to_string(),
        ),
        Orientation::Inward | Orientation::Undetermined => {}
    }
    if result.decoy_only_records * 4 > result.records_examined && result.records_examined > 0 {
        warnings.push(
            "over a quarter of examined records mapped only to decoys, \
             which usually means genomic contamination"
                .to_string(),
        );
    }
    if result.mapping_errors > 0 {
        warnings.push(format!(
            "{} records were shorter than the index k-mer and could not be mapped",
            result.mapping_errors
        ));
    }
    warnings
}
