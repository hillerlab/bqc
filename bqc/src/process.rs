// Copyright (c) 2026 The Hiller Lab at the Senckenberg Gesellschaft für Naturforschung
// Distributed under the terms of the GNU General Public License, Version 3.0.

//! The fused record processor.
//!
//! [`Workflow`] is the single implementation of the `adapter -> trim -> filter`
//! pipeline. The `adapter`, `trim` and `filter` subcommands construct a
//! workflow with one enabled stage, so a standalone command is algorithmically
//! identical to `bqc workflow --steps <stage>`.

use serde::{Deserialize, Serialize};

use crate::adapter::{AdapterHit, AdapterStage};
use crate::error::{Error, Result};
use crate::filter::{FilterReason, FilterStage};
use crate::read::{validate_quality, ReadView, Span};
use crate::trim::{TrimOutcome, TrimStage};

/// Which mate of a record is being processed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Mate {
    R1,
    R2,
}

impl Mate {
    /// Stable name used in reports and sidecars.
    #[must_use]
    pub fn name(self) -> &'static str {
        match self {
            Self::R1 => "R1",
            Self::R2 => "R2",
        }
    }

    /// Index into per-mate statistics arrays.
    #[must_use]
    pub fn index(self) -> usize {
        self as usize
    }
}

/// How rejected records are written to the failed output.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, clap::ValueEnum)]
#[serde(rename_all = "lowercase")]
pub enum FailedMode {
    /// Write the original, untransformed record (default).
    Original,
    /// Write the transformed record that failed filtering.
    Processed,
}

/// Pair retention policy.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, clap::ValueEnum)]
#[serde(rename_all = "lowercase")]
pub enum PairPolicy {
    /// A pair is accepted only when both mates pass.
    Strict,
    /// Pairs with one failing mate contribute the surviving mate to a
    /// single-end orphan output.
    Orphan,
}

/// Where a processed record is routed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PairDisposition {
    /// Both mates passed: the pair enters the accepted output.
    Accepted,
    /// Neither mate passed (or strict pairing with any failure): the pair
    /// enters the rejected output.
    Rejected,
    /// R1 passed and R2 failed (orphan policy): R1 becomes a single-end orphan.
    OrphanR1,
    /// R2 passed and R1 failed (orphan policy): R2 becomes a single-end orphan.
    OrphanR2,
}

impl PairDisposition {
    /// Stable name of the routing decision, shared by the correction log and the
    /// report so the two can never disagree about where a pair went.
    #[must_use]
    pub fn name(self) -> &'static str {
        match self {
            Self::Accepted => "accepted",
            Self::Rejected => "rejected",
            Self::OrphanR1 => "orphan_r1",
            Self::OrphanR2 => "orphan_r2",
        }
    }
}

/// The outcome of processing one fragment of a segmented read.
#[derive(Debug, Clone, Copy)]
pub struct FragmentResult {
    /// Where the fragment came from, and what delimited it.
    pub fragment: crate::segment::Fragment,
    /// Coordinates retained after trimming the fragment, in source-read
    /// coordinates.
    pub retained: Span,
    pub trimming: TrimOutcome,
    pub reasons: FilterReason,
}

impl FragmentResult {
    /// Whether this fragment passed every configured filter.
    #[must_use]
    pub fn passed(&self) -> bool {
        self.reasons.passed()
    }
}

/// The outcome of processing one mate.
#[derive(Debug, Clone, Copy)]
pub struct MateResult {
    /// Coordinates retained from the original read.
    pub retained: Span,
    /// Length of the original read.
    pub original_length: usize,
    /// Length after adapter removal, before trimming.
    pub adapter_trimmed_length: usize,
    pub adapter_hit: Option<AdapterHit>,
    /// What linked segmentation did to this mate, when the stage ran.
    pub linked: Option<crate::linked::LinkedOutcome>,
    /// Bases removed by paired-overlap inference.
    pub overlap_removed: usize,
    pub trimming: TrimOutcome,
    pub reasons: FilterReason,
}

impl MateResult {
    /// Whether this mate passed every configured filter.
    #[must_use]
    pub fn passed(&self) -> bool {
        self.reasons.passed()
    }

    /// Bases removed by adapter trimming.
    #[must_use]
    pub fn adapter_bases_removed(&self) -> usize {
        self.original_length - self.adapter_trimmed_length
    }

    /// Final retained length.
    #[must_use]
    pub fn final_length(&self) -> usize {
        self.retained.len()
    }
}

/// The outcome of processing one record (a single read or a pair).
#[derive(Debug, Clone, Copy)]
pub struct PairResult {
    pub r1: MateResult,
    pub r2: Option<MateResult>,
    /// Where the record is routed.
    pub disposition: PairDisposition,
    /// The overlap alignment, when one was computed and accepted. Shared by
    /// adapter trimming and correction: one pair, one analysis.
    pub overlap: Option<crate::overlap::Overlap>,
    /// What base correction did to this pair.
    pub correction: crate::correct::CorrectionSummary,
}

impl PairResult {
    /// Iterates over the processed mates.
    pub fn mates(&self) -> impl Iterator<Item = (Mate, &MateResult)> {
        std::iter::once((Mate::R1, &self.r1)).chain(self.r2.iter().map(|result| (Mate::R2, result)))
    }

    /// Whether the record enters the accepted output.
    #[must_use]
    pub fn accepted(&self) -> bool {
        self.disposition == PairDisposition::Accepted
    }
}

/// A compiled, ready-to-run processing pipeline.
#[derive(Debug, Clone, Serialize)]
pub struct Workflow {
    pub adapter: Option<AdapterStage>,
    /// Linked segmentation, which runs at the start of the adapter stage.
    pub linked: Option<crate::linked::LinkedStage>,
    pub trim: Option<TrimStage>,
    pub filter: Option<FilterStage>,
    pub correction: Option<crate::correct::CorrectionStage>,
    /// Internal adapter splitting, which replaces the adapter stage: adapters
    /// become delimiters instead of ends to trim from. Only the `segment`
    /// command sets this, and it is mutually exclusive with everything above.
    pub segment: Option<crate::segment::SegmentStage>,
    /// Overlap parameters for the single per-pair analysis. Present when adapter
    /// overlap inference or correction needs it; both read this one copy, so the
    /// thresholds are never duplicated.
    pub overlap: Option<crate::overlap::OverlapParams>,
    pub pair_policy: PairPolicy,
}

impl Workflow {
    /// Builds a workflow, rejecting configurations with no effective operation.
    pub fn new(
        adapter: Option<AdapterStage>,
        linked: Option<crate::linked::LinkedStage>,
        trim: Option<TrimStage>,
        filter: Option<FilterStage>,
        correction: Option<crate::correct::CorrectionStage>,
    ) -> Result<Self> {
        if adapter.is_none()
            && linked.is_none()
            && trim.is_none()
            && filter.is_none()
            && correction.is_none()
        {
            return Err(Error::config(
                "no operation configured; enable at least one adapter, trim, filter \
                 or correction option",
            ));
        }
        // One resolved copy of the overlap parameters, used by whichever stages
        // need it. Correction needs the analysis even when adapter trimming does
        // not act on it.
        let overlap = adapter
            .as_ref()
            .and_then(|stage| stage.paired_overlap)
            .or_else(|| correction.map(|_| crate::overlap::OverlapParams::default()));
        Ok(Self {
            adapter,
            linked,
            trim,
            filter,
            correction,
            segment: None,
            overlap,
            pair_policy: PairPolicy::Strict,
        })
    }

    /// Builds a segmenting workflow: adapters cut reads into fragments, and the
    /// optional trim and filter stages then apply to each fragment.
    ///
    /// Segmentation has its own constructor because it is the one workflow whose
    /// output cardinality differs from its input's, so no other stage that
    /// assumes one-record-in-one-record-out may be enabled alongside it.
    #[must_use]
    pub fn segmenting(
        segment: crate::segment::SegmentStage,
        trim: Option<TrimStage>,
        filter: Option<FilterStage>,
    ) -> Self {
        Self {
            adapter: None,
            linked: None,
            trim,
            filter,
            correction: None,
            segment: Some(segment),
            overlap: None,
            pair_policy: PairPolicy::Strict,
        }
    }

    /// Overrides the resolved overlap parameters.
    ///
    /// Used when correction is enabled without adapter overlap trimming, so the
    /// shared analysis still honours the user's overlap options.
    #[must_use]
    pub fn with_overlap_params(mut self, params: crate::overlap::OverlapParams) -> Self {
        if self.overlap.is_some() {
            self.overlap = Some(params);
        }
        self
    }

    /// Sets the pair retention policy.
    ///
    /// The orphan policy only makes sense with a filter stage and paired
    /// input; both are enforced at the configuration boundary.
    #[must_use]
    pub fn with_pair_policy(mut self, policy: PairPolicy) -> Self {
        self.pair_policy = policy;
        self
    }

    /// Whether any stage reads quality values.
    #[must_use]
    pub fn needs_quality(&self) -> bool {
        self.trim.as_ref().is_some_and(TrimStage::needs_quality)
            || self.filter.as_ref().is_some_and(FilterStage::needs_quality)
            || self.correction.is_some()
    }

    /// Whether any record can be rejected by this workflow.
    #[must_use]
    pub fn can_reject(&self) -> bool {
        self.filter.is_some()
    }

    /// Names of the enabled stages, in execution order.
    #[must_use]
    pub fn stage_order(&self) -> Vec<&'static str> {
        let mut stages = Vec::new();
        if self.correction.is_some() {
            stages.push("correct");
        }
        if self.linked.is_some() {
            stages.push("linked");
        }
        if self.segment.is_some() {
            stages.push("segment");
        }
        if self.adapter.is_some() {
            stages.push("adapter");
        }
        if self.trim.is_some() {
            stages.push("trim");
        }
        if self.filter.is_some() {
            stages.push("filter");
        }
        stages
    }

    /// Processes one record.
    ///
    /// Valid for workflows without base correction, which need no scratch space;
    /// use [`Workflow::process_corrected`] when correction is enabled, because the
    /// corrected bases have to outlive this call to reach the writer.
    pub fn process(
        &self,
        record: u64,
        r1: ReadView<'_>,
        r2: Option<ReadView<'_>>,
    ) -> Result<PairResult> {
        debug_assert!(
            self.correction.is_none(),
            "a correcting workflow must use process_corrected: the corrected bases \
             live in the scratch buffer this call would discard"
        );
        self.process_corrected(
            record,
            r1,
            r2,
            &mut crate::correct::CorrectionScratch::default(),
        )
    }

    /// Processes one record, correcting bases into `scratch`.
    ///
    /// `scratch` holds the correction plan and the copy-on-write buffers for any
    /// corrected mate; it is worker-local and reused across records. After this
    /// returns, `scratch.corrected(mate)` yields the bytes the writer must emit
    /// for a corrected mate, and the original record supplies the rest.
    pub fn process_corrected(
        &self,
        record: u64,
        r1: ReadView<'_>,
        r2: Option<ReadView<'_>>,
        scratch: &mut crate::correct::CorrectionScratch,
    ) -> Result<PairResult> {
        scratch.clear();
        // One overlap analysis per pair, shared by correction and adapter
        // trimming.
        let overlap = self.overlap.and_then(|params| {
            r2.and_then(|r2| crate::overlap::find_overlap(r1.sequence, r2.sequence, params))
        });

        // Correction is planned against the original pair, then materialised
        // before any other stage, so adapter matching, trimming and filtering all
        // see corrected bases and qualities.
        let mut correction = crate::correct::CorrectionSummary::default();
        if let (Some(stage), Some(r2), Some(overlap)) = (self.correction, r2, overlap) {
            if self.needs_quality() {
                crate::read::validate_quality(
                    r1.quality.ok_or(Error::MissingQuality("base correction"))?,
                    record,
                    Mate::R1.name(),
                )?;
                crate::read::validate_quality(
                    r2.quality.ok_or(Error::MissingQuality("base correction"))?,
                    record,
                    Mate::R2.name(),
                )?;
            }
            correction = crate::correct::plan(stage, r1, r2, overlap, &mut scratch.edits);
        }
        let (r1, r2) = match r2 {
            Some(r2) => {
                let (r1, r2) = scratch.apply_pair(r1, r2);
                (r1, Some(r2))
            }
            None => (r1, None),
        };

        let insert = overlap
            .filter(|_| {
                self.adapter
                    .as_ref()
                    .is_some_and(|stage| stage.paired_overlap.is_some())
            })
            .map(|overlap| overlap.insert_length);
        let r1_result = self.process_mate(record, Mate::R1, r1, insert)?;
        let r2_result = r2
            .map(|r2| self.process_mate(record, Mate::R2, r2, insert))
            .transpose()?;
        let r1_passed = r1_result.passed();
        let r2_passed = r2_result.as_ref().is_none_or(MateResult::passed);
        let disposition = match (r1_passed, r2_passed) {
            (true, true) => PairDisposition::Accepted,
            (false, false) => PairDisposition::Rejected,
            _ if self.pair_policy == PairPolicy::Orphan && r2_result.is_some() => {
                if r1_passed {
                    PairDisposition::OrphanR1
                } else {
                    PairDisposition::OrphanR2
                }
            }
            // Strict pair retention: a pair survives only if both mates pass.
            _ => PairDisposition::Rejected,
        };
        Ok(PairResult {
            r1: r1_result,
            r2: r2_result,
            disposition,
            overlap,
            correction,
        })
    }

    /// Segments one single-end record into fragments.
    ///
    /// The fragments are cut in `scratch` and the per-fragment trim and filter
    /// outcomes are appended to `output`; both are worker-owned and reused, so
    /// the record loop allocates nothing. Fragments are appended in ascending
    /// coordinate order, which is the order they must be written in.
    pub fn process_segments(
        &self,
        record: u64,
        read: ReadView<'_>,
        scratch: &mut crate::segment::SegmentScratch,
        output: &mut Vec<FragmentResult>,
    ) -> Result<crate::segment::SegmentOutcome> {
        let Some(stage) = self.segment.as_ref() else {
            return Err(Error::config("no segmentation stage configured"));
        };
        output.clear();
        if self.needs_quality() {
            let Some(quality) = read.quality else {
                return Err(Error::MissingQuality("quality trimming or filtering"));
            };
            validate_quality(quality, record, Mate::R1.name())?;
        }
        let outcome = crate::segment::segment(stage, read.sequence, scratch);
        for fragment in &scratch.fragments {
            // Trimming and filtering see only the fragment, exactly as they would
            // see a whole read: they act on the span, and the span is the cut.
            let mut span = fragment.span;
            let trimming = match &self.trim {
                Some(stage) => crate::trim::apply(stage.mate(Mate::R1), read, &mut span),
                None => TrimOutcome::default(),
            };
            let reasons = match &self.filter {
                Some(stage) => stage.evaluate(span.sequence(read), span.quality(read)),
                None => FilterReason::empty(),
            };
            debug_assert!(span.start <= span.end && span.end <= read.len());
            output.push(FragmentResult {
                fragment: *fragment,
                retained: span,
                trimming,
                reasons,
            });
        }
        Ok(outcome)
    }

    fn process_mate(
        &self,
        record: u64,
        mate: Mate,
        read: ReadView<'_>,
        insert: Option<usize>,
    ) -> Result<MateResult> {
        if self.needs_quality() {
            let Some(quality) = read.quality else {
                return Err(Error::MissingQuality("quality trimming or filtering"));
            };
            validate_quality(quality, record, mate.name())?;
        }

        let mut span = Span::full(read.len());
        let mut overlap_removed = 0usize;

        // Linked segmentation runs first: it retains the insert between two
        // configured flanks. A successful linked match already knows both
        // boundaries, so ordinary adapter matching is skipped for that read —
        // the same precedent paired-overlap inference sets.
        let mut linked = None;
        let mut linked_matched = false;
        if let Some(stage) = self.linked.as_ref() {
            let outcome = crate::linked::find(
                stage.definitions(mate),
                self.adapter
                    .as_ref()
                    .map_or_else(crate::adapter::AdapterParams::default, |a| a.params),
                span.sequence(read),
            );
            if let crate::linked::LinkedOutcome::Matched(found) = outcome {
                // Coordinates come back relative to the searched span.
                let base = span.start;
                span.start = base + found.retained.start;
                span.end = base + found.retained.end;
                linked_matched = true;
            }
            linked = Some(outcome);
        }
        let skip_adapter = linked_matched
            || matches!(
                (self.linked.as_ref().map(|stage| stage.unmatched), linked),
                (
                    Some(crate::linked::Unmatched::Keep | crate::linked::Unmatched::Fail),
                    Some(crate::linked::LinkedOutcome::Unmatched(_))
                )
            );
        let adapter_hit = match (self.adapter.as_ref().filter(|_| !skip_adapter), insert) {
            // A successful overlap alignment fixes the insert boundary; explicit
            // adapter matching is skipped for the pair, because the boundary is
            // already known and an in-insert match would be a false positive.
            (Some(_), Some(insert)) => {
                let before = span.len();
                span.truncate_to(insert);
                overlap_removed = before - span.len();
                None
            }
            (Some(stage), None) => {
                let hit = stage.find(mate, span.sequence(read));
                if let Some(hit) = hit {
                    // `hit.start` is relative to the retained span, which
                    // still covers the whole read at this point.
                    span.end = span.start + hit.start;
                }
                hit
            }
            (None, _) => None,
        };
        let adapter_trimmed_length = span.len();

        let trimming = match &self.trim {
            Some(stage) => crate::trim::apply(stage.mate(mate), read, &mut span),
            None => TrimOutcome::default(),
        };

        let mut reasons = match &self.filter {
            Some(stage) => stage.evaluate(span.sequence(read), span.quality(read)),
            None => FilterReason::empty(),
        };
        // `--linked-unmatched fail` routes an unmatched read through the normal
        // rejection path rather than inventing a second one.
        if let (
            Some(crate::linked::Unmatched::Fail),
            Some(crate::linked::LinkedOutcome::Unmatched(_)),
        ) = (self.linked.as_ref().map(|stage| stage.unmatched), linked)
        {
            reasons |= FilterReason::LINKED_UNMATCHED;
        }

        debug_assert!(span.start <= span.end && span.end <= read.len());
        Ok(MateResult {
            retained: span,
            original_length: read.len(),
            adapter_trimmed_length,
            adapter_hit,
            linked,
            overlap_removed,
            trimming,
            reasons,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::adapter::{Adapter, AdapterParams};
    use crate::trim::{MateTrim, QualityCut};

    const ADAPTER: &[u8] = b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";

    fn adapter_stage() -> AdapterStage {
        AdapterStage::new(
            vec![Adapter::new("r1", ADAPTER).unwrap()],
            vec![Adapter::new("r2", b"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT").unwrap()],
            AdapterParams::default(),
            None,
        )
        .unwrap()
    }

    #[test]
    fn workflow_requires_at_least_one_stage() {
        assert!(Workflow::new(None, None, None, None, None).is_err());
        assert!(Workflow::new(Some(adapter_stage()), None, None, None, None).is_ok());
    }

    #[test]
    fn adapter_stage_shortens_the_read_at_the_match() {
        let workflow = Workflow::new(Some(adapter_stage()), None, None, None, None).unwrap();
        let sequence = [b"ACGTACGTACGTACGTACGT".as_slice(), ADAPTER].concat();
        let quality = vec![b'I'; sequence.len()];
        let read = ReadView::unchecked(&sequence, Some(&quality));
        let result = workflow.process(0, read, None).unwrap();
        assert_eq!(result.r1.retained, Span { start: 0, end: 20 });
        assert_eq!(result.r1.adapter_hit.unwrap().start, 20);
        assert_eq!(result.r1.adapter_bases_removed(), ADAPTER.len());
        assert!(
            result.accepted(),
            "no filter stage means nothing is rejected"
        );
    }

    #[test]
    fn per_mate_adapters_are_applied_independently() {
        let workflow = Workflow::new(Some(adapter_stage()), None, None, None, None).unwrap();
        // R1's adapter appears in R2 and vice versa: neither should match.
        let r2_only = [b"ACGTACGTACGTACGTACGT".as_slice(), ADAPTER].concat();
        let read = ReadView::unchecked(&r2_only, None);
        let result = workflow.process(0, read, Some(read)).unwrap();
        assert!(result.r1.adapter_hit.is_some());
        assert!(result.r2.unwrap().adapter_hit.is_none());
    }

    #[test]
    fn stages_run_in_declaration_order() {
        let trim = TrimStage {
            r1: MateTrim {
                tail: 4,
                ..MateTrim::default()
            },
            r2: MateTrim::default(),
        };
        let filter = FilterStage {
            min_length: Some(20),
            qualified_quality: 15,
            ..FilterStage::default()
        };
        let workflow =
            Workflow::new(Some(adapter_stage()), None, Some(trim), Some(filter), None).unwrap();
        assert_eq!(workflow.stage_order(), vec!["adapter", "trim", "filter"]);

        // 20 bases survive the adapter, then 4 are cut, leaving 16 (< 20).
        let sequence = [b"ACGTACGTACGTACGTACGT".as_slice(), ADAPTER].concat();
        let read = ReadView::unchecked(&sequence, None);
        let result = workflow.process(0, read, None).unwrap();
        assert_eq!(result.r1.adapter_trimmed_length, 20);
        assert_eq!(result.r1.final_length(), 16);
        assert_eq!(result.r1.reasons, FilterReason::TOO_SHORT);
        assert!(!result.accepted());
    }

    #[test]
    fn strict_pairing_rejects_the_pair_when_either_mate_fails() {
        let filter = FilterStage {
            min_length: Some(8),
            qualified_quality: 15,
            ..FilterStage::default()
        };
        let workflow = Workflow::new(None, None, None, Some(filter), None).unwrap();
        let long = ReadView::unchecked(b"ACGTACGTACGT", None);
        let short = ReadView::unchecked(b"ACGT", None);

        assert!(workflow.process(0, long, Some(long)).unwrap().accepted());

        let result = workflow.process(1, long, Some(short)).unwrap();
        assert!(!result.accepted());
        assert!(result.r1.passed());
        assert!(!result.r2.unwrap().passed());

        let result = workflow.process(2, short, Some(long)).unwrap();
        assert!(!result.accepted());
        assert!(!result.r1.passed());
        assert!(result.r2.unwrap().passed());

        let result = workflow.process(3, short, Some(short)).unwrap();
        assert!(!result.accepted());
        assert!(!result.r1.passed() && !result.r2.unwrap().passed());
    }

    #[test]
    fn orphan_pairing_routes_surviving_mates() {
        let filter = FilterStage {
            min_length: Some(8),
            qualified_quality: 15,
            ..FilterStage::default()
        };
        let workflow = Workflow::new(None, None, None, Some(filter), None)
            .unwrap()
            .with_pair_policy(PairPolicy::Orphan);
        let long = ReadView::unchecked(b"ACGTACGTACGT", None);
        let short = ReadView::unchecked(b"ACGT", None);

        let result = workflow.process(0, long, Some(long)).unwrap();
        assert_eq!(result.disposition, PairDisposition::Accepted);

        let result = workflow.process(1, long, Some(short)).unwrap();
        assert_eq!(result.disposition, PairDisposition::OrphanR1);

        let result = workflow.process(2, short, Some(long)).unwrap();
        assert_eq!(result.disposition, PairDisposition::OrphanR2);

        let result = workflow.process(3, short, Some(short)).unwrap();
        assert_eq!(result.disposition, PairDisposition::Rejected);

        // A single-end record can never become an orphan.
        let result = workflow.process(4, short, None).unwrap();
        assert_eq!(result.disposition, PairDisposition::Rejected);

        // The same configuration under strict pairing rejects broken pairs.
        let strict = Workflow::new(None, None, None, Some(filter), None).unwrap();
        let result = strict.process(5, long, Some(short)).unwrap();
        assert_eq!(result.disposition, PairDisposition::Rejected);
    }

    #[test]
    fn quality_dependent_workflows_reject_quality_free_records() {
        let trim = TrimStage {
            r1: MateTrim {
                quality_tail: Some(QualityCut {
                    minimum_phred: 20,
                    window: 4,
                }),
                ..MateTrim::default()
            },
            r2: MateTrim::default(),
        };
        let workflow = Workflow::new(None, None, Some(trim), None, None).unwrap();
        assert!(workflow.needs_quality());
        let read = ReadView::unchecked(b"ACGTACGT", None);
        assert!(matches!(
            workflow.process(0, read, None),
            Err(Error::MissingQuality(_))
        ));
    }

    #[test]
    fn invalid_quality_bytes_are_rejected_with_context() {
        let filter = FilterStage {
            min_mean_quality: Some(20),
            qualified_quality: 15,
            ..FilterStage::default()
        };
        let workflow = Workflow::new(None, None, None, Some(filter), None).unwrap();
        let read = ReadView::unchecked(b"ACGT", Some(b"II\x01I"));
        let err = workflow.process(11, read, None).unwrap_err();
        assert!(matches!(
            err,
            Error::InvalidQualityEncoding {
                record: 11,
                mate: "R1",
                byte: 1
            }
        ));
    }

    #[test]
    fn paired_overlap_trims_read_through_at_the_insert_boundary() {
        // A 100-base insert read by two 130-base mates. The insert is
        // deterministic but non-repetitive so the alignment is unambiguous.
        let mut rng = 0x1234_5678_9ABC_DEF1u64;
        let insert: Vec<u8> = (0..100)
            .map(|_| {
                rng ^= rng >> 12;
                rng ^= rng << 25;
                rng ^= rng >> 27;
                rng = rng.wrapping_mul(0x2545_F491_4F6C_DD1D);
                b"ACGT"[(rng >> 62) as usize]
            })
            .collect();
        let revcomp: Vec<u8> = insert
            .iter()
            .rev()
            .map(|&b| match b {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                _ => b'A',
            })
            .collect();
        let r1_seq = [insert.as_slice(), b"TTGGAACCGGTTCCAAGGTTGGAACTTAGG"].concat();
        let r2_seq = [revcomp.as_slice(), b"AACCTTGGAACCGGTTCCAAGGTTGGAACT"].concat();
        assert_eq!(r1_seq.len(), 130);
        assert_eq!(r2_seq.len(), 130);

        let stage = AdapterStage::new(
            Vec::new(),
            Vec::new(),
            AdapterParams::default(),
            Some(crate::overlap::OverlapParams::default()),
        )
        .unwrap();
        let workflow = Workflow::new(Some(stage), None, None, None, None).unwrap();
        let r1 = ReadView::unchecked(&r1_seq, None);
        let r2 = ReadView::unchecked(&r2_seq, None);
        let result = workflow.process(0, r1, Some(r2)).unwrap();
        assert_eq!(result.r1.final_length(), 100);
        assert_eq!(result.r1.overlap_removed, 30);
        assert_eq!(result.r2.unwrap().final_length(), 100);
        assert_eq!(result.r2.unwrap().overlap_removed, 30);
        // Overlap evidence replaces explicit matching for the pair.
        assert!(result.r1.adapter_hit.is_none());
        assert_eq!(result.r1.adapter_trimmed_length, 100);

        // Without an R2 there is no overlap signal and nothing happens.
        let result = workflow.process(1, r1, None).unwrap();
        assert_eq!(result.r1.final_length(), 130);
        assert_eq!(result.r1.overlap_removed, 0);
    }

    #[test]
    fn explicit_adapters_are_the_fallback_when_no_overlap_is_found() {
        // Unrelated mates: no overlap, but R1 carries the explicit adapter.
        let stage = AdapterStage::new(
            vec![Adapter::new("r1", ADAPTER).unwrap()],
            Vec::new(),
            AdapterParams::default(),
            Some(crate::overlap::OverlapParams::default()),
        )
        .unwrap();
        let workflow = Workflow::new(Some(stage), None, None, None, None).unwrap();
        let r1_seq = [b"ACGTACGTACGTACGTACGT".as_slice(), ADAPTER].concat();
        let r2_seq = vec![b'T'; 60];
        let r1 = ReadView::unchecked(&r1_seq, None);
        let r2 = ReadView::unchecked(&r2_seq, None);
        let result = workflow.process(0, r1, Some(r2)).unwrap();
        assert_eq!(result.r1.final_length(), 20);
        assert!(
            result.r1.adapter_hit.is_some(),
            "explicit match as fallback"
        );
        assert_eq!(result.r1.overlap_removed, 0);
        // R2 has no adapters configured and no overlap: untouched.
        assert_eq!(result.r2.unwrap().final_length(), 60);
    }

    #[test]
    fn spans_never_leave_the_original_read() {
        let trim = TrimStage {
            r1: MateTrim {
                front: 100,
                tail: 100,
                ..MateTrim::default()
            },
            r2: MateTrim::default(),
        };
        let workflow = Workflow::new(None, None, Some(trim), None, None).unwrap();
        let read = ReadView::unchecked(b"ACGT", None);
        let result = workflow.process(0, read, None).unwrap();
        assert_eq!(result.r1.retained.len(), 0);
        assert!(result.r1.retained.end <= 4);
    }
}
