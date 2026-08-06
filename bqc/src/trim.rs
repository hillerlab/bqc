// Copyright (c) 2026 The Hiller Lab at the Senckenberg Gesellschaft für Naturforschung
// Distributed under the terms of the GNU General Public License, Version 3.0.

//! Non-adapter shortening operations.
//!
//! Operations are applied in a fixed order, which is part of the `bqc`
//! output contract:
//!
//! ```text
//! fixed front/tail
//!   -> quality front / right / tail
//!   -> terminal N
//!   -> poly-G
//!   -> poly-X
//!   -> maximum-length truncation
//! ```
//!
//! Every operation updates the retained [`Span`] only; no sequence or quality
//! data is copied.

use serde::Serialize;

use crate::error::{Error, Result};
use crate::read::{phred_sum, ReadView, Span};

/// Number of tracked trimming operations.
pub const TRIM_OPS: usize = 6;

/// Identifies a trimming operation for statistics purposes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TrimOp {
    Fixed = 0,
    Quality = 1,
    TerminalN = 2,
    PolyG = 3,
    PolyX = 4,
    MaxLength = 5,
}

impl TrimOp {
    /// All operations in reporting order.
    pub const ALL: [TrimOp; TRIM_OPS] = [
        TrimOp::Fixed,
        TrimOp::Quality,
        TrimOp::TerminalN,
        TrimOp::PolyG,
        TrimOp::PolyX,
        TrimOp::MaxLength,
    ];

    /// Stable name used in reports.
    #[must_use]
    pub fn name(self) -> &'static str {
        match self {
            Self::Fixed => "fixed",
            Self::Quality => "quality",
            Self::TerminalN => "terminal_n",
            Self::PolyG => "poly_g",
            Self::PolyX => "poly_x",
            Self::MaxLength => "max_length",
        }
    }
}

/// Bases removed by each trimming operation for one mate.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct TrimOutcome {
    pub removed: [u32; TRIM_OPS],
}

impl TrimOutcome {
    /// Total bases removed by the trim stage.
    #[must_use]
    pub fn total(&self) -> u32 {
        self.removed.iter().sum()
    }

    #[inline]
    fn record(&mut self, op: TrimOp, before: usize, after: usize) {
        debug_assert!(after <= before);
        self.removed[op as usize] += (before - after) as u32;
    }
}

/// A windowed quality cut.
#[derive(Debug, Clone, Copy, Serialize)]
pub struct QualityCut {
    pub minimum_phred: u8,
    pub window: usize,
}

impl QualityCut {
    /// Validates a quality cut configuration.
    pub fn validate(self, flag: &str) -> Result<Self> {
        if self.window == 0 {
            return Err(Error::config(format!("--{flag}-window must be at least 1")));
        }
        Ok(self)
    }

    /// Whether `quality` meets the threshold on average.
    ///
    /// Uses an integer sum comparison rather than a floating-point mean:
    /// `sum >= threshold * len`.
    #[inline]
    #[must_use]
    fn window_qualifies(&self, quality: &[u8]) -> bool {
        phred_sum(quality) >= u32::from(self.minimum_phred) * quality.len() as u32
    }
}

/// Homopolymer tail trimming parameters.
#[derive(Debug, Clone, Copy, Serialize)]
pub struct PolyParams {
    pub min_length: usize,
    pub max_mismatch_rate: f64,
}

impl Default for PolyParams {
    fn default() -> Self {
        Self {
            min_length: 10,
            max_mismatch_rate: 0.10,
        }
    }
}

impl PolyParams {
    /// Validates poly-tail parameters.
    pub fn validate(self, flag: &str) -> Result<Self> {
        if self.min_length == 0 {
            return Err(Error::config(format!(
                "--{flag}-min-length must be at least 1"
            )));
        }
        if !(0.0..=1.0).contains(&self.max_mismatch_rate) {
            return Err(Error::config(format!(
                "--{flag}-max-mismatch-rate must be within 0.0..=1.0 (got {})",
                self.max_mismatch_rate
            )));
        }
        Ok(self)
    }
}

/// Resolved trimming configuration for one mate.
#[derive(Debug, Clone, Copy, Default, Serialize)]
pub struct MateTrim {
    pub front: usize,
    pub tail: usize,
    pub quality_front: Option<QualityCut>,
    pub quality_right: Option<QualityCut>,
    pub quality_tail: Option<QualityCut>,
    pub terminal_n: bool,
    pub poly_g: Option<PolyParams>,
    pub poly_x: Option<PolyParams>,
    pub max_length: Option<usize>,
}

impl MateTrim {
    /// Whether any operation is configured for this mate.
    #[must_use]
    pub fn is_noop(&self) -> bool {
        self.front == 0
            && self.tail == 0
            && self.quality_front.is_none()
            && self.quality_right.is_none()
            && self.quality_tail.is_none()
            && !self.terminal_n
            && self.poly_g.is_none()
            && self.poly_x.is_none()
            && self.max_length.is_none()
    }

    /// Whether this mate's configuration reads quality values.
    #[must_use]
    pub fn needs_quality(&self) -> bool {
        self.quality_front.is_some() || self.quality_right.is_some() || self.quality_tail.is_some()
    }
}

/// Compiled trim stage.
#[derive(Debug, Clone, Copy, Serialize)]
pub struct TrimStage {
    pub r1: MateTrim,
    pub r2: MateTrim,
}

impl TrimStage {
    /// Whether either mate reads quality values.
    #[must_use]
    pub fn needs_quality(&self) -> bool {
        self.r1.needs_quality() || self.r2.needs_quality()
    }

    /// Resolved configuration for `mate`.
    #[must_use]
    pub fn mate(&self, mate: crate::process::Mate) -> &MateTrim {
        match mate {
            crate::process::Mate::R1 => &self.r1,
            crate::process::Mate::R2 => &self.r2,
        }
    }
}

/// Applies every configured trimming operation to `span`, in contract order.
///
/// `read.quality` must be `Some` when [`MateTrim::needs_quality`] is true; the
/// caller is responsible for that check.
pub fn apply(config: &MateTrim, read: ReadView<'_>, span: &mut Span) -> TrimOutcome {
    let mut outcome = TrimOutcome::default();

    // 1. Fixed front/tail.
    if config.front > 0 || config.tail > 0 {
        let before = span.len();
        span.trim_front(config.front);
        span.trim_back(config.tail);
        outcome.record(TrimOp::Fixed, before, span.len());
    }

    // 2. Quality trimming: front, then right, then tail.
    if let Some(quality) = read.quality {
        let before = span.len();
        if let Some(cut) = config.quality_front {
            quality_front(&cut, quality, span);
        }
        if let Some(cut) = config.quality_right {
            quality_right(&cut, quality, span);
        }
        if let Some(cut) = config.quality_tail {
            quality_tail(&cut, quality, span);
        }
        outcome.record(TrimOp::Quality, before, span.len());
    }

    // 3. Terminal ambiguous bases.
    if config.terminal_n {
        let before = span.len();
        terminal_n(read.sequence, span);
        outcome.record(TrimOp::TerminalN, before, span.len());
    }

    // 4. Poly-G, then poly-X over whatever remains.
    if let Some(params) = config.poly_g {
        let before = span.len();
        poly_tail(&params, read.sequence, span, Some(b'G'));
        outcome.record(TrimOp::PolyG, before, span.len());
    }
    if let Some(params) = config.poly_x {
        let before = span.len();
        poly_tail(&params, read.sequence, span, None);
        outcome.record(TrimOp::PolyX, before, span.len());
    }

    // 5. Maximum-length truncation.
    if let Some(max_length) = config.max_length {
        let before = span.len();
        span.truncate_to(max_length);
        outcome.record(TrimOp::MaxLength, before, span.len());
    }

    outcome
}

/// Slides a window forward from the 5' end, dropping one base at a time until
/// the leading window qualifies.
fn quality_front(cut: &QualityCut, quality: &[u8], span: &mut Span) {
    while !span.is_empty() {
        let end = (span.start + cut.window).min(span.end);
        if cut.window_qualifies(&quality[span.start..end]) {
            break;
        }
        span.start += 1;
    }
}

/// Slides a window backward from the 3' end, dropping one base at a time until
/// the trailing window qualifies.
fn quality_tail(cut: &QualityCut, quality: &[u8], span: &mut Span) {
    while !span.is_empty() {
        let start = span.end.saturating_sub(cut.window).max(span.start);
        if cut.window_qualifies(&quality[start..span.end]) {
            break;
        }
        span.end -= 1;
    }
}

/// Scans windows left to right and truncates the read at the first window that
/// falls below the threshold.
fn quality_right(cut: &QualityCut, quality: &[u8], span: &mut Span) {
    if span.is_empty() {
        return;
    }
    let window = cut.window.min(span.len());
    let mut start = span.start;
    while start + window <= span.end {
        if !cut.window_qualifies(&quality[start..start + window]) {
            span.end = start;
            return;
        }
        start += 1;
    }
}

/// Removes contiguous `N` bases from both boundaries of the span.
fn terminal_n(sequence: &[u8], span: &mut Span) {
    while span.start < span.end && sequence[span.start] == b'N' {
        span.start += 1;
    }
    while span.end > span.start && sequence[span.end - 1] == b'N' {
        span.end -= 1;
    }
}

/// Canonical bases considered by poly-X trimming, in tie-breaking order.
const CANONICAL_BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];

/// Trims a homopolymer-like 3' tail.
///
/// The longest suffix of length `i` is removed for which `i >= min_length` and
/// `mismatches <= floor(i * max_mismatch_rate)`. When `base` is `None` every
/// canonical base is evaluated and the longest qualifying tail wins — which base
/// achieved it is not observable, since the outcome is a length.
fn poly_tail(params: &PolyParams, sequence: &[u8], span: &mut Span, base: Option<u8>) {
    let candidates: &[u8] = match &base {
        Some(base) => std::slice::from_ref(base),
        None => &CANONICAL_BASES,
    };
    let best = candidates
        .iter()
        .map(|&candidate| longest_poly_tail(params, sequence, *span, candidate))
        .max()
        .unwrap_or(0);
    span.trim_back(best);
}

/// Length of the longest qualifying homopolymer tail of `base`.
///
/// The scan stops as soon as no longer suffix can qualify. Mismatches only grow
/// as the suffix lengthens, and the largest allowance any suffix can have is
/// `span_len * rate`, so once the count passes that the remaining positions
/// cannot produce a match. On reads without a homopolymer tail this ends the scan
/// after a couple of dozen bases instead of walking the whole read, four times
/// over for poly-X.
///
/// The tolerance test is written as `mismatches <= length * rate` rather than
/// `mismatches <= floor(length * rate)`. For an integer `m`, `m <= floor(x)`
/// holds exactly when `m <= x`, so this is the same predicate without a libm
/// `floor` call, which does not inline.
fn longest_poly_tail(params: &PolyParams, sequence: &[u8], span: Span, base: u8) -> usize {
    let ceiling = span.len() as f64 * params.max_mismatch_rate;
    let mut mismatches = 0usize;
    let mut best = 0usize;
    for (examined, position) in (span.start..span.end).rev().enumerate() {
        if sequence[position] != base {
            mismatches += 1;
            if mismatches as f64 > ceiling {
                break;
            }
        }
        let length = examined + 1;
        if length >= params.min_length
            && mismatches as f64 <= length as f64 * params.max_mismatch_rate
        {
            best = length;
        }
    }
    best
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::read::phred;

    fn view<'a>(sequence: &'a [u8], quality: &'a [u8]) -> ReadView<'a> {
        ReadView::unchecked(sequence, Some(quality))
    }

    fn trimmed<'a>(config: &MateTrim, read: ReadView<'a>) -> (&'a [u8], TrimOutcome) {
        let mut span = Span::full(read.len());
        let outcome = apply(config, read, &mut span);
        (span.sequence(read), outcome)
    }

    #[test]
    fn fixed_trimming_removes_both_ends() {
        let read = view(b"AAAACCCCGGGG", b"IIIIIIIIIIII");
        let config = MateTrim {
            front: 4,
            tail: 4,
            ..MateTrim::default()
        };
        let (sequence, outcome) = trimmed(&config, read);
        assert_eq!(sequence, b"CCCC");
        assert_eq!(outcome.removed[TrimOp::Fixed as usize], 8);
    }

    #[test]
    fn fixed_trimming_exceeding_length_yields_empty_read() {
        let read = view(b"ACGT", b"IIII");
        let config = MateTrim {
            front: 10,
            tail: 10,
            ..MateTrim::default()
        };
        let (sequence, outcome) = trimmed(&config, read);
        assert!(sequence.is_empty());
        assert_eq!(
            outcome.removed[TrimOp::Fixed as usize],
            4,
            "never over-counts"
        );
    }

    #[test]
    fn fixed_front_and_tail_meeting_at_one_coordinate() {
        let read = view(b"ACGTAC", b"IIIIII");
        let config = MateTrim {
            front: 3,
            tail: 3,
            ..MateTrim::default()
        };
        let (sequence, _) = trimmed(&config, read);
        assert!(sequence.is_empty());
    }

    #[test]
    fn quality_front_stops_at_first_qualifying_window() {
        // '!' == Q0, 'I' == Q40. Window 4, threshold Q20 (window sum 80).
        // [0..4] = "!!!!" -> 0, [1..5] = "!!!I" -> 40, [2..6] = "!!II" -> 80,
        // which meets the threshold exactly, so two bases are removed.
        let read = view(b"AAAACCCCGGGG", b"!!!!IIIIIIII");
        let config = MateTrim {
            quality_front: Some(QualityCut {
                minimum_phred: 20,
                window: 4,
            }),
            ..MateTrim::default()
        };
        let (sequence, outcome) = trimmed(&config, read);
        assert_eq!(sequence, b"AACCCCGGGG");
        assert_eq!(outcome.removed[TrimOp::Quality as usize], 2);
    }

    #[test]
    fn quality_tail_stops_at_first_qualifying_window() {
        // Mirror image of the front case: [6..10] = "II!!" -> 80 qualifies.
        let read = view(b"AAAACCCCGGGG", b"IIIIIIII!!!!");
        let config = MateTrim {
            quality_tail: Some(QualityCut {
                minimum_phred: 20,
                window: 4,
            }),
            ..MateTrim::default()
        };
        let (sequence, outcome) = trimmed(&config, read);
        assert_eq!(sequence, b"AAAACCCCGG");
        assert_eq!(outcome.removed[TrimOp::Quality as usize], 2);
    }

    #[test]
    fn quality_window_at_exact_threshold_is_retained() {
        // '5' == Q20; the window mean equals the threshold exactly.
        let read = view(b"ACGT", b"5555");
        let config = MateTrim {
            quality_tail: Some(QualityCut {
                minimum_phred: 20,
                window: 4,
            }),
            ..MateTrim::default()
        };
        assert_eq!(trimmed(&config, read).0, b"ACGT");

        // One point below the threshold fails.
        let read = view(b"ACGT", b"5554");
        assert!(trimmed(&config, read).0.len() < 4);
    }

    #[test]
    fn quality_window_shorter_than_configured_window_uses_remaining_bases() {
        // Only two bases remain, both Q40, with a configured window of 8.
        let read = view(b"AC", b"II");
        let config = MateTrim {
            quality_tail: Some(QualityCut {
                minimum_phred: 20,
                window: 8,
            }),
            ..MateTrim::default()
        };
        assert_eq!(trimmed(&config, read).0, b"AC");

        // Two low-quality bases are fully removed.
        let read = view(b"AC", b"!!");
        assert!(trimmed(&config, read).0.is_empty());
    }

    #[test]
    fn quality_right_truncates_at_the_failing_window() {
        // Window 4, threshold Q20 (window sum 80). Scanning left to right the
        // first failing window is [3..7] = "I!!!" -> 40, so the read is cut at
        // offset 3 even though the low-quality run starts at offset 4.
        let read = view(b"AAAACCCCGGGG", b"IIII!!!!IIII");
        let config = MateTrim {
            quality_right: Some(QualityCut {
                minimum_phred: 20,
                window: 4,
            }),
            ..MateTrim::default()
        };
        let (sequence, _) = trimmed(&config, read);
        assert_eq!(sequence, b"AAA");
    }

    #[test]
    fn quality_right_keeps_reads_that_never_fail() {
        let read = view(b"AAAACCCC", b"IIIIIIII");
        let config = MateTrim {
            quality_right: Some(QualityCut {
                minimum_phred: 20,
                window: 4,
            }),
            ..MateTrim::default()
        };
        assert_eq!(trimmed(&config, read).0, b"AAAACCCC");
    }

    #[test]
    fn quality_trimming_can_remove_the_entire_read() {
        let read = view(b"ACGTACGT", b"!!!!!!!!");
        let config = MateTrim {
            quality_tail: Some(QualityCut {
                minimum_phred: 20,
                window: 4,
            }),
            ..MateTrim::default()
        };
        let (sequence, outcome) = trimmed(&config, read);
        assert!(sequence.is_empty());
        assert_eq!(outcome.removed[TrimOp::Quality as usize], 8);
    }

    #[test]
    fn terminal_n_trims_both_ends_but_not_internal_ns() {
        let read = view(b"NNACGNTACNN", b"IIIIIIIIIII");
        let config = MateTrim {
            terminal_n: true,
            ..MateTrim::default()
        };
        let (sequence, outcome) = trimmed(&config, read);
        assert_eq!(sequence, b"ACGNTAC");
        assert_eq!(outcome.removed[TrimOp::TerminalN as usize], 4);
    }

    #[test]
    fn terminal_n_only_read_is_fully_removed() {
        let read = view(b"NNNN", b"IIII");
        let config = MateTrim {
            terminal_n: true,
            ..MateTrim::default()
        };
        assert!(trimmed(&config, read).0.is_empty());
    }

    #[test]
    fn poly_g_trims_a_g_tail_with_a_tolerated_mismatch() {
        // 12-base tail with one mismatch: floor(12 * 0.1) == 1 is tolerated.
        let read = view(b"ACGTACGTACGTGGGGGGAGGGGG", b"IIIIIIIIIIIIIIIIIIIIIIII");
        let config = MateTrim {
            poly_g: Some(PolyParams::default()),
            ..MateTrim::default()
        };
        let (sequence, outcome) = trimmed(&config, read);
        assert_eq!(sequence, b"ACGTACGTACGT");
        assert_eq!(outcome.removed[TrimOp::PolyG as usize], 12);
    }

    #[test]
    fn poly_g_respects_the_minimum_length() {
        let read = view(b"ACGTACGTGGGGG", b"IIIIIIIIIIIII");
        let config = MateTrim {
            poly_g: Some(PolyParams {
                min_length: 10,
                max_mismatch_rate: 0.1,
            }),
            ..MateTrim::default()
        };
        assert_eq!(
            trimmed(&config, read).0,
            b"ACGTACGTGGGGG",
            "5 G's is below min_length"
        );

        let config = MateTrim {
            poly_g: Some(PolyParams {
                min_length: 5,
                max_mismatch_rate: 0.0,
            }),
            ..MateTrim::default()
        };
        assert_eq!(trimmed(&config, read).0, b"ACGTACGT");
    }

    #[test]
    fn poly_x_examines_every_canonical_base() {
        let read = view(b"GGGGCCCCAAAAAAAAAA", b"IIIIIIIIIIIIIIIIII");
        let config = MateTrim {
            poly_x: Some(PolyParams {
                min_length: 10,
                max_mismatch_rate: 0.0,
            }),
            ..MateTrim::default()
        };
        assert_eq!(trimmed(&config, read).0, b"GGGGCCCC");
    }

    #[test]
    fn poly_x_ties_prefer_the_earlier_base_but_pick_the_longest_tail() {
        // An all-T read: only T qualifies, and the whole read is the tail.
        let read = view(b"TTTTTTTTTT", b"IIIIIIIIII");
        let config = MateTrim {
            poly_x: Some(PolyParams {
                min_length: 10,
                max_mismatch_rate: 0.0,
            }),
            ..MateTrim::default()
        };
        assert!(trimmed(&config, read).0.is_empty());

        // With a 50% mismatch budget an alternating tail qualifies for both
        // bases at the same length; the outcome is identical either way.
        let read = view(b"ACACACACAC", b"IIIIIIIIII");
        let config = MateTrim {
            poly_x: Some(PolyParams {
                min_length: 10,
                max_mismatch_rate: 0.5,
            }),
            ..MateTrim::default()
        };
        assert!(trimmed(&config, read).0.is_empty());
    }

    #[test]
    fn poly_g_runs_before_poly_x_without_double_counting() {
        // G tail followed by nothing; poly-x then sees the A tail beneath it.
        let read = view(b"CCCCAAAAAAAAAAGGGGGGGGGG", b"IIIIIIIIIIIIIIIIIIIIIIII");
        let config = MateTrim {
            poly_g: Some(PolyParams {
                min_length: 10,
                max_mismatch_rate: 0.0,
            }),
            poly_x: Some(PolyParams {
                min_length: 10,
                max_mismatch_rate: 0.0,
            }),
            ..MateTrim::default()
        };
        let (sequence, outcome) = trimmed(&config, read);
        assert_eq!(sequence, b"CCCC");
        assert_eq!(outcome.removed[TrimOp::PolyG as usize], 10);
        assert_eq!(outcome.removed[TrimOp::PolyX as usize], 10);
        assert_eq!(outcome.total(), 20);
    }

    #[test]
    fn max_length_truncation_keeps_the_leading_bases() {
        let read = view(b"ACGTACGTACGT", b"IIIIIIIIIIII");
        let config = MateTrim {
            max_length: Some(5),
            ..MateTrim::default()
        };
        let (sequence, outcome) = trimmed(&config, read);
        assert_eq!(sequence, b"ACGTA");
        assert_eq!(outcome.removed[TrimOp::MaxLength as usize], 7);
    }

    #[test]
    fn max_length_truncation_is_applied_after_other_operations() {
        // The front cut happens first, so truncation counts from the new start.
        let read = view(b"AAAACCCCGGGG", b"IIIIIIIIIIII");
        let config = MateTrim {
            front: 4,
            max_length: Some(4),
            ..MateTrim::default()
        };
        assert_eq!(trimmed(&config, read).0, b"CCCC");
    }

    #[test]
    fn operations_compose_in_contract_order() {
        // Layout: 2 bases cut by --front, ACGT, a 10 base poly-G tail, 2 N's
        // and 2 low-quality bases. Each stage exposes work for the next one.
        let read = view(b"TTACGTGGGGGGGGGGNNAA", b"IIIIIIIIIIIIIIIIII!!");
        let config = MateTrim {
            front: 2,
            quality_tail: Some(QualityCut {
                minimum_phred: 20,
                window: 1,
            }),
            terminal_n: true,
            poly_g: Some(PolyParams {
                min_length: 10,
                max_mismatch_rate: 0.0,
            }),
            max_length: Some(4),
            ..MateTrim::default()
        };
        let (sequence, outcome) = trimmed(&config, read);
        assert_eq!(sequence, b"ACGT");
        assert_eq!(outcome.removed[TrimOp::Fixed as usize], 2);
        assert_eq!(outcome.removed[TrimOp::Quality as usize], 2);
        assert_eq!(outcome.removed[TrimOp::TerminalN as usize], 2);
        assert_eq!(outcome.removed[TrimOp::PolyG as usize], 10);
        assert_eq!(outcome.removed[TrimOp::MaxLength as usize], 0);
    }

    #[test]
    fn quality_operations_are_skipped_when_quality_is_absent() {
        let read = ReadView::unchecked(b"ACGTACGT", None);
        let config = MateTrim {
            quality_tail: Some(QualityCut {
                minimum_phred: 40,
                window: 4,
            }),
            ..MateTrim::default()
        };
        let mut span = Span::full(read.len());
        let outcome = apply(&config, read, &mut span);
        assert_eq!(span.len(), 8);
        assert_eq!(outcome.total(), 0);
    }

    #[test]
    fn empty_reads_are_handled_by_every_operation() {
        let read = view(b"", b"");
        let config = MateTrim {
            front: 3,
            tail: 3,
            quality_front: Some(QualityCut {
                minimum_phred: 30,
                window: 4,
            }),
            quality_right: Some(QualityCut {
                minimum_phred: 30,
                window: 4,
            }),
            terminal_n: true,
            poly_g: Some(PolyParams::default()),
            poly_x: Some(PolyParams::default()),
            max_length: Some(2),
            ..MateTrim::default()
        };
        let mut span = Span::full(0);
        let outcome = apply(&config, read, &mut span);
        assert_eq!(span, Span { start: 0, end: 0 });
        assert_eq!(outcome.total(), 0);
    }

    #[test]
    fn validation_rejects_degenerate_parameters() {
        assert!(QualityCut {
            minimum_phred: 20,
            window: 0
        }
        .validate("quality-tail")
        .is_err());
        assert!(PolyParams {
            min_length: 0,
            max_mismatch_rate: 0.1
        }
        .validate("poly-g")
        .is_err());
        assert!(PolyParams {
            min_length: 5,
            max_mismatch_rate: 1.5
        }
        .validate("poly-g")
        .is_err());
        assert!(PolyParams {
            min_length: 5,
            max_mismatch_rate: f64::NAN
        }
        .validate("poly-g")
        .is_err());
    }

    #[test]
    fn phred_is_the_only_decoding_path() {
        // Guards against accidental re-implementation of the offset.
        assert_eq!(phred(b'5'), 20);
        let cut = QualityCut {
            minimum_phred: 20,
            window: 4,
        };
        assert!(cut.window_qualifies(b"5555"));
        assert!(!cut.window_qualifies(b"4444"));
    }
}
