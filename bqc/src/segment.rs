// Copyright (c) 2026 The Hiller Lab at the Senckenberg Gesellschaft für Naturforschung
// Distributed under the terms of the GNU General Public License, Version 3.0.

//! Arbitrary internal adapter splitting.
//!
//! Adapters act as delimiters rather than as ends to trim from:
//!
//! ```text
//! PREFIX [ADAPTER A] SEGMENT 1 [ADAPTER B] SEGMENT 2 [ADAPTER C] SUFFIX
//!    ->  PREFIX, SEGMENT 1, SEGMENT 2, SUFFIX
//! ```
//!
//! This changes output cardinality — one input record becomes zero, one or many
//! output records — which is why it lives behind its own `segment` command
//! instead of an option on `adapter`.
//!
//! # Pipeline
//!
//! ```text
//! find every candidate adapter match
//!   -> select a non-overlapping subset
//!   -> cut the read at those boundaries
//!   -> trim each fragment
//!   -> filter each fragment
//! ```
//!
//! # Candidate selection
//!
//! Candidates are ordered by start coordinate, then fewer edits, then lower error
//! rate, then greater overlap, then adapter declaration order. They are then
//! accepted greedily from the left, and any candidate overlapping an already
//! accepted delimiter is suppressed. That is deliberately the simplest rule that
//! is deterministic; weighted interval scheduling is not used unless a real case
//! shows the simple rule segmenting wrongly.

use serde::{Deserialize, Serialize};

use crate::adapter::{Adapter, AdapterHit, AdapterParams, verify_at};
use crate::error::{Error, Result};
use crate::read::Span;

/// Safety limit on fragments emitted from one read.
pub const DEFAULT_MAX_SEGMENTS: usize = 64;

/// What happens to the fragments at the ends of a read.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize, clap::ValueEnum)]
#[serde(rename_all = "lowercase")]
pub enum Terminal {
    /// Keep them (default).
    #[default]
    Keep,
    /// Keep only fragments with an adapter on both sides.
    Discard,
}

/// Compiled segmentation configuration.
#[derive(Debug, Clone, Serialize)]
pub struct SegmentStage {
    pub adapters: Vec<Adapter>,
    pub params: AdapterParams,
    pub terminal: Terminal,
    pub min_segment_length: usize,
    pub max_segments: usize,
}

impl SegmentStage {
    /// Builds the stage, rejecting configurations that cannot segment.
    pub fn new(
        adapters: Vec<Adapter>,
        params: AdapterParams,
        terminal: Terminal,
        min_segment_length: usize,
        max_segments: usize,
    ) -> Result<Self> {
        if adapters.is_empty() {
            return Err(Error::config(
                "segmentation requires --adapter-r1 or --adapter-fasta",
            ));
        }
        if max_segments == 0 {
            return Err(Error::config("--max-segments-per-read must be at least 1"));
        }
        for adapter in &adapters {
            if adapter.len() < params.min_overlap {
                return Err(Error::InvalidAdapter(format!(
                    "adapter '{}' is {} bases long, shorter than --min-overlap {}",
                    adapter.name,
                    adapter.len(),
                    params.min_overlap
                )));
            }
        }
        Ok(Self {
            adapters,
            params,
            terminal,
            min_segment_length,
            max_segments,
        })
    }
}

/// One emitted fragment.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Fragment {
    /// Position within the source read.
    pub span: Span,
    /// Zero-based index among the fragments this read emitted.
    pub index: usize,
    /// Adapter delimiting the fragment on the left, if any.
    pub left_adapter: Option<usize>,
    /// Adapter delimiting the fragment on the right, if any.
    pub right_adapter: Option<usize>,
}

impl Fragment {
    /// Whether the fragment has an adapter on both sides.
    #[must_use]
    pub fn internal(&self) -> bool {
        self.left_adapter.is_some() && self.right_adapter.is_some()
    }
}

/// What segmenting one read produced.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct SegmentOutcome {
    /// Candidate matches found before suppression.
    pub candidates: usize,
    /// Delimiters accepted.
    pub boundaries: usize,
    /// Candidates dropped for overlapping an accepted delimiter.
    pub suppressed: usize,
    /// Fragments dropped for being empty.
    pub empty: usize,
    /// Fragments dropped for being shorter than the minimum.
    pub too_short: usize,
    /// Fragments dropped because they were terminal and terminals are discarded.
    pub terminal_discarded: usize,
    /// Fragments dropped because the read hit the safety limit.
    pub over_limit: usize,
    /// Terminal fragments emitted.
    pub terminal: usize,
    /// Internal fragments emitted.
    pub internal: usize,
}

impl SegmentOutcome {
    /// Fragments emitted for filtering.
    #[must_use]
    pub fn emitted(&self) -> usize {
        self.terminal + self.internal
    }
}

/// Reused worker buffers. Segmentation never allocates in the record loop.
#[derive(Debug, Default)]
pub struct SegmentScratch {
    // `hits` holds the accepted delimiters once `segment` returns.
    hits: Vec<AdapterHit>,
    /// Banded-alignment workspace, only ever filled under `--allow-indels`.
    dp: Vec<u32>,
    pub fragments: Vec<Fragment>,
}

impl SegmentScratch {
    /// The delimiters accepted for the last segmented read.
    #[must_use]
    pub fn delimiters(&self) -> &[AdapterHit] {
        &self.hits
    }
}

/// Collects every candidate adapter match in `sequence`.
///
/// Unlike [`crate::adapter::find_three_prime`], which returns one best hit, this
/// reports every coordinate any adapter accepts, because a delimiter can occur
/// anywhere and several can occur in one read. `output` and `scratch` are
/// worker-owned buffers, so no allocation happens per read. The best-hit path is
/// untouched: ordinary trimming never collects candidates.
///
/// A match away from the read's 3' end is necessarily a whole adapter, because
/// [`verify_at`] only shortens the compared window when the read runs out. That
/// is the property that makes an accepted hit usable as a delimiter.
pub fn find_all(
    adapters: &[Adapter],
    params: AdapterParams,
    sequence: &[u8],
    scratch: &mut Vec<u32>,
    output: &mut Vec<AdapterHit>,
) {
    output.clear();
    let len = sequence.len();
    if len < params.min_overlap {
        return;
    }
    if params.allow_indels {
        scratch.clear();
        scratch.resize(2 * (len + 1), u32::MAX);
    }
    for start in 0..=len - params.min_overlap {
        if params.allow_indels {
            // The banded matcher reports the best alignment at a coordinate, not
            // one per adapter; ties there are already broken deterministically.
            if let Some(hit) =
                crate::adapter::indel_best_at(adapters, params, sequence, start, scratch)
            {
                output.push(hit);
            }
            continue;
        }
        for (index, adapter) in adapters.iter().enumerate() {
            if let Some(mut hit) = verify_at(adapter, params, sequence, start) {
                hit.adapter = index;
                output.push(hit);
            }
        }
    }
}

/// Orders candidates and suppresses those overlapping an accepted delimiter.
///
/// Returns the number suppressed. The surviving delimiters are left in `hits`,
/// sorted by coordinate and mutually non-overlapping.
pub fn select(hits: &mut Vec<AdapterHit>) -> usize {
    let before = hits.len();
    // Earliest coordinate first, then the strongest match at that coordinate.
    hits.sort_by(|a, b| {
        let key = |hit: &AdapterHit| {
            (
                hit.start,
                hit.errors,
                (hit.errors as f64 / hit.overlap.max(1) as f64 * 1e9) as u64,
                std::cmp::Reverse(hit.overlap),
                hit.adapter,
            )
        };
        key(a).cmp(&key(b))
    });
    let mut accepted_end = 0usize;
    let mut first = true;
    hits.retain(|hit| {
        // A delimiter occupies `start .. start + consumed`; anything reaching into
        // that interval describes the same adapter copy and is dropped.
        if first || hit.start >= accepted_end {
            accepted_end = hit.start + hit.consumed;
            first = false;
            true
        } else {
            false
        }
    });
    before - hits.len()
}

/// Cuts `length` bases at the accepted delimiters, appending the fragments.
fn cut(
    stage: &SegmentStage,
    length: usize,
    hits: &[AdapterHit],
    fragments: &mut Vec<Fragment>,
) -> SegmentOutcome {
    let mut outcome = SegmentOutcome {
        boundaries: hits.len(),
        ..SegmentOutcome::default()
    };
    let mut cursor = 0usize;
    let mut left_adapter = None;
    // Each delimiter closes the fragment before it and opens the next.
    for hit in hits.iter().chain(std::iter::once(&AdapterHit {
        adapter: usize::MAX,
        start: length,
        errors: 0,
        overlap: 0,
        consumed: 0,
    })) {
        let terminator = hit.adapter != usize::MAX;
        let span = Span {
            start: cursor,
            end: hit.start,
        };
        let candidate = Fragment {
            span,
            index: fragments.len(),
            left_adapter,
            right_adapter: terminator.then_some(hit.adapter),
        };
        push_fragment(stage, candidate, fragments, &mut outcome);
        cursor = hit.start + hit.consumed;
        left_adapter = Some(hit.adapter);
    }
    outcome
}

/// Applies the emission rules to one candidate fragment.
fn push_fragment(
    stage: &SegmentStage,
    mut candidate: Fragment,
    fragments: &mut Vec<Fragment>,
    outcome: &mut SegmentOutcome,
) {
    if candidate.span.is_empty() {
        // Adjacent adapters, or one at a read end.
        outcome.empty += 1;
        return;
    }
    if stage.terminal == Terminal::Discard && !candidate.internal() {
        outcome.terminal_discarded += 1;
        return;
    }
    if candidate.span.len() < stage.min_segment_length {
        outcome.too_short += 1;
        return;
    }
    if fragments.len() >= stage.max_segments {
        outcome.over_limit += 1;
        return;
    }
    candidate.index = fragments.len();
    if candidate.internal() {
        outcome.internal += 1;
    } else {
        outcome.terminal += 1;
    }
    fragments.push(candidate);
}

/// Segments one read, leaving the fragments in `scratch`.
pub fn segment(
    stage: &SegmentStage,
    sequence: &[u8],
    scratch: &mut SegmentScratch,
) -> SegmentOutcome {
    scratch.fragments.clear();
    find_all(
        &stage.adapters,
        stage.params,
        sequence,
        &mut scratch.dp,
        &mut scratch.hits,
    );
    let candidates = scratch.hits.len();
    let suppressed = select(&mut scratch.hits);
    let mut outcome = cut(stage, sequence.len(), &scratch.hits, &mut scratch.fragments);
    outcome.candidates = candidates;
    outcome.suppressed = suppressed;
    outcome
}

/// The header a fragment carries, derived from its source header.
///
/// Centralised so provenance formatting exists in exactly one place.
pub fn fragment_header(out: &mut Vec<u8>, source: &[u8], fragment: Fragment) {
    use std::io::Write as _;
    out.clear();
    out.extend_from_slice(source);
    let _ = write!(
        out,
        "|segment={}|span={}-{}",
        fragment.index, fragment.span.start, fragment.span.end
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    const A: &[u8] = b"ACGTACGTAA";
    const B: &[u8] = b"TTGGCCAATT";

    fn params() -> AdapterParams {
        AdapterParams {
            min_overlap: 8,
            max_error_rate: 0.10,
            max_errors: None,
            allow_indels: false,
        }
    }

    fn stage(adapters: &[&[u8]]) -> SegmentStage {
        let adapters = adapters
            .iter()
            .enumerate()
            .map(|(index, sequence)| Adapter::new(format!("a{index}"), sequence).unwrap())
            .collect();
        SegmentStage::new(adapters, params(), Terminal::Keep, 1, DEFAULT_MAX_SEGMENTS).unwrap()
    }

    fn run(stage: &SegmentStage, sequence: &[u8]) -> (Vec<Fragment>, SegmentOutcome) {
        let mut scratch = SegmentScratch::default();
        let outcome = segment(stage, sequence, &mut scratch);
        (scratch.fragments.clone(), outcome)
    }

    fn pieces(sequence: &[u8], fragments: &[Fragment]) -> Vec<Vec<u8>> {
        fragments
            .iter()
            .map(|fragment| sequence[fragment.span.start..fragment.span.end].to_vec())
            .collect()
    }

    #[test]
    fn a_read_without_an_adapter_yields_one_fragment() {
        let sequence = b"GGGGGTTTTTCCCCC";
        let (fragments, outcome) = run(&stage(&[A]), sequence);
        assert_eq!(pieces(sequence, &fragments), vec![sequence.to_vec()]);
        assert_eq!(outcome.boundaries, 0);
        assert_eq!(outcome.terminal, 1);
        assert_eq!(outcome.internal, 0);
    }

    #[test]
    fn one_adapter_splits_the_read_in_two() {
        let sequence = [b"GGGGG".as_slice(), A, b"TTTTT"].concat();
        let (fragments, outcome) = run(&stage(&[A]), &sequence);
        assert_eq!(
            pieces(&sequence, &fragments),
            vec![b"GGGGG".to_vec(), b"TTTTT".to_vec()]
        );
        assert_eq!(outcome.boundaries, 1);
        assert_eq!(outcome.terminal, 2, "both pieces touch a read end");
        assert_eq!(fragments[0].right_adapter, Some(0));
        assert_eq!(fragments[1].left_adapter, Some(0));
        assert_eq!(fragments[0].index, 0);
        assert_eq!(fragments[1].index, 1);
    }

    #[test]
    fn several_adapters_produce_prefix_middles_and_suffix() {
        let sequence = [
            b"GGGGG".as_slice(),
            A,
            b"TTTTTTT",
            B,
            b"CCCCCCCC",
            A,
            b"AAAAA",
        ]
        .concat();
        let (fragments, outcome) = run(&stage(&[A, B]), &sequence);
        assert_eq!(
            pieces(&sequence, &fragments),
            vec![
                b"GGGGG".to_vec(),
                b"TTTTTTT".to_vec(),
                b"CCCCCCCC".to_vec(),
                b"AAAAA".to_vec(),
            ]
        );
        assert_eq!(outcome.boundaries, 3);
        assert_eq!(outcome.internal, 2, "the two middles");
        assert_eq!(outcome.terminal, 2, "prefix and suffix");
        assert_eq!(fragments[1].left_adapter, Some(0));
        assert_eq!(fragments[1].right_adapter, Some(1));
    }

    #[test]
    fn an_adapter_at_the_first_or_last_base_yields_no_empty_fragment() {
        let leading = [A, b"TTTTT".as_slice()].concat();
        let (fragments, outcome) = run(&stage(&[A]), &leading);
        assert_eq!(pieces(&leading, &fragments), vec![b"TTTTT".to_vec()]);
        assert_eq!(outcome.empty, 1, "the zero-length prefix is discarded");

        let trailing = [b"TTTTT".as_slice(), A].concat();
        let (fragments, outcome) = run(&stage(&[A]), &trailing);
        assert_eq!(pieces(&trailing, &fragments), vec![b"TTTTT".to_vec()]);
        assert_eq!(outcome.empty, 1);
    }

    #[test]
    fn adjacent_adapters_produce_no_empty_fragment_between_them() {
        let sequence = [b"GGGGG".as_slice(), A, B, b"TTTTT"].concat();
        let (fragments, outcome) = run(&stage(&[A, B]), &sequence);
        assert_eq!(
            pieces(&sequence, &fragments),
            vec![b"GGGGG".to_vec(), b"TTTTT".to_vec()]
        );
        assert_eq!(outcome.boundaries, 2);
        assert_eq!(outcome.empty, 1, "nothing lies between the two adapters");
    }

    #[test]
    fn a_repeated_adapter_is_found_at_every_occurrence() {
        let sequence = [b"GG".as_slice(), A, b"TTTTT", A, b"CCC"].concat();
        let (fragments, outcome) = run(&stage(&[A]), &sequence);
        assert_eq!(
            pieces(&sequence, &fragments),
            vec![b"GG".to_vec(), b"TTTTT".to_vec(), b"CCC".to_vec()]
        );
        assert_eq!(outcome.boundaries, 2);
    }

    #[test]
    fn overlapping_candidates_are_suppressed_deterministically() {
        // A self-similar adapter inside a longer repeat matches at several
        // coordinates. The earliest wins; the rest describe the same copy shifted
        // by a period and are suppressed.
        let stage = stage(&[b"ATATATATATAT"]);
        let sequence = [b"GGGG".as_slice(), b"ATATATATATATATAT", b"CCCC"].concat();
        let (fragments, outcome) = run(&stage, &sequence);
        assert_eq!(outcome.candidates, 3, "matches at offsets 4, 6 and 8");
        assert_eq!(outcome.boundaries, 1, "one delimiter survives");
        assert_eq!(outcome.suppressed, 2);
        assert_eq!(
            pieces(&sequence, &fragments),
            vec![b"GGGG".to_vec(), b"ATATCCCC".to_vec()],
            "the leftmost copy is cut out; the repeat's tail stays in the suffix"
        );
    }

    #[test]
    fn competing_adapters_at_one_coordinate_break_by_errors_then_order() {
        let mut noisy = A.to_vec();
        noisy[9] = if noisy[9] == b'A' { b'C' } else { b'A' };
        let sequence = [b"GGGGG".as_slice(), A, b"TTTTT"].concat();
        // The exact adapter is declared second but matches with no errors.
        let stage = SegmentStage::new(
            vec![
                Adapter::new("noisy", &noisy).unwrap(),
                Adapter::new("exact", A).unwrap(),
            ],
            params(),
            Terminal::Keep,
            1,
            DEFAULT_MAX_SEGMENTS,
        )
        .unwrap();
        let (fragments, _) = run(&stage, &sequence);
        assert_eq!(fragments[0].right_adapter, Some(1), "fewest errors wins");
    }

    #[test]
    fn a_delimiter_with_one_mismatch_still_cuts() {
        // The budget over ten bases is one substitution.
        let mut damaged = A.to_vec();
        damaged[5] = if damaged[5] == b'C' { b'G' } else { b'C' };
        let sequence = [b"GGGGG".as_slice(), &damaged, b"TTTTT"].concat();
        let (fragments, outcome) = run(&stage(&[A]), &sequence);
        assert_eq!(outcome.boundaries, 1);
        assert_eq!(
            pieces(&sequence, &fragments),
            vec![b"GGGGG".to_vec(), b"TTTTT".to_vec()]
        );

        // Two substitutions exceed it, and the read stays whole.
        damaged[6] = if damaged[6] == b'C' { b'G' } else { b'C' };
        let sequence = [b"GGGGG".as_slice(), &damaged, b"TTTTT"].concat();
        let (fragments, outcome) = run(&stage(&[A]), &sequence);
        assert_eq!(outcome.boundaries, 0);
        assert_eq!(pieces(&sequence, &fragments), vec![sequence.clone()]);
    }

    #[test]
    fn indel_aware_matching_cuts_at_the_bases_the_alignment_consumed() {
        let mut stage = stage(&[A]);
        stage.params.allow_indels = true;
        // The planted copy is missing the adapter's fifth base, so the delimiter
        // occupies nine read bases while aligning ten adapter bases. Cutting by
        // the aligned length rather than the consumed length would leave one base
        // of adapter at the head of the suffix.
        let mut damaged = A.to_vec();
        damaged.remove(4);
        let sequence = [b"GGGGG".as_slice(), &damaged, b"TTTTT"].concat();
        let (fragments, outcome) = run(&stage, &sequence);
        assert_eq!(outcome.boundaries, 1, "a deletion is within the budget");
        assert_eq!(
            pieces(&sequence, &fragments),
            vec![b"GGGGG".to_vec(), b"TTTTT".to_vec()]
        );
    }

    #[test]
    fn terminal_fragments_can_be_discarded() {
        let sequence = [b"GGGGG".as_slice(), A, b"TTTTTTT", B, b"CCCCC"].concat();
        let mut stage = stage(&[A, B]);
        stage.terminal = Terminal::Discard;
        let (fragments, outcome) = run(&stage, &sequence);
        assert_eq!(
            pieces(&sequence, &fragments),
            vec![b"TTTTTTT".to_vec()],
            "only the doubly-delimited middle survives"
        );
        assert_eq!(outcome.terminal_discarded, 2);
        assert_eq!(outcome.internal, 1);
        assert_eq!(
            fragments[0].index, 0,
            "indices are assigned after filtering"
        );
    }

    #[test]
    fn the_minimum_segment_length_is_enforced() {
        let sequence = [b"GG".as_slice(), A, b"TTTTTTTT"].concat();
        let mut stage = stage(&[A]);
        stage.min_segment_length = 3;
        let (fragments, outcome) = run(&stage, &sequence);
        assert_eq!(pieces(&sequence, &fragments), vec![b"TTTTTTTT".to_vec()]);
        assert_eq!(outcome.too_short, 1);
        // Exactly at the minimum is kept.
        stage.min_segment_length = 2;
        let (fragments, _) = run(&stage, &sequence);
        assert_eq!(fragments.len(), 2);
    }

    #[test]
    fn the_segment_limit_bounds_a_pathological_read() {
        // Twelve delimiters would give thirteen fragments; the limit caps it.
        let mut sequence = Vec::new();
        for _ in 0..12 {
            sequence.extend_from_slice(b"GGGG");
            sequence.extend_from_slice(A);
        }
        sequence.extend_from_slice(b"CCCC");
        let mut stage = stage(&[A]);
        stage.max_segments = 5;
        let (fragments, outcome) = run(&stage, &sequence);
        assert_eq!(fragments.len(), 5);
        assert_eq!(outcome.over_limit, 8);
        assert!(fragments.iter().enumerate().all(|(i, f)| f.index == i));
    }

    #[test]
    fn fragments_never_overlap_and_exclude_the_adapters() {
        let sequence = [b"GGGGG".as_slice(), A, b"TTTTTTT", B, b"CCCCC"].concat();
        let (fragments, _) = run(&stage(&[A, B]), &sequence);
        for window in fragments.windows(2) {
            assert!(
                window[0].span.end <= window[1].span.start,
                "fragments must be ordered and disjoint"
            );
        }
        for fragment in &fragments {
            let piece = &sequence[fragment.span.start..fragment.span.end];
            assert!(!piece.windows(A.len()).any(|w| w == A));
            assert!(!piece.windows(B.len()).any(|w| w == B));
        }
    }

    #[test]
    fn header_suffixes_are_formatted_in_one_place() {
        let mut out = Vec::new();
        fragment_header(
            &mut out,
            b"SRR1.42 1:N:0",
            Fragment {
                span: Span { start: 10, end: 40 },
                index: 2,
                left_adapter: Some(0),
                right_adapter: None,
            },
        );
        assert_eq!(out, b"SRR1.42 1:N:0|segment=2|span=10-40");
    }

    #[test]
    fn validation_rejects_unusable_configurations() {
        assert!(
            SegmentStage::new(Vec::new(), params(), Terminal::Keep, 1, 8).is_err(),
            "no adapters"
        );
        let adapters = vec![Adapter::new("a", A).unwrap()];
        assert!(
            SegmentStage::new(adapters.clone(), params(), Terminal::Keep, 1, 0).is_err(),
            "zero segment limit"
        );
        let short = vec![Adapter::new("a", b"ACGT").unwrap()];
        assert!(SegmentStage::new(short, params(), Terminal::Keep, 1, 8).is_err());
    }

    #[test]
    fn find_all_reuses_its_buffer() {
        let adapters = vec![Adapter::new("a", A).unwrap()];
        let sequence = [b"GG".as_slice(), A, b"TT"].concat();
        let mut output = Vec::new();
        let mut dp = Vec::new();
        find_all(&adapters, params(), &sequence, &mut dp, &mut output);
        let first = output.len();
        assert!(first >= 1);
        // A second call must not accumulate onto the first.
        find_all(&adapters, params(), &sequence, &mut dp, &mut output);
        assert_eq!(output.len(), first);
        // A read too short for any match empties the buffer.
        find_all(&adapters, params(), b"AC", &mut dp, &mut output);
        assert!(output.is_empty());
    }
}
