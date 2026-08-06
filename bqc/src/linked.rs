// Copyright (c) 2026 The Hiller Lab at the Senckenberg Gesellschaft für Naturforschung
// Distributed under the terms of the GNU General Public License, Version 3.0.

//! Adapter-linked segmentation.
//!
//! An amplicon read looks like this:
//!
//! ```text
//! NNNNN [5' ADAPTER] INSERT [3' ADAPTER] NNNNN
//! ```
//!
//! A linked adapter definition names both flanks, and the stage retains what lies
//! between them. It stays a one-record-in, one-record-out transformation: only the
//! retained span changes. Emitting several fragments from one read is the separate
//! `segment` command.
//!
//! # Matching
//!
//! For each definition, in order:
//!
//! 1. Find the 5' adapter within `max_five_prime_offset` bases of the read start.
//! 2. Find the 3' adapter downstream of the 5' boundary.
//! 3. Reject pairs whose adapters occur in the wrong order.
//! 4. Take the span between the two adapter boundaries.
//! 5. Enforce `minimum_insert_length`.
//! 6. Apply the `both` or `either` requirement.
//!
//! Both searches reuse the ordinary adapter matcher, so the error budget, the
//! `N`-never-matches rule and indel awareness all behave exactly as they do for a
//! plain 3' adapter.
//!
//! # Selection
//!
//! Among the definitions that produced a valid candidate, the winner is chosen by,
//! in order:
//!
//! 1. Satisfying both adapters, over satisfying only one.
//! 2. Fewer total edit operations.
//! 3. Lower combined error rate.
//! 4. Greater total adapter overlap.
//! 5. Longer retained insert.
//! 6. Declaration order.

use serde::{Deserialize, Serialize};

use crate::adapter::{Adapter, AdapterHit, AdapterParams, find_five_prime, find_three_prime};
use crate::error::{Error, Result};
use crate::read::Span;

/// Default number of leading bases the 5' adapter may start within.
pub const DEFAULT_FIVE_PRIME_OFFSET: usize = 3;

/// How many of a definition's two adapters must match.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize, clap::ValueEnum)]
#[serde(rename_all = "lowercase")]
pub enum Require {
    /// Both adapters must match (default).
    #[default]
    Both,
    /// Either adapter alone is enough.
    Either,
}

/// What happens to a read no definition matched.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize, clap::ValueEnum)]
#[serde(rename_all = "kebab-case")]
pub enum Unmatched {
    /// Continue through the ordinary adapter options (default).
    #[default]
    Continue,
    /// Leave the read exactly as it is.
    Keep,
    /// Mark the read `LINKED_UNMATCHED` so filtering rejects it.
    Fail,
}

/// One linked adapter definition.
#[derive(Debug, Clone, Serialize)]
pub struct LinkedAdapter {
    pub name: String,
    pub five_prime: Adapter,
    pub three_prime: Adapter,
    pub require: Require,
    /// Bases the 5' adapter may start within.
    pub max_five_prime_offset: usize,
    /// Bases of 3' adapter that may hang past the read's end. `None` leaves the
    /// matcher's own partial-match rule in charge.
    pub max_three_prime_overhang: Option<usize>,
    pub minimum_insert_length: usize,
}

/// The compiled stage: definitions per mate plus the shared policy.
#[derive(Debug, Clone, Serialize)]
pub struct LinkedStage {
    pub r1: Vec<LinkedAdapter>,
    pub r2: Vec<LinkedAdapter>,
    pub unmatched: Unmatched,
}

impl LinkedStage {
    /// Builds the stage, rejecting definitions that can never match.
    pub fn new(
        r1: Vec<LinkedAdapter>,
        r2: Vec<LinkedAdapter>,
        unmatched: Unmatched,
        params: AdapterParams,
    ) -> Result<Self> {
        if r1.is_empty() && r2.is_empty() {
            return Err(Error::config(
                "linked segmentation requires --linked-5p-r1 and --linked-3p-r1 \
                 (or the R2 equivalents)",
            ));
        }
        for definition in r1.iter().chain(&r2) {
            for adapter in [&definition.five_prime, &definition.three_prime] {
                if adapter.len() < params.min_overlap {
                    return Err(Error::InvalidAdapter(format!(
                        "linked adapter '{}' side '{}' is {} bases long, shorter than \
                         --min-overlap {}",
                        definition.name,
                        adapter.name,
                        adapter.len(),
                        params.min_overlap
                    )));
                }
            }
        }
        Ok(Self { r1, r2, unmatched })
    }

    /// Definitions declared for `mate`.
    #[must_use]
    pub fn definitions(&self, mate: crate::process::Mate) -> &[LinkedAdapter] {
        match mate {
            crate::process::Mate::R1 => &self.r1,
            crate::process::Mate::R2 => &self.r2,
        }
    }
}

/// Why a linked candidate was not usable.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Rejected {
    /// Neither adapter matched.
    Neither,
    /// Only the 5' adapter matched, and `both` was required.
    FiveOnly,
    /// Only the 3' adapter matched, and `both` was required.
    ThreeOnly,
    /// The 3' adapter came before the 5' adapter.
    OutOfOrder,
    /// The insert was shorter than the minimum.
    InsertTooShort,
}

/// One accepted linked match.
#[derive(Debug, Clone, Copy)]
pub struct LinkedMatch {
    /// Index of the winning definition within its mate's list.
    pub definition: usize,
    pub five_prime: Option<AdapterHit>,
    pub three_prime: Option<AdapterHit>,
    /// The insert, in coordinates of the span that was searched.
    pub retained: Span,
}

impl LinkedMatch {
    /// Total edit operations across both flanks.
    #[must_use]
    pub fn errors(&self) -> usize {
        self.five_prime.map_or(0, |hit| hit.errors) + self.three_prime.map_or(0, |hit| hit.errors)
    }

    /// Total adapter bases aligned across both flanks.
    #[must_use]
    pub fn overlap(&self) -> usize {
        self.five_prime.map_or(0, |hit| hit.overlap) + self.three_prime.map_or(0, |hit| hit.overlap)
    }

    /// Whether both flanks matched.
    #[must_use]
    pub fn both(&self) -> bool {
        self.five_prime.is_some() && self.three_prime.is_some()
    }

    /// Combined error rate over the matched adapter bases.
    #[must_use]
    pub fn error_rate(&self) -> f64 {
        let overlap = self.overlap();
        if overlap == 0 {
            return 0.0;
        }
        self.errors() as f64 / overlap as f64
    }

    /// Bases removed ahead of the insert.
    #[must_use]
    pub fn leading_removed(&self, searched: Span) -> usize {
        self.retained.start - searched.start
    }

    /// Bases removed after the insert.
    #[must_use]
    pub fn trailing_removed(&self, searched: Span) -> usize {
        searched.end - self.retained.end
    }
}

/// The outcome of the linked stage for one mate.
#[derive(Debug, Clone, Copy)]
pub enum LinkedOutcome {
    /// A definition matched; the span was narrowed to the insert.
    Matched(LinkedMatch),
    /// Nothing matched. Carries the most informative reason, for statistics.
    Unmatched(Rejected),
}

/// Finds the winning linked match in `sequence`, narrowing `span` on success.
///
/// `sequence` is the span's own slice, so returned coordinates are relative to it.
#[must_use]
pub fn find(
    definitions: &[LinkedAdapter],
    params: AdapterParams,
    sequence: &[u8],
) -> LinkedOutcome {
    let mut best: Option<LinkedMatch> = None;
    // The least dismissive reason seen, so statistics report the most specific
    // thing that went wrong rather than "neither".
    let mut reason = Rejected::Neither;
    for (index, definition) in definitions.iter().enumerate() {
        match evaluate(definition, params, sequence, index) {
            Ok(candidate) => {
                let improves = best.is_none_or(|current| better(&candidate, &current));
                if improves {
                    best = Some(candidate);
                }
            }
            Err(rejected) => {
                if severity(rejected) > severity(reason) {
                    reason = rejected;
                }
            }
        }
    }
    match best {
        Some(candidate) => LinkedOutcome::Matched(candidate),
        None => LinkedOutcome::Unmatched(reason),
    }
}

/// How specific a rejection reason is; the most specific one is reported.
fn severity(reason: Rejected) -> u8 {
    match reason {
        Rejected::Neither => 0,
        Rejected::FiveOnly | Rejected::ThreeOnly => 1,
        Rejected::OutOfOrder => 2,
        Rejected::InsertTooShort => 3,
    }
}

/// Whether `candidate` beats `current` under the documented ordering.
fn better(candidate: &LinkedMatch, current: &LinkedMatch) -> bool {
    let key = |m: &LinkedMatch| {
        (
            !m.both(),                           // both-adapter candidates first
            m.errors(),                          // fewer edits
            (m.error_rate() * 1e9) as u64,       // lower combined error rate
            std::cmp::Reverse(m.overlap()),      // greater adapter overlap
            std::cmp::Reverse(m.retained.len()), // longer insert
        )
    };
    key(candidate) < key(current)
}

/// Evaluates one definition against `sequence`.
fn evaluate(
    definition: &LinkedAdapter,
    params: AdapterParams,
    sequence: &[u8],
    index: usize,
) -> std::result::Result<LinkedMatch, Rejected> {
    let five = find_five_prime(
        std::slice::from_ref(&definition.five_prime),
        params,
        sequence,
        definition.max_five_prime_offset,
    );
    // The insert can only start after the 5' adapter, so the 3' search begins
    // there — which also avoids re-scanning what the 5' match already claimed.
    // `consumed`, not `overlap`: an indel makes the read bases the flank occupies
    // differ from the adapter bases it aligned.
    let insert_start = five.map_or(0, |hit| hit.start + hit.consumed);
    if insert_start > sequence.len() {
        return Err(Rejected::OutOfOrder);
    }
    let three = find_three_prime(
        std::slice::from_ref(&definition.three_prime),
        params,
        &sequence[insert_start..],
    )
    .map(|mut hit| {
        hit.start += insert_start; // back into the searched span's coordinates
        hit
    });

    if let (Some(three), Some(overhang)) = (three, definition.max_three_prime_overhang) {
        // A 3' adapter that runs off the read end is only accepted if the part
        // hanging over is within the configured limit.
        let visible = sequence.len() - three.start;
        if definition.three_prime.len() > visible
            && definition.three_prime.len() - visible > overhang
        {
            return Err(Rejected::ThreeOnly);
        }
    }

    match (five, three, definition.require) {
        (None, None, _) => return Err(Rejected::Neither),
        (Some(_), None, Require::Both) => return Err(Rejected::FiveOnly),
        (None, Some(_), Require::Both) => return Err(Rejected::ThreeOnly),
        _ => {}
    }

    let end = three.map_or(sequence.len(), |hit| hit.start);
    if end < insert_start {
        return Err(Rejected::OutOfOrder);
    }
    let retained = Span {
        start: insert_start,
        end,
    };
    if retained.len() < definition.minimum_insert_length {
        return Err(Rejected::InsertTooShort);
    }
    Ok(LinkedMatch {
        definition: index,
        five_prime: five,
        three_prime: three,
        retained,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    const FIVE: &[u8] = b"ACGTACGTACGT";
    const THREE: &[u8] = b"TGCATGCATGCA";
    const INSERT: &[u8] = b"GGACTTGACCTTAAGGCATTCCAAGT";

    fn params() -> AdapterParams {
        AdapterParams {
            min_overlap: 8,
            max_error_rate: 0.15,
            max_errors: None,
            allow_indels: false,
        }
    }

    fn definition(name: &str) -> LinkedAdapter {
        LinkedAdapter {
            name: name.to_string(),
            five_prime: Adapter::new("5p", FIVE).unwrap(),
            three_prime: Adapter::new("3p", THREE).unwrap(),
            require: Require::Both,
            max_five_prime_offset: DEFAULT_FIVE_PRIME_OFFSET,
            max_three_prime_overhang: None,
            minimum_insert_length: 5,
        }
    }

    /// `prefix + 5' + insert + 3' + suffix`
    fn read(prefix: &[u8], insert: &[u8], suffix: &[u8]) -> Vec<u8> {
        [prefix, FIVE, insert, THREE, suffix].concat()
    }

    fn find_one(definition: &LinkedAdapter, sequence: &[u8]) -> LinkedOutcome {
        find(std::slice::from_ref(definition), params(), sequence)
    }

    fn matched(outcome: LinkedOutcome) -> LinkedMatch {
        match outcome {
            LinkedOutcome::Matched(m) => m,
            LinkedOutcome::Unmatched(reason) => panic!("expected a match, got {reason:?}"),
        }
    }

    fn rejected(outcome: LinkedOutcome) -> Rejected {
        match outcome {
            LinkedOutcome::Matched(m) => panic!("expected no match, got {m:?}"),
            LinkedOutcome::Unmatched(reason) => reason,
        }
    }

    #[test]
    fn exact_linked_adapters_retain_the_insert() {
        let sequence = read(b"", INSERT, b"");
        let outcome = matched(find_one(&definition("amplicon"), &sequence));
        assert!(outcome.both());
        assert_eq!(outcome.errors(), 0);
        assert_eq!(
            &sequence[outcome.retained.start..outcome.retained.end],
            INSERT
        );
        assert_eq!(outcome.retained.len(), INSERT.len());
    }

    #[test]
    fn flanking_bases_are_removed_from_both_ends() {
        let sequence = read(b"AAT", INSERT, b"CCCCC");
        let searched = Span::full(sequence.len());
        let outcome = matched(find_one(&definition("amplicon"), &sequence));
        assert_eq!(
            &sequence[outcome.retained.start..outcome.retained.end],
            INSERT
        );
        assert_eq!(outcome.leading_removed(searched), 3 + FIVE.len());
        assert_eq!(outcome.trailing_removed(searched), THREE.len() + 5);
    }

    #[test]
    fn the_five_prime_offset_limit_is_enforced() {
        let within = read(b"AAT", INSERT, b"");
        assert!(matches!(
            find_one(&definition("amplicon"), &within),
            LinkedOutcome::Matched(_)
        ));
        // One base too far: the 5' adapter is no longer found, so `both` fails.
        let beyond = read(b"AATG", INSERT, b"");
        assert_eq!(
            rejected(find_one(&definition("amplicon"), &beyond)),
            Rejected::ThreeOnly
        );
    }

    #[test]
    fn substitutions_inside_the_budget_still_match() {
        let mut sequence = read(b"", INSERT, b"");
        sequence[2] = if sequence[2] == b'A' { b'C' } else { b'A' };
        let outcome = matched(find_one(&definition("amplicon"), &sequence));
        assert_eq!(outcome.errors(), 1);
        assert!(outcome.both());
    }

    #[test]
    fn indels_are_matched_when_enabled() {
        let mut flank = FIVE.to_vec();
        flank.insert(5, b'G');
        let sequence = [flank.as_slice(), INSERT, THREE].concat();
        let indel_params = AdapterParams {
            allow_indels: true,
            ..params()
        };
        let outcome = find(
            std::slice::from_ref(&definition("amplicon")),
            indel_params,
            &sequence,
        );
        let outcome = matched(outcome);
        assert!(outcome.both(), "an inserted base must not break the flank");
        assert_eq!(
            &sequence[outcome.retained.start..outcome.retained.end],
            INSERT
        );
    }

    #[test]
    fn a_missing_five_prime_adapter_is_reported() {
        let sequence = [INSERT, THREE].concat();
        assert_eq!(
            rejected(find_one(&definition("amplicon"), &sequence)),
            Rejected::ThreeOnly
        );
    }

    #[test]
    fn a_missing_three_prime_adapter_is_reported() {
        let sequence = [FIVE, INSERT].concat();
        assert_eq!(
            rejected(find_one(&definition("amplicon"), &sequence)),
            Rejected::FiveOnly
        );
    }

    #[test]
    fn neither_adapter_present_is_reported() {
        assert_eq!(
            rejected(find_one(&definition("amplicon"), INSERT)),
            Rejected::Neither
        );
    }

    #[test]
    fn either_accepts_one_flank_alone() {
        let mut definition = definition("amplicon");
        definition.require = Require::Either;
        // Only the 5' flank: the insert runs to the read end.
        let sequence = [FIVE, INSERT].concat();
        let outcome = matched(find_one(&definition, &sequence));
        assert!(!outcome.both());
        assert_eq!(
            &sequence[outcome.retained.start..outcome.retained.end],
            INSERT
        );
        // Only the 3' flank: the insert starts at the read start.
        let sequence = [INSERT, THREE].concat();
        let outcome = matched(find_one(&definition, &sequence));
        assert_eq!(outcome.retained.start, 0);
        assert_eq!(
            &sequence[outcome.retained.start..outcome.retained.end],
            INSERT
        );
    }

    #[test]
    fn reversed_adapter_order_is_rejected() {
        // The 3' adapter appears first, and the 5' adapter after it. Searching
        // downstream of the 5' match finds no 3' adapter, so `both` fails.
        let sequence = [THREE, INSERT, FIVE].concat();
        let outcome = find_one(&definition("amplicon"), &sequence);
        assert!(
            matches!(outcome, LinkedOutcome::Unmatched(_)),
            "{outcome:?}"
        );
    }

    #[test]
    fn an_empty_insert_is_rejected_by_the_minimum() {
        let sequence = [FIVE, THREE].concat();
        let mut definition = definition("amplicon");
        definition.minimum_insert_length = 1;
        assert_eq!(
            rejected(find_one(&definition, &sequence)),
            Rejected::InsertTooShort
        );
    }

    #[test]
    fn an_insert_exactly_at_the_minimum_is_kept() {
        let insert = &INSERT[..10];
        let sequence = [FIVE, insert, THREE].concat();
        let mut definition = definition("amplicon");
        definition.minimum_insert_length = 10;
        let outcome = matched(find_one(&definition, &sequence));
        assert_eq!(outcome.retained.len(), 10);
        definition.minimum_insert_length = 11;
        assert_eq!(
            rejected(find_one(&definition, &sequence)),
            Rejected::InsertTooShort
        );
    }

    #[test]
    fn both_adapters_beat_one_and_ties_fall_to_declaration_order() {
        // Definition 0 matches only its 5' flank; definition 1 matches both.
        let other = LinkedAdapter {
            name: "other".to_string(),
            three_prime: Adapter::new("3p-other", b"GGGGGGGGGGGG").unwrap(),
            require: Require::Either,
            ..definition("other")
        };
        let sequence = read(b"", INSERT, b"");
        let outcome = matched(find(&[other, definition("amplicon")], params(), &sequence));
        assert_eq!(outcome.definition, 1, "the both-adapter candidate wins");
        assert!(outcome.both());

        // Two identical definitions: the first declared wins.
        let outcome = matched(find(
            &[definition("first"), definition("second")],
            params(),
            &sequence,
        ));
        assert_eq!(outcome.definition, 0);
    }

    #[test]
    fn fewer_errors_wins_between_two_matching_definitions() {
        let mut noisy = definition("noisy");
        let mut flank = FIVE.to_vec();
        flank[1] = if flank[1] == b'A' { b'C' } else { b'A' };
        noisy.five_prime = Adapter::new("5p-noisy", &flank).unwrap();
        let sequence = read(b"", INSERT, b"");
        // `noisy` matches with one error, the clean definition with none.
        let outcome = matched(find(&[noisy, definition("clean")], params(), &sequence));
        assert_eq!(outcome.definition, 1);
        assert_eq!(outcome.errors(), 0);
    }

    #[test]
    fn the_three_prime_overhang_limit_is_enforced() {
        // Only four bases of the 3' adapter are visible at the read end.
        let sequence = [FIVE, INSERT, &THREE[..4]].concat();
        let mut definition = definition("amplicon");
        definition.max_five_prime_offset = 0;
        // The default matcher needs `min_overlap` bases, so four never match.
        assert_eq!(
            rejected(find_one(&definition, &sequence)),
            Rejected::FiveOnly
        );
        // With a shorter minimum overlap the partial flank matches, and the
        // overhang limit then decides.
        let lenient = AdapterParams {
            min_overlap: 4,
            ..params()
        };
        definition.max_three_prime_overhang = Some(8);
        let outcome = find(std::slice::from_ref(&definition), lenient, &sequence);
        assert!(matches!(outcome, LinkedOutcome::Matched(_)), "{outcome:?}");
        definition.max_three_prime_overhang = Some(7);
        assert_eq!(
            rejected(find(std::slice::from_ref(&definition), lenient, &sequence)),
            Rejected::ThreeOnly
        );
    }

    #[test]
    fn stage_validation_rejects_unusable_definitions() {
        assert!(LinkedStage::new(Vec::new(), Vec::new(), Unmatched::Continue, params()).is_err());
        let mut short = definition("short");
        short.five_prime = Adapter::new("5p", b"ACGT").unwrap();
        let error =
            LinkedStage::new(vec![short], Vec::new(), Unmatched::Continue, params()).unwrap_err();
        assert!(
            format!("{error}").contains("shorter than --min-overlap"),
            "{error}"
        );
    }
}
