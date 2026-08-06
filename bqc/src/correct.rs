// Copyright (c) 2026 The Hiller Lab at the Senckenberg Gesellschaft für Naturforschung
// Distributed under the terms of the GNU General Public License, Version 3.0.

//! Paired-overlap base correction.
//!
//! Where R1 and R2 genuinely overlap they sequence the same molecule twice, from
//! opposite ends. A disagreement between them is therefore a sequencing error in
//! one of the two reads, and when one side is confident and the other is not, the
//! confident base can replace the doubtful one.
//!
//! # Contract
//!
//! For every mismatch inside an accepted overlap:
//!
//! ```text
//! R1 >= donor quality and R2 <= recipient quality  -> correct R2 from R1
//! R2 >= donor quality and R1 <= recipient quality  -> correct R1 from R2
//! otherwise                                        -> leave both unchanged
//! ```
//!
//! Thresholds are inclusive, and validation requires `donor > recipient`, which
//! makes the two rules mutually exclusive: no mismatch can be corrected in both
//! directions. Defaults are Q30 and Q14.
//!
//! * The donor base is written into the recipient's own orientation, so it is
//!   complemented on the way across.
//! * The donor's raw quality byte is copied exactly. Qualities are never raised
//!   on bases that already agree.
//! * A low-quality `N` may be corrected, but an `N` — or any non-canonical base —
//!   is never used as donor evidence.
//! * Only ungapped overlaps are corrected. `bqc`'s overlap inference is
//!   ungapped by construction, so every accepted overlap qualifies.
//!
//! # Zero-copy
//!
//! Every edit is planned against the original pair before anything is mutated,
//! and only a corrected mate is copied: its retained sequence and quality are
//! written into reused worker buffers, and the other mate stays borrowed from the
//! memory-mapped block. Pairs with no corrections copy nothing at all.

use serde::{Deserialize, Serialize};

use crate::error::{Error, Result};
use crate::overlap::{Overlap, bases_disagree, complement};
use crate::process::Mate;
use crate::read::{ReadView, phred};

/// How much detail the correction log records.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize, clap::ValueEnum)]
#[serde(rename_all = "lowercase")]
pub enum LogDetail {
    /// One row per corrected pair.
    #[default]
    Reads,
    /// One row per corrected base.
    Bases,
}

/// Compiled correction configuration.
#[derive(Debug, Clone, Copy, Serialize)]
pub struct CorrectionStage {
    /// Lowest Phred score accepted as donor evidence, inclusive.
    pub donor_quality: u8,
    /// Highest Phred score that may be overwritten, inclusive.
    pub recipient_quality: u8,
    pub log_detail: LogDetail,
}

impl Default for CorrectionStage {
    fn default() -> Self {
        Self {
            donor_quality: 30,
            recipient_quality: 14,
            log_detail: LogDetail::Reads,
        }
    }
}

impl CorrectionStage {
    /// Validates the thresholds.
    ///
    /// `donor > recipient` is required, which is what keeps the two correction
    /// rules mutually exclusive.
    pub fn validate(self) -> Result<Self> {
        if self.donor_quality <= self.recipient_quality {
            return Err(Error::config(format!(
                "--donor-quality ({}) must be greater than --recipient-quality ({}); \
                 otherwise a base could be both donor and recipient",
                self.donor_quality, self.recipient_quality
            )));
        }
        Ok(self)
    }
}

/// One planned base correction, in original mate coordinates.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CorrectionEdit {
    /// Mate being corrected.
    pub target: Mate,
    /// Position corrected, in the target mate's own coordinates.
    pub target_position: usize,
    /// Position the evidence came from, in the donor mate's own coordinates.
    pub donor_position: usize,
    pub original_base: u8,
    pub corrected_base: u8,
    /// Raw quality bytes, not decoded Phred scores.
    pub original_quality: u8,
    pub corrected_quality: u8,
}

impl CorrectionEdit {
    /// The mate the evidence came from.
    #[must_use]
    pub fn donor(&self) -> Mate {
        match self.target {
            Mate::R1 => Mate::R2,
            Mate::R2 => Mate::R1,
        }
    }
}

/// What correcting one pair produced.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct CorrectionSummary {
    /// Mismatches inside the accepted overlap.
    pub mismatches: u32,
    pub corrected_r1: u32,
    pub corrected_r2: u32,
    /// Mismatches the quality rule left alone.
    pub unresolved: u32,
    /// Mismatches the quality rule would have corrected, had the donor base been
    /// canonical.
    pub skipped_noncanonical: u32,
}

impl CorrectionSummary {
    /// Total bases corrected in this pair.
    #[must_use]
    pub fn corrected(&self) -> u32 {
        self.corrected_r1 + self.corrected_r2
    }

    /// Whether both mates were corrected.
    #[must_use]
    pub fn both_mates(&self) -> bool {
        self.corrected_r1 > 0 && self.corrected_r2 > 0
    }

    /// Bases corrected in `mate`.
    #[must_use]
    pub fn corrected_in(&self, mate: Mate) -> u32 {
        match mate {
            Mate::R1 => self.corrected_r1,
            Mate::R2 => self.corrected_r2,
        }
    }
}

/// Reused worker buffers. Correction never allocates in the record loop.
#[derive(Debug, Default)]
pub struct CorrectionScratch {
    pub edits: Vec<CorrectionEdit>,
    r1_sequence: Vec<u8>,
    r1_quality: Vec<u8>,
    r2_sequence: Vec<u8>,
    r2_quality: Vec<u8>,
}

impl CorrectionScratch {
    /// Clears the planned edits, keeping the allocations.
    pub fn clear(&mut self) {
        self.edits.clear();
    }

    /// The corrected bytes for `mate`, if it was corrected.
    ///
    /// Valid until the next pair is planned. The writer reads the corrected
    /// sequence and quality from here rather than from the CBQ record.
    #[must_use]
    pub fn corrected(&self, mate: Mate) -> Option<(&[u8], &[u8])> {
        if !self.edits.iter().any(|edit| edit.target == mate) {
            return None;
        }
        Some(match mate {
            Mate::R1 => (&self.r1_sequence, &self.r1_quality),
            Mate::R2 => (&self.r2_sequence, &self.r2_quality),
        })
    }

    /// Applies the planned edits to both mates, returning views of the results.
    ///
    /// A mate with no edits comes back borrowed from the original record, so
    /// uncorrected reads — the overwhelming majority — copy nothing. Both mates
    /// are produced by one call because they live in separate buffers: that is
    /// what lets two views coexist under a single mutable borrow of the scratch.
    pub fn apply_pair<'scratch, 'read: 'scratch>(
        &'scratch mut self,
        r1: ReadView<'read>,
        r2: ReadView<'read>,
    ) -> (ReadView<'scratch>, ReadView<'scratch>) {
        let Self {
            edits,
            r1_sequence,
            r1_quality,
            r2_sequence,
            r2_quality,
        } = self;
        (
            apply_into(edits, Mate::R1, r1, r1_sequence, r1_quality),
            apply_into(edits, Mate::R2, r2, r2_sequence, r2_quality),
        )
    }
}

/// Writes `read` with `mate`'s edits applied into the given buffers.
///
/// Returns the original view untouched when no edit targets `mate`.
fn apply_into<'a>(
    edits: &[CorrectionEdit],
    mate: Mate,
    read: ReadView<'a>,
    sequence: &'a mut Vec<u8>,
    quality: &'a mut Vec<u8>,
) -> ReadView<'a> {
    if !edits.iter().any(|edit| edit.target == mate) {
        return read;
    }
    sequence.clear();
    sequence.extend_from_slice(read.sequence);
    quality.clear();
    if let Some(original) = read.quality {
        quality.extend_from_slice(original);
    }
    for edit in edits.iter().filter(|edit| edit.target == mate) {
        sequence[edit.target_position] = edit.corrected_base;
        if !quality.is_empty() {
            quality[edit.target_position] = edit.corrected_quality;
        }
    }
    ReadView::unchecked(
        sequence,
        read.quality.is_some().then_some(quality.as_slice()),
    )
}

/// Whether a base can act as donor evidence.
#[inline]
fn is_canonical(base: u8) -> bool {
    matches!(base, b'A' | b'C' | b'G' | b'T')
}

/// Plans the corrections for one pair, appending them to `edits`.
///
/// Nothing is mutated: every decision is made against the original pair. Edits
/// come back ordered by mate and then position, which is the order the log
/// requires.
pub fn plan(
    stage: CorrectionStage,
    r1: ReadView<'_>,
    r2: ReadView<'_>,
    overlap: Overlap,
    edits: &mut Vec<CorrectionEdit>,
) -> CorrectionSummary {
    let mut summary = CorrectionSummary::default();
    let (Some(quality1), Some(quality2)) = (r1.quality, r2.quality) else {
        // Correction requires stored qualities; configuration rejects the run
        // before this point, so this is only a guard.
        return summary;
    };

    for (position1, position2) in overlap.aligned_positions(r2.sequence.len()) {
        let base1 = r1.sequence[position1];
        let base2 = r2.sequence[position2];
        if !bases_disagree(base1, base2) {
            continue;
        }
        summary.mismatches += 1;

        let phred1 = phred(quality1[position1]);
        let phred2 = phred(quality2[position2]);
        let r1_donates = phred1 >= stage.donor_quality && phred2 <= stage.recipient_quality;
        let r2_donates = phred2 >= stage.donor_quality && phred1 <= stage.recipient_quality;
        debug_assert!(
            !(r1_donates && r2_donates),
            "donor > recipient makes the rules exclusive"
        );

        if r1_donates && is_canonical(base1) {
            edits.push(CorrectionEdit {
                target: Mate::R2,
                target_position: position2,
                donor_position: position1,
                original_base: base2,
                // The mates are in opposite orientations.
                corrected_base: complement(base1),
                original_quality: quality2[position2],
                corrected_quality: quality1[position1],
            });
            summary.corrected_r2 += 1;
        } else if r2_donates && is_canonical(base2) {
            edits.push(CorrectionEdit {
                target: Mate::R1,
                target_position: position1,
                donor_position: position2,
                original_base: base1,
                corrected_base: complement(base2),
                original_quality: quality1[position1],
                corrected_quality: quality2[position2],
            });
            summary.corrected_r1 += 1;
        } else if r1_donates || r2_donates {
            // The qualities permitted a correction but the donor base was not
            // canonical, which is a different fact from "the qualities said no".
            summary.skipped_noncanonical += 1;
        } else {
            summary.unresolved += 1;
        }
    }

    // The log is ordered by mate then position. R1 targets accumulate in
    // ascending order and R2 targets in descending order, so a sort is the
    // clearest way to state the requirement.
    edits.sort_unstable_by_key(|edit| (edit.target.index(), edit.target_position));
    summary
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::overlap::{OverlapParams, find_overlap, tests::revcomp};

    /// A quality string from Phred scores.
    fn quality(scores: &[u8]) -> Vec<u8> {
        scores.iter().map(|&score| b'!' + score).collect()
    }

    fn stage() -> CorrectionStage {
        CorrectionStage::default()
    }

    /// A 40 base insert sequenced from both ends with a full overlap.
    fn pair() -> (Vec<u8>, Vec<u8>) {
        let insert = b"ACGTTGCAAGGCCTTAACGTACGTTGCAAGGCCTTAACGT".to_vec();
        let mate2 = revcomp(&insert);
        (insert, mate2)
    }

    /// Plans corrections for a pair whose overlap is found by the real matcher.
    fn plan_pair(
        stage: CorrectionStage,
        sequence1: &[u8],
        quality1: &[u8],
        sequence2: &[u8],
        quality2: &[u8],
    ) -> (Vec<CorrectionEdit>, CorrectionSummary) {
        let params = OverlapParams {
            min_overlap: 20,
            max_error_rate: 0.20,
        };
        let overlap = find_overlap(sequence1, sequence2, params).expect("overlap");
        let r1 = ReadView::unchecked(sequence1, Some(quality1));
        let r2 = ReadView::unchecked(sequence2, Some(quality2));
        let mut edits = Vec::new();
        let summary = plan(stage, r1, r2, overlap, &mut edits);
        (edits, summary)
    }

    #[test]
    fn validation_requires_a_donor_above_the_recipient() {
        assert!(stage().validate().is_ok());
        assert!(
            CorrectionStage {
                donor_quality: 20,
                recipient_quality: 20,
                ..stage()
            }
            .validate()
            .is_err()
        );
        assert!(
            CorrectionStage {
                donor_quality: 10,
                recipient_quality: 20,
                ..stage()
            }
            .validate()
            .is_err()
        );
    }

    #[test]
    fn a_perfect_overlap_needs_no_corrections() {
        let (mate1, mate2) = pair();
        let good = quality(&vec![35u8; mate1.len()]);
        let (edits, summary) = plan_pair(stage(), &mate1, &good, &mate2, &good);
        assert!(edits.is_empty());
        assert_eq!(summary, CorrectionSummary::default());
    }

    #[test]
    fn r1_corrects_r2_when_only_r1_is_confident() {
        let (mate1, mut mate2) = pair();
        // Break R2 at the position matching R1[10] and make it low quality.
        let position2 = mate2.len() - 1 - 10;
        let original = mate2[position2];
        mate2[position2] = if original == b'A' { b'C' } else { b'A' };
        let mut scores2 = vec![35u8; mate2.len()];
        scores2[position2] = 10;
        let quality1 = quality(&vec![35u8; mate1.len()]);
        let quality2 = quality(&scores2);

        let (edits, summary) = plan_pair(stage(), &mate1, &quality1, &mate2, &quality2);
        assert_eq!(edits.len(), 1);
        let edit = edits[0];
        assert_eq!(edit.target, Mate::R2);
        assert_eq!(edit.donor(), Mate::R1);
        assert_eq!(edit.target_position, position2);
        assert_eq!(edit.donor_position, 10);
        assert_eq!(edit.original_base, mate2[position2]);
        assert_eq!(edit.corrected_base, complement(mate1[10]));
        assert_eq!(edit.corrected_base, original, "the true base is restored");
        assert_eq!(edit.corrected_quality, quality1[10], "donor byte copied");
        assert_eq!(edit.original_quality, quality2[position2]);
        assert_eq!(summary.corrected_r2, 1);
        assert_eq!(summary.corrected_r1, 0);
        assert_eq!(summary.mismatches, 1);
    }

    #[test]
    fn r2_corrects_r1_in_the_other_direction() {
        let (mut mate1, mate2) = pair();
        let truth = mate1[7];
        mate1[7] = if truth == b'G' { b'T' } else { b'G' };
        let mut scores1 = vec![35u8; mate1.len()];
        scores1[7] = 5;
        let quality1 = quality(&scores1);
        let quality2 = quality(&vec![35u8; mate2.len()]);

        let (edits, summary) = plan_pair(stage(), &mate1, &quality1, &mate2, &quality2);
        assert_eq!(edits.len(), 1);
        assert_eq!(edits[0].target, Mate::R1);
        assert_eq!(edits[0].target_position, 7);
        assert_eq!(edits[0].corrected_base, truth);
        assert_eq!(summary.corrected_r1, 1);
    }

    #[test]
    fn thresholds_are_inclusive_at_both_ends() {
        let (mate1, mut mate2) = pair();
        let position2 = mate2.len() - 1 - 12;
        mate2[position2] = if mate2[position2] == b'A' { b'C' } else { b'A' };
        let mut scores2 = vec![35u8; mate2.len()];
        let mut scores1 = vec![35u8; mate1.len()];

        // Exactly Q30 donor and exactly Q14 recipient: corrected.
        scores1[12] = 30;
        scores2[position2] = 14;
        let (edits, _) = plan_pair(
            stage(),
            &mate1,
            &quality(&scores1),
            &mate2,
            &quality(&scores2),
        );
        assert_eq!(edits.len(), 1, "Q30 donor and Q14 recipient are inclusive");

        // Q29 donor: refused.
        scores1[12] = 29;
        scores2[position2] = 14;
        let (edits, summary) = plan_pair(
            stage(),
            &mate1,
            &quality(&scores1),
            &mate2,
            &quality(&scores2),
        );
        assert!(edits.is_empty(), "Q29 is below the donor threshold");
        assert_eq!(summary.unresolved, 1);

        // Q15 recipient: refused.
        scores1[12] = 30;
        scores2[position2] = 15;
        let (edits, summary) = plan_pair(
            stage(),
            &mate1,
            &quality(&scores1),
            &mate2,
            &quality(&scores2),
        );
        assert!(edits.is_empty(), "Q15 is above the recipient threshold");
        assert_eq!(summary.unresolved, 1);
    }

    #[test]
    fn both_high_and_both_low_are_left_alone() {
        let (mate1, mut mate2) = pair();
        let position2 = mate2.len() - 1 - 3;
        mate2[position2] = if mate2[position2] == b'A' { b'C' } else { b'A' };
        for score in [35u8, 5] {
            let scores = vec![score; mate1.len()];
            let (edits, summary) = plan_pair(
                stage(),
                &mate1,
                &quality(&scores),
                &mate2,
                &quality(&scores),
            );
            assert!(edits.is_empty(), "equal qualities never correct (Q{score})");
            assert_eq!(summary.unresolved, 1);
            assert_eq!(summary.mismatches, 1);
        }
    }

    #[test]
    fn a_low_quality_recipient_n_is_corrected() {
        let (mate1, mut mate2) = pair();
        let position2 = mate2.len() - 1 - 15;
        mate2[position2] = b'N';
        let mut scores2 = vec![35u8; mate2.len()];
        scores2[position2] = 2;
        let (edits, summary) = plan_pair(
            stage(),
            &mate1,
            &quality(&vec![35u8; mate1.len()]),
            &mate2,
            &quality(&scores2),
        );
        assert_eq!(edits.len(), 1);
        assert_eq!(edits[0].original_base, b'N');
        assert_eq!(edits[0].corrected_base, complement(mate1[15]));
        assert_eq!(summary.corrected_r2, 1);
    }

    #[test]
    fn an_n_is_never_donor_evidence() {
        let (mut mate1, mate2) = pair();
        // R1 carries a high-quality N; R2 is low quality at that position.
        mate1[9] = b'N';
        let mut scores2 = vec![35u8; mate2.len()];
        scores2[mate2.len() - 1 - 9] = 3;
        let (edits, summary) = plan_pair(
            stage(),
            &mate1,
            &quality(&vec![40u8; mate1.len()]),
            &mate2,
            &quality(&scores2),
        );
        assert!(edits.is_empty(), "a high-quality N must not donate");
        assert_eq!(summary.skipped_noncanonical, 1);
        assert_eq!(
            summary.unresolved, 0,
            "counted separately from a quality veto"
        );
    }

    #[test]
    fn multiple_edits_in_both_mates_are_ordered_by_mate_then_position() {
        let (mut mate1, mut mate2) = pair();
        let mut scores1 = vec![35u8; mate1.len()];
        let mut scores2 = vec![35u8; mate2.len()];
        // Two positions where R2 wins, two where R1 wins.
        for &position1 in &[4usize, 20] {
            mate1[position1] = if mate1[position1] == b'A' { b'C' } else { b'A' };
            scores1[position1] = 4;
        }
        for &position1 in &[11usize, 30] {
            let position2 = mate2.len() - 1 - position1;
            mate2[position2] = if mate2[position2] == b'A' { b'C' } else { b'A' };
            scores2[position2] = 4;
        }
        let (edits, summary) = plan_pair(
            stage(),
            &mate1,
            &quality(&scores1),
            &mate2,
            &quality(&scores2),
        );
        assert_eq!(summary.corrected_r1, 2);
        assert_eq!(summary.corrected_r2, 2);
        assert_eq!(summary.mismatches, 4);
        assert!(summary.both_mates());
        assert_eq!(summary.corrected(), 4);

        let keys: Vec<(usize, usize)> = edits
            .iter()
            .map(|edit| (edit.target.index(), edit.target_position))
            .collect();
        let mut sorted = keys.clone();
        sorted.sort_unstable();
        assert_eq!(keys, sorted, "log order is mate then position");
        assert_eq!(edits[0].target, Mate::R1);
        assert_eq!(edits[0].target_position, 4);
    }

    #[test]
    fn matching_bases_never_change_and_qualities_never_rise() {
        let (mate1, mate2) = pair();
        // R1 confident, R2 doubtful, but the bases agree.
        let mut scores2 = vec![2u8; mate2.len()];
        scores2[0] = 2;
        let (edits, summary) = plan_pair(
            stage(),
            &mate1,
            &quality(&vec![40u8; mate1.len()]),
            &mate2,
            &quality(&scores2),
        );
        assert!(edits.is_empty(), "agreeing bases are never touched");
        assert_eq!(summary.mismatches, 0);
    }

    #[test]
    fn correction_applies_only_to_the_targeted_mate() {
        let (mate1, mut mate2) = pair();
        let position2 = mate2.len() - 1 - 6;
        mate2[position2] = if mate2[position2] == b'A' { b'C' } else { b'A' };
        let mut scores2 = vec![35u8; mate2.len()];
        scores2[position2] = 6;
        let quality1 = quality(&vec![35u8; mate1.len()]);
        let quality2 = quality(&scores2);
        let (edits, _) = plan_pair(stage(), &mate1, &quality1, &mate2, &quality2);

        let mut scratch = CorrectionScratch {
            edits,
            ..CorrectionScratch::default()
        };
        let original1 = ReadView::unchecked(&mate1, Some(&quality1));
        let original2 = ReadView::unchecked(&mate2, Some(&quality2));

        let (borrowed1, corrected2) = scratch.apply_pair(original1, original2);
        assert!(
            std::ptr::eq(borrowed1.sequence, mate1.as_slice()),
            "an uncorrected mate stays borrowed"
        );
        assert_eq!(corrected2.sequence[position2], complement(mate1[6]));
        assert_eq!(corrected2.quality.unwrap()[position2], quality1[6]);
        assert_eq!(corrected2.sequence.len(), mate2.len(), "length unchanged");
        assert_eq!(
            corrected2.quality.unwrap().len(),
            quality2.len(),
            "quality length unchanged"
        );
        // Every other position is untouched.
        for index in 0..mate2.len() {
            if index != position2 {
                assert_eq!(corrected2.sequence[index], mate2[index]);
                assert_eq!(corrected2.quality.unwrap()[index], quality2[index]);
            }
        }
    }

    #[test]
    fn negative_offsets_are_corrected_with_the_right_coordinates() {
        // R2 reads through past R1's 5' end, so the overlap offset is negative.
        let insert = b"ACGTTGCAAGGCCTTAACGTACGTTGCAAGG".to_vec();
        let mate1 = insert[4..].to_vec(); // R1 starts inside the insert
        let mut mate2 = revcomp(&insert);
        let params = OverlapParams {
            min_overlap: 20,
            max_error_rate: 0.2,
        };
        let overlap = find_overlap(&mate1, &mate2, params).expect("overlap");
        assert!(overlap.offset < 0, "offset should be negative: {overlap:?}");

        // Break R2 at the position aligned to R1[2].
        let (_, position2) = overlap
            .aligned_positions(mate2.len())
            .find(|&(position1, _)| position1 == 2)
            .expect("R1[2] is aligned");
        let truth = mate2[position2];
        mate2[position2] = if truth == b'A' { b'C' } else { b'A' };
        let mut scores2 = vec![35u8; mate2.len()];
        scores2[position2] = 8;
        let quality1 = quality(&vec![35u8; mate1.len()]);
        let quality2 = quality(&scores2);

        let overlap = find_overlap(&mate1, &mate2, params).expect("overlap");
        let mut edits = Vec::new();
        let summary = plan(
            stage(),
            ReadView::unchecked(&mate1, Some(&quality1)),
            ReadView::unchecked(&mate2, Some(&quality2)),
            overlap,
            &mut edits,
        );
        assert_eq!(summary.corrected_r2, 1, "{summary:?}");
        assert_eq!(edits[0].target_position, position2);
        assert_eq!(edits[0].donor_position, 2);
        assert_eq!(edits[0].corrected_base, truth);
    }

    #[test]
    fn first_and_last_overlap_positions_are_included() {
        let (mate1, mut mate2) = pair();
        let last_r1 = mate1.len() - 1;
        for &position1 in &[0usize, last_r1] {
            let position2 = mate2.len() - 1 - position1;
            mate2[position2] = if mate2[position2] == b'A' { b'C' } else { b'A' };
        }
        let mut scores2 = vec![35u8; mate2.len()];
        scores2[mate2.len() - 1] = 4;
        scores2[0] = 4;
        let (edits, summary) = plan_pair(
            stage(),
            &mate1,
            &quality(&vec![35u8; mate1.len()]),
            &mate2,
            &quality(&scores2),
        );
        assert_eq!(summary.corrected_r2, 2, "both boundary positions corrected");
        assert_eq!(edits.len(), 2);
    }
}
