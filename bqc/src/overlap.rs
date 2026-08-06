// Copyright (c) 2026 The Hiller Lab at the Senckenberg Gesellschaft für Naturforschung
// Distributed under the terms of the GNU General Public License, Version 3.0.

//! Paired-end overlap analysis.
//!
//! R1 and R2 are read from opposite ends of the same insert, so the reverse
//! complement of R2 must align to R1 wherever the two mates genuinely overlap.
//! A successful alignment reveals the insert length `L`: any bases past `L` in
//! either mate are adapter read-through and can be trimmed without knowing the
//! adapter sequence.
//!
//! # Contract
//!
//! Every placement of reverse-complemented R2 relative to R1 is scored by its
//! Hamming distance over the aligned region (`N` never matches). A placement is
//! accepted when
//!
//! ```text
//! overlap    >= min_overlap            (default 30)
//! mismatches <= floor(overlap * max_error_rate)
//! ```
//!
//! Among all accepted placements the winner is chosen by, in order: longest
//! overlap, fewest mismatches, smallest offset. Length wins first because a
//! short, perfect alignment is cheap to fake — an adapter's homopolymer tail
//! aligning against itself — while a long, nearly-perfect alignment is almost
//! always the true shared insert. The inferred insert length is
//! `offset + len(R2)`: the reverse complement of R2 ends at the last insert
//! base, so its 3' coordinate is exactly one past the insert.

use serde::Serialize;

use crate::error::{Error, Result};

/// Overlap alignment parameters.
#[derive(Debug, Clone, Copy, Serialize)]
pub struct OverlapParams {
    /// Minimum aligned region accepted as evidence.
    pub min_overlap: usize,
    /// Maximum mismatch fraction over the aligned region.
    pub max_error_rate: f64,
}

impl Default for OverlapParams {
    fn default() -> Self {
        Self {
            // Thirty aligned bases make a chance alignment with at most 10%
            // mismatches unlikely; this matches fastp's overlap requirement.
            min_overlap: 30,
            max_error_rate: 0.10,
        }
    }
}

impl OverlapParams {
    /// Validates overlap parameters.
    pub fn validate(self) -> Result<Self> {
        if self.min_overlap == 0 {
            return Err(Error::config(
                "--paired-overlap-min-overlap must be at least 1",
            ));
        }
        if !(0.0..=1.0).contains(&self.max_error_rate) {
            return Err(Error::config(format!(
                "--max-error-rate must be within 0.0..=1.0 (got {})",
                self.max_error_rate
            )));
        }
        Ok(self)
    }
}

/// A successful overlap alignment between R1 and reverse-complemented R2.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Overlap {
    /// Position of reverse-complemented R2 relative to the start of R1.
    /// Negative when R2 reads through the insert and past R1's 5' end.
    pub offset: isize,
    /// Number of aligned bases.
    pub overlap: usize,
    /// Mismatches inside the aligned region.
    pub mismatches: usize,
    /// Inferred insert length: `offset + len(R2)`.
    pub insert_length: usize,
}

/// The Watson–Crick complement of a base; ambiguous bases complement to `N`.
#[inline]
#[must_use]
pub fn complement(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        _ => b'N',
    }
}

/// Finds the winning overlap alignment between R1 and R2, if any.
#[must_use]
pub fn find_overlap(r1: &[u8], r2: &[u8], params: OverlapParams) -> Option<Overlap> {
    let (len1, len2) = (r1.len(), r2.len());
    if len1 < params.min_overlap || len2 < params.min_overlap {
        return None;
    }
    // Offsets leaving at least `min_overlap` aligned bases: reverse-complemented
    // R2 sits at `offset`, so the aligned region is [max(0, o), min(len1, o+len2)).
    // Read lengths are bounded by physical memory, so they always fit an
    // isize; `cast_signed` documents the intent.
    let (len1s, len2s) = (len1.cast_signed(), len2.cast_signed());
    let min_overlap = params.min_overlap.cast_signed();
    let first = min_overlap - len2s;
    let last = len1s - min_overlap;
    let mut best: Option<Overlap> = None;
    for offset in first..=last {
        let lo = offset.max(0) as usize;
        let hi = (offset + len2s).min(len1s) as usize;
        let overlap = hi - lo;
        if overlap < params.min_overlap {
            continue;
        }
        // The cast truncates towards zero, which for a non-negative product is
        // `floor`; see `AdapterParams::error_limit`.
        let limit = (overlap as f64 * params.max_error_rate) as usize;
        let (first_r1, first_r2, length) = aligned_region(offset, overlap, len2);
        let mut mismatches = 0usize;
        for step in 0..length {
            if bases_disagree(r1[first_r1 + step], r2[first_r2 - step]) {
                mismatches += 1;
                if mismatches > limit {
                    break;
                }
            }
        }
        if mismatches > limit {
            continue;
        }
        // Longest overlap wins, then fewest mismatches; scan order keeps the
        // smallest offset on a full tie.
        let improves = best
            .is_none_or(|current| (overlap, current.mismatches) > (current.overlap, mismatches));
        if improves {
            best = Some(Overlap {
                offset,
                overlap,
                mismatches,
                insert_length: (offset + len2s) as usize,
            });
        }
    }
    best
}

/// The aligned region of an overlap, as
/// `(first R1 position, first R2 position, length)`.
///
/// R1 is walked forwards from the first coordinate and R2 *backwards* from the
/// second, because the mates are sequenced in opposite orientations. This is the
/// only place the offset arithmetic lives: both the overlap search and base
/// correction derive their coordinates from it.
#[inline]
#[must_use]
fn aligned_region(offset: isize, overlap: usize, r2_len: usize) -> (usize, usize, usize) {
    // R1 coordinate `c` pairs with `r2[r2_len - 1 - (c - offset)]`. The aligned
    // region starts at `max(offset, 0)` in R1, so `max(-offset, 0)` bases of
    // reverse-complemented R2 lie before R1's start and are skipped.
    let first_r1 = offset.max(0) as usize;
    let skipped = (-offset).max(0) as usize;
    debug_assert!(skipped + overlap <= r2_len);
    (first_r1, r2_len - 1 - skipped, overlap)
}

impl Overlap {
    /// Iterates the aligned `(r1_position, r2_position)` coordinate pairs.
    ///
    /// Allocation free, and the single source of the offset arithmetic.
    pub fn aligned_positions(&self, r2_len: usize) -> impl Iterator<Item = (usize, usize)> + '_ {
        let (first_r1, first_r2, length) = aligned_region(self.offset, self.overlap, r2_len);
        (0..length).map(move |step| (first_r1 + step, first_r2 - step))
    }
}

/// Whether two aligned bases disagree.
///
/// `N` never agrees with anything, including another `N`, which matches the rule
/// the adapter matcher uses.
#[inline]
#[must_use]
pub fn bases_disagree(r1_base: u8, r2_base: u8) -> bool {
    let expected = complement(r2_base);
    r1_base != expected || r1_base == b'N'
}

/// The 3' trim coordinate implied by an overlap for a read of `length` bases.
///
/// Bases beyond the insert are adapter read-through; a read no longer than the
/// insert is kept whole.
#[must_use]
pub fn insert_boundary(overlap: Overlap, length: usize) -> usize {
    overlap.insert_length.min(length)
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;

    /// Reverse-complement helper for building fixtures.
    pub(crate) fn revcomp(sequence: &[u8]) -> Vec<u8> {
        sequence.iter().rev().map(|&b| complement(b)).collect()
    }

    /// A deterministic non-repetitive sequence (periodic inserts would align
    /// at many phase-shifted offsets and break the fixtures' geometry).
    fn sequence(seed: u64, length: usize) -> Vec<u8> {
        let mut state = seed | 1;
        (0..length)
            .map(|_| {
                state ^= state >> 12;
                state ^= state << 25;
                state ^= state >> 27;
                state = state.wrapping_mul(0x2545_F491_4F6C_DD1D);
                b"ACGT"[(state >> 62) as usize]
            })
            .collect()
    }

    /// R1 and R2 reads of length `insert + tail` for an insert of `insert`
    /// bases, with `tail` bases of adapter read-through on each mate.
    fn pair(insert: usize, tail: usize) -> (Vec<u8>, Vec<u8>) {
        let insert_seq = sequence(0xFACE, insert);
        let r1 = [insert_seq.as_slice(), &sequence(0xBEEF, tail)].concat();
        let r2 = [revcomp(&insert_seq).as_slice(), &sequence(0xDEAD, tail)].concat();
        (r1, r2)
    }

    #[test]
    fn read_through_pairs_yield_the_insert_length() {
        // 100-base insert, both reads 130 bases: 30 bases of read-through.
        let (r1, r2) = pair(100, 30);
        let overlap = find_overlap(&r1, &r2, OverlapParams::default()).unwrap();
        assert_eq!(overlap.insert_length, 100);
        assert_eq!(overlap.mismatches, 0);
        assert_eq!(overlap.overlap, 100);
        assert_eq!(insert_boundary(overlap, r1.len()), 100);
        assert_eq!(insert_boundary(overlap, r2.len()), 100);
    }

    #[test]
    fn inserts_longer_than_the_reads_trim_nothing() {
        // 200-base insert, 130-base reads: the overlap is 60 bases, and both
        // reads are entirely insert.
        let (r1, r2) = pair(200, 0);
        let r1 = &r1[..130];
        let r2 = &r2[..130];
        let overlap = find_overlap(r1, r2, OverlapParams::default()).unwrap();
        assert_eq!(overlap.insert_length, 200);
        assert_eq!(overlap.overlap, 60);
        assert_eq!(insert_boundary(overlap, 130), 130);
    }

    #[test]
    fn a_few_mismatches_are_tolerated() {
        let (mut r1, r2) = pair(100, 30);
        r1[10] = b'N'; // one ambiguous base inside the overlap
        let params = OverlapParams {
            max_error_rate: 0.05,
            ..OverlapParams::default()
        };
        let overlap = find_overlap(&r1, &r2, params).unwrap();
        assert_eq!(overlap.insert_length, 100);
        assert_eq!(overlap.mismatches, 1);

        let strict = OverlapParams {
            max_error_rate: 0.0,
            ..params
        };
        assert!(find_overlap(&r1, &r2, strict).is_none());
    }

    #[test]
    fn unrelated_reads_do_not_overlap() {
        let r1 = sequence(1, 150);
        let r2 = sequence(99, 150);
        assert!(find_overlap(&r1, &r2, OverlapParams::default()).is_none());
    }

    #[test]
    fn reads_shorter_than_the_minimum_overlap_are_skipped() {
        let (r1, r2) = pair(20, 0);
        assert!(find_overlap(&r1, &r2, OverlapParams::default()).is_none());
    }

    #[test]
    fn the_exact_minimum_overlap_is_accepted() {
        // A 100-base insert read by two 65-base mates overlaps by exactly 30.
        let (r1_full, r2_full) = pair(100, 0);
        let r1 = &r1_full[..65];
        let r2 = &r2_full[..65];
        let params = OverlapParams {
            min_overlap: 30,
            ..OverlapParams::default()
        };
        let overlap = find_overlap(r1, r2, params).unwrap();
        assert_eq!(overlap.overlap, 30);
        assert_eq!(overlap.insert_length, 100);
    }

    #[test]
    fn a_short_perfect_alignment_does_not_beat_the_true_longer_one() {
        // The true alignment has one mismatch over 100 bases; a homopolymer
        // tail offers a spurious perfect 30-base alignment. Length wins.
        let (mut r1, r2) = pair(100, 30);
        // G/C homopolymer read-through tails create the spurious candidate.
        for base in r1.iter_mut().skip(100) {
            *base = b'G';
        }
        let mut r2 = r2;
        for base in r2.iter_mut().skip(100) {
            *base = b'C';
        }
        r1[10] = if r1[10] == b'A' { b'C' } else { b'A' };
        let overlap = find_overlap(&r1, &r2, OverlapParams::default()).unwrap();
        assert_eq!(overlap.insert_length, 100);
        assert_eq!(overlap.mismatches, 1);
    }

    #[test]
    fn asymmetric_read_through_trims_each_mate_independently() {
        // 100-base insert; R1 120 bases (20 of read-through), R2 110 (10).
        let (r1_full, r2_full) = pair(100, 0);
        let r1 = [r1_full.as_slice(), &[b'G'; 20]].concat();
        let r2 = [r2_full.as_slice(), &[b'C'; 10]].concat();
        let overlap = find_overlap(&r1, &r2, OverlapParams::default()).unwrap();
        assert_eq!(overlap.insert_length, 100);
        assert_eq!(insert_boundary(overlap, r1.len()), 100);
        assert_eq!(insert_boundary(overlap, r2.len()), 100);
    }

    #[test]
    fn validation_rejects_degenerate_parameters() {
        assert!(
            OverlapParams {
                min_overlap: 0,
                ..OverlapParams::default()
            }
            .validate()
            .is_err()
        );
        assert!(
            OverlapParams {
                min_overlap: 30,
                max_error_rate: 1.5,
            }
            .validate()
            .is_err()
        );
        assert!(
            OverlapParams {
                min_overlap: 30,
                max_error_rate: f64::NAN,
            }
            .validate()
            .is_err()
        );
    }
}
