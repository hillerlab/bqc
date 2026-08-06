// Copyright (c) 2026 The Hiller Lab at the Senckenberg Gesellschaft für Naturforschung
// Distributed under the terms of the GNU General Public License, Version 3.0.

//! Three-prime adapter matching.
//!
//! # Contract
//!
//! An adapter is trimmed from the 3' end of a read. For every candidate start
//! coordinate `p` in the read, the adapter prefix is compared against
//! `read[p..p + overlap]` where `overlap = min(adapter_len, read_len - p)`. A
//! candidate matches when
//!
//! ```text
//! overlap >= min_overlap
//! errors  <= floor(overlap * max_error_rate)
//! errors  <= max_errors                      (when configured)
//! ```
//!
//! Among all matches the winner is chosen by, in order: earliest coordinate,
//! fewest errors, longest overlap, adapter declaration order. The read is then
//! shortened to `p`.
//!
//! `N` never matches — in the read or in the adapter — and counts as one error.
//! This is deliberately conservative: a poly-N tail would otherwise match any
//! adapter.
//!
//! # Indel-aware matching
//!
//! By default only substitutions are counted. With `--allow-indels`, each
//! candidate start is evaluated by a banded edit-distance alignment between an
//! adapter prefix `adapter[..i]` and a read window `read[p..p + j]`: insertions
//! and deletions each cost one error, the overlap reported is the number of
//! aligned adapter bases `i`, and the same error budget applies. Every
//! substitution-only alignment is also an edit-distance alignment with the same
//! cost, so the banded scan subsumes the cheaper one and replaces it.
//!
//! Selection keeps the earliest coordinate, with one necessary refinement.
//! Edit distance casts *shadows*: an occurrence at `p` is usually also
//! reachable at `p - d` by consuming the `d` extra read bases as insertions, at
//! a cost of at least `d` errors. Blindly taking the earliest match would
//! therefore shave a base or two off otherwise clean matches. Since a shadow at
//! distance `d` costs at least `d` errors, the true occurrence lies within
//! `errors` bases of the earliest accepted match; that window is scanned and
//! the candidate with the most matched adapter bases wins, ties going to the
//! earliest coordinate. Matches beyond the window are genuinely separate
//! occurrences, so the earliest-coordinate rule still removes the longest
//! contaminated suffix.

use std::fs;
use std::path::Path;

use serde::Serialize;

use crate::error::{Error, Result};

/// One named adapter sequence, normalized to uppercase `ACGTN`.
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct Adapter {
    pub name: String,
    #[serde(serialize_with = "crate::report::serialize_bytes")]
    pub sequence: Vec<u8>,
}

impl Adapter {
    /// Validates and normalizes an adapter sequence.
    pub fn new(name: impl Into<String>, sequence: &[u8]) -> Result<Self> {
        let name = name.into();
        if sequence.is_empty() {
            return Err(Error::InvalidAdapter(format!("adapter '{name}' is empty")));
        }
        let mut normalized = Vec::with_capacity(sequence.len());
        for (offset, &base) in sequence.iter().enumerate() {
            let base = base.to_ascii_uppercase();
            if !matches!(base, b'A' | b'C' | b'G' | b'T' | b'N') {
                return Err(Error::InvalidAdapter(format!(
                    "adapter '{name}' contains unsupported symbol {:?} at position {}; \
                     only A, C, G, T and N are supported",
                    char::from(base),
                    offset + 1
                )));
            }
            normalized.push(base);
        }
        Ok(Self {
            name,
            sequence: normalized,
        })
    }

    #[must_use]
    pub fn len(&self) -> usize {
        self.sequence.len()
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }
}

/// Matching thresholds shared by every adapter.
#[derive(Debug, Clone, Copy, Serialize)]
pub struct AdapterParams {
    pub min_overlap: usize,
    pub max_error_rate: f64,
    pub max_errors: Option<usize>,
    /// Count insertions and deletions in addition to substitutions.
    pub allow_indels: bool,
}

impl Default for AdapterParams {
    fn default() -> Self {
        Self {
            min_overlap: 8,
            max_error_rate: 0.10,
            max_errors: None,
            allow_indels: false,
        }
    }
}

impl AdapterParams {
    /// Maximum number of substitutions tolerated over `overlap` bases.
    ///
    /// Casting `f64` to an integer in Rust already truncates towards zero, which
    /// for a non-negative product is exactly `floor`. Writing the cast alone
    /// rather than `.floor() as usize` gives bit-identical results and keeps this
    /// out of libm: `floor` does not inline, and this is evaluated once per
    /// candidate coordinate per adapter, which made it 6–11% of a run.
    #[inline]
    #[must_use]
    pub fn error_limit(&self, overlap: usize) -> usize {
        let rate_limit = (overlap as f64 * self.max_error_rate) as usize;
        match self.max_errors {
            Some(explicit) => explicit.min(rate_limit),
            None => rate_limit,
        }
    }
}

/// A single accepted adapter match.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AdapterHit {
    /// Index of the matching adapter within its mate's declaration list.
    pub adapter: usize,
    /// Coordinate at which the adapter starts, relative to the retained span.
    pub start: usize,
    pub errors: usize,
    /// Number of aligned adapter bases. With `--allow-indels` this can differ
    /// from the number of read bases the alignment consumed.
    pub overlap: usize,
    /// Number of *read* bases the alignment covered. Equal to `overlap` unless
    /// indels shifted the two apart. A five-prime caller needs this one: the
    /// insert begins at `start + consumed`.
    pub consumed: usize,
}

/// Compiled adapter configuration for both mates.
#[derive(Debug, Clone, Serialize)]
pub struct AdapterStage {
    pub r1: Vec<Adapter>,
    pub r2: Vec<Adapter>,
    pub params: AdapterParams,
    /// Paired-overlap inference, enabled independently of explicit adapters.
    pub paired_overlap: Option<crate::overlap::OverlapParams>,
}

impl AdapterStage {
    /// Builds the stage, rejecting adapters that can never match.
    pub fn new(
        r1: Vec<Adapter>,
        r2: Vec<Adapter>,
        params: AdapterParams,
        paired_overlap: Option<crate::overlap::OverlapParams>,
    ) -> Result<Self> {
        if params.min_overlap == 0 {
            return Err(Error::config("--min-overlap must be at least 1"));
        }
        if !(0.0..=1.0).contains(&params.max_error_rate) {
            return Err(Error::config(format!(
                "--max-error-rate must be within 0.0..=1.0 (got {})",
                params.max_error_rate
            )));
        }
        if let Some(overlap) = paired_overlap {
            overlap.validate()?;
        }
        if r1.is_empty() && r2.is_empty() && paired_overlap.is_none() {
            return Err(Error::config(
                "adapter removal requires --adapter-r1, --adapter-r2, --adapter-fasta, \
                 --paired-overlap or --auto-detect",
            ));
        }
        for adapter in r1.iter().chain(&r2) {
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
            r1,
            r2,
            params,
            paired_overlap,
        })
    }

    /// Adapters declared for `mate`.
    #[must_use]
    pub fn adapters(&self, mate: crate::process::Mate) -> &[Adapter] {
        match mate {
            crate::process::Mate::R1 => &self.r1,
            crate::process::Mate::R2 => &self.r2,
        }
    }

    /// Finds the winning adapter match in `sequence`, if any.
    #[must_use]
    pub fn find(&self, mate: crate::process::Mate, sequence: &[u8]) -> Option<AdapterHit> {
        find_three_prime(self.adapters(mate), self.params, sequence)
    }
}

/// Counts substitutions between an adapter prefix and a read window.
///
/// Returns `None` as soon as the count exceeds `limit`.
///
/// Specialising this on whether the adapter contains `N` — hoisting the `N` test
/// out of the loop — was measured and made no difference on any fixture, so the
/// single readable loop stays.
#[inline]
fn count_errors(adapter: &[u8], read: &[u8], limit: usize) -> Option<usize> {
    debug_assert_eq!(adapter.len(), read.len());
    let mut errors = 0;
    for (&a, &r) in adapter.iter().zip(read) {
        // `N` is never a match, even against another `N`.
        if a != r || a == b'N' {
            errors += 1;
            if errors > limit {
                return None;
            }
        }
    }
    Some(errors)
}

/// Verifies one adapter at an exact coordinate under the three-prime rules.
///
/// This is the per-candidate step of [`find_three_prime`], exposed for callers
/// that already know where an adapter would start — adapter detection looks the
/// coordinate up from a seed index rather than scanning for it. The acceptance
/// rules are identical, so a hit here is a hit there.
#[must_use]
pub fn verify_at(
    adapter: &Adapter,
    params: AdapterParams,
    sequence: &[u8],
    start: usize,
) -> Option<AdapterHit> {
    let overlap = adapter.len().min(sequence.len().checked_sub(start)?);
    if overlap < params.min_overlap {
        return None;
    }
    let errors = count_errors(
        &adapter.sequence[..overlap],
        &sequence[start..start + overlap],
        params.error_limit(overlap),
    )?;
    Some(AdapterHit {
        adapter: 0,
        start,
        errors,
        overlap,
        consumed: overlap,
    })
}

/// Finds a five-prime adapter near the start of `sequence`.
///
/// The mirror image of [`find_three_prime`]'s acceptance rule: the adapter must
/// begin within `max_offset` bases of the read's 5' end and be consumed whole,
/// unless the read's own 5' end cuts it short. Everything else — the error
/// budget, the `N` rule, indel awareness, tie-breaking — is the three-prime
/// machinery reused unchanged.
///
/// The returned hit's `start` is where the adapter begins; the insert therefore
/// begins at `start + overlap` for a fully consumed adapter.
#[must_use]
pub fn find_five_prime(
    adapters: &[Adapter],
    params: AdapterParams,
    sequence: &[u8],
    max_offset: usize,
) -> Option<AdapterHit> {
    let len = sequence.len();
    let last_start = max_offset.min(len);
    let mut best: Option<AdapterHit> = None;
    // The banded alignment needs scratch space; allocated once, and only when
    // indel-aware matching is enabled.
    let mut scratch = if params.allow_indels {
        vec![u32::MAX; 2 * (len + 1)]
    } else {
        Vec::new()
    };
    for start in 0..=last_start {
        for (index, adapter) in adapters.iter().enumerate() {
            // The adapter must fit in what remains, or be truncated by the 5' end
            // of the read itself — which only happens at `start == 0`.
            let available = len - start;
            let candidate = if params.allow_indels {
                // The banded alignment's acceptance rule already covers a fully
                // consumed adapter, which is what a 5' flank is.
                indel_best_at(
                    std::slice::from_ref(adapter),
                    params,
                    sequence,
                    start,
                    &mut scratch,
                )
            } else if available >= adapter.len() || start == 0 {
                verify_at(adapter, params, sequence, start)
            } else {
                None
            };
            let Some(mut candidate) = candidate else {
                continue;
            };
            candidate.adapter = index;
            // Earliest offset wins: the loop returns as soon as a start
            // produced any hit. Within a start, `supersedes` decides.
            let improves = best.is_none_or(|current| supersedes(&candidate, &current));
            if improves {
                best = Some(candidate);
            }
        }
        if best.is_some() {
            return best;
        }
    }
    None
}

/// Whether `candidate` beats `current` at the same coordinate: fewer errors
/// wins, then longer overlap, then declaration order (preserved by replacing
/// only on improvement).
fn supersedes(candidate: &AdapterHit, current: &AdapterHit) -> bool {
    (candidate.errors, current.overlap) < (current.errors, candidate.overlap)
}

/// Scans `sequence` for the winning three-prime adapter match.
#[must_use]
pub fn find_three_prime(
    adapters: &[Adapter],
    params: AdapterParams,
    sequence: &[u8],
) -> Option<AdapterHit> {
    if params.allow_indels {
        scan_with_indels(adapters, params, sequence)
    } else {
        scan_substitutions(adapters, params, sequence)
    }
}

/// Aligned adapter bases that matched: the evidence a candidate carries.
#[inline]
fn matched_bases(hit: &AdapterHit) -> usize {
    hit.overlap - hit.errors
}

/// Whether `candidate` beats `current` under the indel-mode ordering:
/// most matched adapter bases, then earliest coordinate, then fewest errors,
/// then longest overlap, then declaration order.
pub(crate) fn indel_improves(candidate: &AdapterHit, current: &AdapterHit) -> bool {
    (
        std::cmp::Reverse(matched_bases(candidate)),
        candidate.start,
        candidate.errors,
        std::cmp::Reverse(candidate.overlap),
        candidate.adapter,
    ) < (
        std::cmp::Reverse(matched_bases(current)),
        current.start,
        current.errors,
        std::cmp::Reverse(current.overlap),
        current.adapter,
    )
}

/// Scans `sequence` allowing insertions and deletions.
///
/// The earliest coordinate still wins, exactly as for substitutions — a read
/// carrying two adapter copies must be trimmed at the first one. But edit
/// distance casts *shadows*: a true occurrence at `p` is usually also reachable
/// at `p - d` by consuming the `d` extra read bases as insertions, which costs
/// at least `d` errors. Taking the earliest match blindly would therefore
/// over-trim by a base or two — worst of all on perfectly clean adapter
/// matches, where a shadow is the only reason the read would lose extra bases.
///
/// Because a shadow at distance `d` pays at least `d` errors, the true
/// occurrence lies within `errors` bases of the earliest accepted match. That
/// window is scanned and the candidate carrying the most matched adapter bases
/// wins; ties fall back to the earliest coordinate. Everything past the window
/// is a genuinely different occurrence and is left alone.
fn scan_with_indels(
    adapters: &[Adapter],
    params: AdapterParams,
    sequence: &[u8],
) -> Option<AdapterHit> {
    let len = sequence.len();
    if adapters.is_empty() || len < params.min_overlap {
        return None;
    }
    let last_start = len - params.min_overlap;
    let mut scratch = vec![u32::MAX; 2 * (len + 1)];
    let mut best: Option<AdapterHit> = None;
    let mut window_end = last_start;
    let mut start = 0;
    while start <= window_end {
        if let Some(candidate) = indel_best_at(adapters, params, sequence, start, &mut scratch) {
            if best.is_none() {
                // First accepted match: it bounds how far the true occurrence
                // can be from here, which ends the scan early.
                window_end = (start + candidate.errors).min(last_start);
            }
            if best
                .as_ref()
                .is_none_or(|current| indel_improves(&candidate, current))
            {
                best = Some(candidate);
            }
        }
        start += 1;
    }
    best
}

/// Scans `sequence` counting only substitutions.
fn scan_substitutions(
    adapters: &[Adapter],
    params: AdapterParams,
    sequence: &[u8],
) -> Option<AdapterHit> {
    let len = sequence.len();
    if adapters.is_empty() || len < params.min_overlap {
        return None;
    }
    // Positions beyond this point leave fewer than `min_overlap` bases.
    let last_start = len - params.min_overlap;
    for start in 0..=last_start {
        let mut best: Option<AdapterHit> = None;
        for (index, adapter) in adapters.iter().enumerate() {
            let overlap = adapter.len().min(len - start);
            if overlap < params.min_overlap {
                continue; // adapter shorter than the required overlap
            }
            let limit = params.error_limit(overlap);
            let Some(errors) = count_errors(
                &adapter.sequence[..overlap],
                &sequence[start..start + overlap],
                limit,
            ) else {
                continue;
            };
            let candidate = AdapterHit {
                adapter: index,
                start,
                errors,
                overlap,
                consumed: overlap,
            };
            let improves = best.is_none_or(|current| supersedes(&candidate, &current));
            if improves {
                best = Some(candidate);
            }
        }
        if best.is_some() {
            return best;
        }
    }
    None
}

/// Aligns every adapter against the read suffix at `start`, returning the best
/// edit-distance hit at that coordinate.
///
/// The DP aligns an adapter prefix `adapter[..i]` with a read window
/// `read[start..start + j]`. Exactly as in the substitution matcher, the read
/// must end in the adapter: a cell is a candidate only when it consumes the
/// whole read suffix (`j == m`, the adapter overhangs the read's 3' end) or
/// the whole adapter (`i == adapter_len`, the read overhangs the adapter's 3'
/// end). A candidate matches when `i >= min_overlap` and
/// `distance <= error_limit(i)`, where `i` — the number of aligned adapter
/// bases — is the reported overlap. Candidates are compared by fewest errors,
/// then longest adapter prefix, then least read consumed.
pub(crate) fn indel_best_at(
    adapters: &[Adapter],
    params: AdapterParams,
    sequence: &[u8],
    start: usize,
    scratch: &mut [u32],
) -> Option<AdapterHit> {
    let read = &sequence[start..];
    let m = read.len();
    let mut best: Option<AdapterHit> = None;
    for (index, adapter) in adapters.iter().enumerate() {
        // No acceptable alignment can have more errors than the budget of the
        // longest candidate, and edit distance never falls below |i - j|, so
        // the band |i - j| <= k loses nothing.
        let k = params.error_limit(adapter.len());
        let rows = adapter.len().min(m + k);
        if rows < params.min_overlap {
            continue;
        }
        let (mut prev_row, mut curr_row) = scratch.split_at_mut(m + 1);
        // Row 0: only the origin is reachable. Leading insertions are
        // forbidden — the adapter's 5' end is anchored at `start`, so the
        // first aligned pair must involve adapter[0]. (Otherwise a hit at
        // `start + d` would masquerade as a hit at `start` with d errors.)
        for (j, cell) in prev_row.iter_mut().enumerate() {
            *cell = if j == 0 { 0 } else { u32::MAX };
        }
        // Best accepted cell so far: (errors, adapter bases, read bases).
        let mut cell_best: Option<(usize, usize, usize)> = None;
        let consider =
            |errors: u32, i: usize, j: usize, cell_best: &mut Option<(usize, usize, usize)>| {
                let errors = errors as usize;
                if errors > params.error_limit(i) {
                    return;
                }
                // Fewest errors, then the longest adapter prefix, then the *most*
                // read bases consumed. The last term is invisible to three-prime
                // trimming, which only needs the coordinate, but it decides where
                // a five-prime flank ends: preferring a shorter read span would
                // leave adapter bases at the start of the retained insert.
                let improves = cell_best.is_none_or(|(e, ai, rj)| {
                    errors < e || (errors == e && (i > ai || (i == ai && j > rj)))
                });
                if improves {
                    *cell_best = Some((errors, i, j));
                }
            };
        for i in 1..=rows {
            let a = adapter.sequence[i - 1];
            let j_lo = i.saturating_sub(k).max(1);
            let j_hi = (i + k).min(m);
            curr_row[0] = if i <= k { i as u32 } else { u32::MAX };
            // Only the band is computed. The cells immediately outside it are
            // set unreachable because the next row reads them; everything
            // further out is never read again.
            if j_lo > 1 {
                curr_row[j_lo - 1] = u32::MAX;
            }
            for j in j_lo..=j_hi {
                let substitution = u32::from(a != read[j - 1] || a == b'N');
                curr_row[j] = (prev_row[j].saturating_add(1))
                    .min(curr_row[j - 1].saturating_add(1))
                    .min(prev_row[j - 1].saturating_add(substitution));
            }
            if j_hi < m {
                curr_row[j_hi + 1] = u32::MAX;
            }
            // Cells on the read's 3' edge: the adapter may hang over.
            if i >= params.min_overlap && i.abs_diff(m) <= k {
                consider(curr_row[m], i, m, &mut cell_best);
            }
            // Cells on the adapter's 3' edge: the read may hang over. Only
            // in-band cells hold a value from this row.
            if i == adapter.len() {
                for (offset, &distance) in curr_row[j_lo..=j_hi].iter().enumerate() {
                    consider(distance, i, j_lo + offset, &mut cell_best);
                }
            }
            std::mem::swap(&mut prev_row, &mut curr_row);
        }
        if let Some((errors, overlap, consumed)) = cell_best {
            let candidate = AdapterHit {
                adapter: index,
                start,
                errors,
                overlap,
                consumed,
            };
            let improves = best.is_none_or(|current| {
                (candidate.errors, current.overlap) < (current.errors, candidate.overlap)
            });
            if improves {
                best = Some(candidate);
            }
        }
    }
    best
}

/// Reads adapters from a FASTA file.
///
/// Record names are the identifier up to the first whitespace. Sequences may
/// span multiple lines and are normalized to uppercase.
pub fn read_adapter_fasta(path: &Path) -> Result<Vec<Adapter>> {
    let text = fs::read(path).map_err(|e| Error::read(path, e))?;
    let mut adapters = Vec::new();
    let mut name: Option<String> = None;
    let mut sequence: Vec<u8> = Vec::new();

    for line in text.split(|&b| b == b'\n') {
        let line = line.strip_suffix(b"\r").unwrap_or(line);
        if line.is_empty() {
            continue;
        }
        if line[0] == b'>' {
            if let Some(name) = name.take() {
                adapters.push(Adapter::new(name, &sequence)?);
                sequence.clear();
            }
            let header = String::from_utf8_lossy(&line[1..]);
            let id = header.split_whitespace().next().unwrap_or_default();
            if id.is_empty() {
                return Err(Error::InvalidAdapter(format!(
                    "{}: FASTA record {} has an empty identifier",
                    path.display(),
                    adapters.len() + 1
                )));
            }
            name = Some(id.to_string());
        } else if name.is_some() {
            sequence.extend_from_slice(line);
        } else {
            return Err(Error::InvalidAdapter(format!(
                "{}: sequence data before the first '>' header",
                path.display()
            )));
        }
    }
    if let Some(name) = name {
        adapters.push(Adapter::new(name, &sequence)?);
    }
    if adapters.is_empty() {
        return Err(Error::InvalidAdapter(format!(
            "{} contains no FASTA records",
            path.display()
        )));
    }
    Ok(adapters)
}

#[cfg(test)]
mod tests {
    use super::*;

    const ILLUMINA: &[u8] = b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";

    fn adapters(seqs: &[&[u8]]) -> Vec<Adapter> {
        seqs.iter()
            .enumerate()
            .map(|(i, s)| Adapter::new(format!("a{i}"), s).unwrap())
            .collect()
    }

    fn find(seqs: &[&[u8]], params: AdapterParams, read: &[u8]) -> Option<AdapterHit> {
        find_three_prime(&adapters(seqs), params, read)
    }

    fn default_params() -> AdapterParams {
        AdapterParams::default()
    }

    #[test]
    fn adapter_normalizes_case_and_rejects_symbols() {
        assert_eq!(Adapter::new("a", b"acgtn").unwrap().sequence, b"ACGTN");
        let err = Adapter::new("a", b"ACGU").unwrap_err();
        assert!(format!("{err}").contains("position 4"), "{err}");
        assert!(Adapter::new("a", b"").is_err());
    }

    #[test]
    fn truncating_the_cast_equals_an_explicit_floor() {
        // The hot paths drop `.floor()` and rely on the cast truncating towards
        // zero. Check that over every rate and overlap they can be given,
        // including the values a validator rejects, so the optimisation can
        // never change a decision.
        for rate_step in 0..=100u32 {
            let rate = f64::from(rate_step) / 100.0;
            for overlap in 0..=512usize {
                let with_floor = (overlap as f64 * rate).floor() as usize;
                let truncated = (overlap as f64 * rate) as usize;
                assert_eq!(with_floor, truncated, "rate {rate}, overlap {overlap}");
            }
        }
        for rate in [f64::NAN, f64::INFINITY, -0.5, 1e300] {
            for overlap in [0usize, 1, 33, 512] {
                assert_eq!(
                    (overlap as f64 * rate).floor() as usize,
                    (overlap as f64 * rate) as usize,
                    "rate {rate}, overlap {overlap}"
                );
            }
        }
    }

    #[test]
    fn integer_comparison_matches_a_floored_threshold() {
        // The poly-tail test drops `floor` differently: `m <= floor(x)` is
        // equivalent to `m <= x` for integer `m`.
        for rate_step in 0..=100u32 {
            let rate = f64::from(rate_step) / 100.0;
            for length in 0..=256usize {
                for mismatches in 0..=length.min(24) {
                    let floored = mismatches <= (length as f64 * rate).floor() as usize;
                    let direct = mismatches as f64 <= length as f64 * rate;
                    assert_eq!(floored, direct, "rate {rate}, len {length}, m {mismatches}");
                }
            }
        }
    }

    #[test]
    fn error_limit_uses_floor_of_rate_and_explicit_cap() {
        let params = AdapterParams {
            min_overlap: 8,
            max_error_rate: 0.10,
            max_errors: None,
            allow_indels: false,
        };
        assert_eq!(params.error_limit(8), 0, "floor(0.8) == 0");
        assert_eq!(params.error_limit(10), 1);
        assert_eq!(params.error_limit(33), 3);
        let capped = AdapterParams {
            max_errors: Some(1),
            ..params
        };
        assert_eq!(capped.error_limit(33), 1);
    }

    #[test]
    fn exact_complete_adapter_is_trimmed_at_its_start() {
        let read = [b"ACGTACGTACGTACGTACGT".as_slice(), ILLUMINA].concat();
        let hit = find(&[ILLUMINA], default_params(), &read).unwrap();
        assert_eq!(hit.start, 20);
        assert_eq!(hit.errors, 0);
        assert_eq!(hit.overlap, ILLUMINA.len());
    }

    #[test]
    fn partial_adapter_at_read_end_is_trimmed() {
        let read = [b"ACGTACGTACGTACGTACGT".as_slice(), &ILLUMINA[..10]].concat();
        let hit = find(&[ILLUMINA], default_params(), &read).unwrap();
        assert_eq!(hit.start, 20);
        assert_eq!(hit.overlap, 10);
    }

    #[test]
    fn exact_minimum_overlap_matches_but_one_less_does_not() {
        let params = AdapterParams {
            min_overlap: 8,
            ..default_params()
        };
        let read = [b"TTTTTTTTTTTT".as_slice(), &ILLUMINA[..8]].concat();
        assert_eq!(find(&[ILLUMINA], params, &read).unwrap().start, 12);

        let read = [b"TTTTTTTTTTTT".as_slice(), &ILLUMINA[..7]].concat();
        assert!(find(&[ILLUMINA], params, &read).is_none());
    }

    #[test]
    fn mismatch_inside_budget_matches_and_outside_does_not() {
        // 20 bases of adapter -> floor(20 * 0.1) == 2 errors allowed.
        let mut tail = ILLUMINA[..20].to_vec();
        tail[5] = if tail[5] == b'A' { b'C' } else { b'A' };
        let read = [b"TTTTTTTTTTTT".as_slice(), &tail].concat();
        let hit = find(&[ILLUMINA], default_params(), &read).unwrap();
        assert_eq!((hit.start, hit.errors), (12, 1));

        // Zero-tolerance configuration rejects the same read.
        let strict = AdapterParams {
            max_errors: Some(0),
            ..default_params()
        };
        assert!(find(&[ILLUMINA], strict, &read).is_none());
    }

    #[test]
    fn earliest_coordinate_wins_over_fewer_errors() {
        // The adapter occurs with one mismatch at offset 12 and exactly at 53.
        let mut early = ILLUMINA.to_vec();
        early[3] = if early[3] == b'A' { b'C' } else { b'A' };
        let read = [
            b"TTTTTTTTTTTT".as_slice(),
            &early,
            b"AAAAAAAA".as_slice(),
            ILLUMINA,
        ]
        .concat();
        let hit = find(&[ILLUMINA], default_params(), &read).unwrap();
        assert_eq!(hit.start, 12, "the longest contaminated suffix is removed");
        assert_eq!(hit.errors, 1);
    }

    #[test]
    fn partial_adapters_only_match_at_the_three_prime_end() {
        // A truncated adapter copy in the middle of the read must be compared
        // against the whole adapter, so it is not a match. Only a suffix may
        // match a prefix of the adapter.
        let read = [
            b"TTTTTTTTTTTT".as_slice(),
            &ILLUMINA[..20],
            b"ACGTACGTACGTACGTACGT".as_slice(),
        ]
        .concat();
        assert!(find(&[ILLUMINA], default_params(), &read).is_none());
    }

    #[test]
    fn tie_break_prefers_fewer_errors_then_longer_overlap_then_declaration_order() {
        // Two adapters match at the same coordinate; the second has no errors.
        let mut mismatched = ILLUMINA[..12].to_vec();
        mismatched[0] = b'T';
        let read = [b"CCCCCCCCCCCC".as_slice(), &ILLUMINA[..12]].concat();
        let hit = find(&[&mismatched, &ILLUMINA[..12]], default_params(), &read).unwrap();
        assert_eq!((hit.adapter, hit.errors), (1, 0));

        // Equal errors: the longer overlap wins.
        let hit = find(&[&ILLUMINA[..8], &ILLUMINA[..12]], default_params(), &read).unwrap();
        assert_eq!((hit.adapter, hit.overlap), (1, 12));

        // Fully identical candidates: declaration order wins.
        let hit = find(&[&ILLUMINA[..12], &ILLUMINA[..12]], default_params(), &read).unwrap();
        assert_eq!(hit.adapter, 0);
    }

    #[test]
    fn adapter_longer_than_read_still_matches_a_prefix() {
        let read = ILLUMINA[..12].to_vec();
        let hit = find(&[ILLUMINA], default_params(), &read).unwrap();
        assert_eq!((hit.start, hit.overlap), (0, 12));
    }

    #[test]
    fn reads_shorter_than_minimum_overlap_never_match() {
        assert!(find(&[ILLUMINA], default_params(), b"AGATCGG").is_none());
        assert!(find(&[ILLUMINA], default_params(), b"").is_none());
    }

    #[test]
    fn n_never_matches_in_read_or_adapter() {
        // A read with N at the adapter start: the N costs one error, which
        // exceeds the budget of floor(8 * 0.1) == 0.
        let mut tail = ILLUMINA[..8].to_vec();
        tail[0] = b'N';
        let read = [b"TTTTTTTTTTTT".as_slice(), &tail].concat();
        assert!(find(&[ILLUMINA], default_params(), &read).is_none());

        // With a wider budget the match is accepted and the N counted.
        let lenient = AdapterParams {
            max_error_rate: 0.5,
            ..default_params()
        };
        assert_eq!(find(&[ILLUMINA], lenient, &read).unwrap().errors, 1);

        // N in the adapter behaves identically, including N-against-N.
        let n_adapter = Adapter::new("n", b"NNNNNNNN").unwrap();
        assert!(find_three_prime(&[n_adapter], lenient, b"AAAANNNNNNNN").is_none());
    }

    #[test]
    fn stage_rejects_adapters_shorter_than_minimum_overlap() {
        let params = AdapterParams {
            min_overlap: 8,
            ..default_params()
        };
        let err = AdapterStage::new(adapters(&[b"ACGTACG"]), Vec::new(), params, None).unwrap_err();
        assert!(
            format!("{err}").contains("shorter than --min-overlap"),
            "{err}"
        );
        assert!(AdapterStage::new(Vec::new(), Vec::new(), params, None).is_err());
    }

    #[test]
    fn stage_rejects_invalid_thresholds() {
        let base = default_params();
        assert!(
            AdapterStage::new(
                adapters(&[ILLUMINA]),
                Vec::new(),
                AdapterParams {
                    min_overlap: 0,
                    ..base
                },
                None
            )
            .is_err()
        );
        assert!(
            AdapterStage::new(
                adapters(&[ILLUMINA]),
                Vec::new(),
                AdapterParams {
                    max_error_rate: 1.5,
                    ..base
                },
                None
            )
            .is_err()
        );
        assert!(
            AdapterStage::new(
                adapters(&[ILLUMINA]),
                Vec::new(),
                AdapterParams {
                    max_error_rate: f64::NAN,
                    ..base
                },
                None
            )
            .is_err()
        );
    }

    #[test]
    fn fasta_parsing_handles_multiline_records() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("adapters.fa");
        fs::write(&path, ">first adapter one\nACGT\nACGT\n>second\nagatcgga\n").unwrap();
        let adapters = read_adapter_fasta(&path).unwrap();
        assert_eq!(adapters.len(), 2);
        assert_eq!(adapters[0].name, "first");
        assert_eq!(adapters[0].sequence, b"ACGTACGT");
        assert_eq!(adapters[1].name, "second");
        assert_eq!(adapters[1].sequence, b"AGATCGGA");
    }

    #[test]
    fn fasta_parsing_rejects_malformed_files() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("bad.fa");
        fs::write(&path, "ACGT\n").unwrap();
        assert!(read_adapter_fasta(&path).is_err());
        fs::write(&path, "").unwrap();
        assert!(read_adapter_fasta(&path).is_err());
        fs::write(&path, ">a\nACGTZ\n").unwrap();
        assert!(read_adapter_fasta(&path).is_err());
    }

    // ------------------------------------------------------------ indel matching

    fn indel_params() -> AdapterParams {
        AdapterParams {
            allow_indels: true,
            ..AdapterParams::default()
        }
    }

    #[test]
    fn indels_match_an_insertion_in_the_read() {
        // The read carries the full adapter with one extra base inserted.
        let mut tail = ILLUMINA.to_vec();
        tail.insert(10, b'C');
        let read = [b"TTGGCCAATTGGCCAATTGG".as_slice(), &tail].concat();
        assert!(find(&[ILLUMINA], default_params(), &read).is_none());
        let hit = find(&[ILLUMINA], indel_params(), &read).unwrap();
        assert_eq!(
            (hit.start, hit.errors, hit.overlap),
            (20, 1, ILLUMINA.len())
        );
    }

    #[test]
    fn indels_match_a_deletion_in_the_read() {
        // The read carries the adapter with one base missing.
        let mut tail = ILLUMINA.to_vec();
        tail.remove(10);
        let read = [b"TTGGCCAATTGGCCAATTGG".as_slice(), &tail].concat();
        assert!(find(&[ILLUMINA], default_params(), &read).is_none());
        let hit = find(&[ILLUMINA], indel_params(), &read).unwrap();
        assert_eq!(
            (hit.start, hit.errors, hit.overlap),
            (20, 1, ILLUMINA.len())
        );
    }

    #[test]
    fn indel_matches_respect_the_error_budget() {
        // Five scattered insertions cost five errors: over the default budget
        // of floor(34 * 0.1) == 3.
        let mut tail = ILLUMINA.to_vec();
        for (n, position) in [4, 10, 16, 22, 28].into_iter().enumerate() {
            tail.insert(position + n, b'C');
        }
        let read = [b"TTGGCCAATTGGCCAATTGG".as_slice(), &tail].concat();
        assert!(find(&[ILLUMINA], indel_params(), &read).is_none());

        let lenient = AdapterParams {
            max_error_rate: 0.2,
            ..indel_params()
        };
        let hit = find(&[ILLUMINA], lenient, &read).unwrap();
        assert_eq!(hit.start, 20);
        assert!(hit.errors <= lenient.error_limit(hit.overlap));
    }

    #[test]
    fn an_earlier_indel_hit_beats_a_later_substitution_hit() {
        // One adapter copy with an insertion at offset 12, an exact copy
        // further downstream. The substitution scan can only recover the
        // exact copy; the indel scan finds the earlier contaminated one.
        let mut early = ILLUMINA.to_vec();
        early.insert(10, b'C');
        let read = [
            b"TTGGCCAATTGG".as_slice(),
            &early,
            b"AATTGGCCAATTGGCCAA".as_slice(),
            ILLUMINA,
        ]
        .concat();
        let substitution = find(&[ILLUMINA], default_params(), &read).unwrap();
        assert_eq!(substitution.start, 12 + early.len() + 18);
        let hit = find(&[ILLUMINA], indel_params(), &read).unwrap();
        assert_eq!(
            (hit.start, hit.errors),
            (12, 1),
            "the earliest coordinate wins"
        );
    }

    #[test]
    fn indel_shadows_never_over_trim_a_clean_match() {
        // Regression: with indels enabled, an exact adapter preceded by random
        // bases used to be trimmed one or two bases early, because a shifted
        // alignment is reachable by paying insertions and the earliest-
        // coordinate rule preferred it. Every insert prefix must trim at the
        // true adapter start.
        let mut state = 0x1234_5678_9abc_def1u64;
        let mut random = move || {
            state ^= state >> 12;
            state ^= state << 25;
            state ^= state >> 27;
            b"ACGT"[(state.wrapping_mul(0x2545_F491_4F6C_DD1D) >> 62) as usize]
        };
        for insert_length in 8..60 {
            let insert: Vec<u8> = (0..insert_length).map(|_| random()).collect();
            let read = [insert.as_slice(), ILLUMINA].concat();
            let hit = find(&[ILLUMINA], indel_params(), &read).unwrap();
            assert_eq!(
                (hit.start, hit.errors),
                (insert_length, 0),
                "insert of {insert_length} bases: {}",
                String::from_utf8_lossy(&read)
            );
        }
    }

    #[test]
    fn a_second_adapter_copy_does_not_rescue_the_first() {
        // Two occurrences far apart: the earliest still wins even though the
        // later one is cleaner. This is the counterpart of the shadow rule —
        // the window must not swallow a genuinely separate occurrence.
        let mut noisy = ILLUMINA.to_vec();
        noisy[20] = if noisy[20] == b'A' { b'C' } else { b'A' };
        noisy.insert(10, b'C');
        let read = [
            b"TTGGCCAATTGG".as_slice(),
            &noisy,
            b"CCAATTGGCCAA".as_slice(),
            ILLUMINA,
        ]
        .concat();
        let hit = find(&[ILLUMINA], indel_params(), &read).unwrap();
        assert_eq!(hit.start, 12, "the earliest occurrence is trimmed");
        assert!(hit.errors >= 2, "the earliest occurrence is the noisy one");
    }

    #[test]
    fn exact_matches_are_identical_with_and_without_indels() {
        let read = [b"ACGTACGTACGTACGTACGT".as_slice(), ILLUMINA].concat();
        let plain = find(&[ILLUMINA], default_params(), &read).unwrap();
        let indel = find(&[ILLUMINA], indel_params(), &read).unwrap();
        assert_eq!(
            (plain.start, plain.errors, plain.overlap),
            (20, 0, ILLUMINA.len())
        );
        assert_eq!(
            (indel.start, indel.errors, indel.overlap),
            (20, 0, ILLUMINA.len())
        );
    }

    #[test]
    fn a_short_exact_prefix_inside_a_longer_read_is_not_a_hit() {
        // The substitution contract requires the read to end in the adapter;
        // an exact 8-mer of adapter prefix followed by unrelated bases must
        // not match even with indels enabled.
        let read = [
            b"TTGGCCAATTGG".as_slice(),
            &ILLUMINA[..8],
            b"CCAATTGGCCAATTGGCCAATTGG".as_slice(),
        ]
        .concat();
        assert!(find(&[ILLUMINA], indel_params(), &read).is_none());
    }

    #[test]
    fn indels_still_trim_a_partial_adapter_at_the_read_end() {
        // Only the first 12 adapter bases are present, one of them inserted.
        let mut tail = ILLUMINA[..12].to_vec();
        tail.insert(5, b'C');
        let read = [b"TTGGCCAATTGG".as_slice(), &tail].concat();
        let hit = find(&[ILLUMINA], indel_params(), &read).unwrap();
        assert_eq!((hit.start, hit.errors, hit.overlap), (12, 1, 12));
    }

    #[test]
    fn indels_choose_the_best_adapter_at_a_coordinate() {
        // Two adapters: the second matches with one insertion, the first only
        // with several substitutions.
        let mut other = ILLUMINA.to_vec();
        for (i, base) in other.iter_mut().enumerate().take(8) {
            if i % 2 == 0 {
                *base = if *base == b'A' { b'C' } else { b'A' };
            }
        }
        let mut tail = ILLUMINA.to_vec();
        tail.insert(10, b'C');
        let read = [b"TTGGCCAATTGG".as_slice(), &tail].concat();
        let hit = find(&[&other, ILLUMINA], indel_params(), &read).unwrap();
        assert_eq!((hit.adapter, hit.errors), (1, 1));
    }

    #[test]
    fn indel_hits_are_bounded_by_the_substitution_hit_coordinate() {
        // With a clean exact match present, the indel scan must not roam past
        // it looking for something better: the exact hit wins on errors.
        let mut noisy = ILLUMINA.to_vec();
        noisy[3] = if noisy[3] == b'A' { b'C' } else { b'A' };
        let read = [
            b"TTGGCCAATTGG".as_slice(),
            &noisy,
            b"CCAATTGG".as_slice(),
            ILLUMINA,
        ]
        .concat();
        let hit = find(&[ILLUMINA], indel_params(), &read).unwrap();
        assert_eq!(hit.start, 12, "the earliest acceptable coordinate wins");
    }

    #[test]
    fn indel_matching_handles_degenerate_reads() {
        assert!(find(&[ILLUMINA], indel_params(), b"").is_none());
        assert!(find(&[ILLUMINA], indel_params(), b"AGATCGG").is_none());
        // A read that is exactly a noisy adapter prefix.
        let hit = find(&[ILLUMINA], indel_params(), &ILLUMINA[..12]);
        assert_eq!(hit.map(|h| (h.start, h.errors)), Some((0, 0)));
    }
}
