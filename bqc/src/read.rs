// Copyright (c) 2026 The Hiller Lab at the Senckenberg Gesellschaft für Naturforschung
// Distributed under the terms of the GNU General Public License, Version 3.0.

//! Borrowed read views and retained coordinate spans.
//!
//! Every transformation in `bqc` is expressed as an update to a [`Span`]
//! over the original, borrowed record. No transformation allocates or rewrites
//! sequence or quality data.

use crate::error::{Error, Result};

/// Offset of the Phred+33 (Sanger / Illumina 1.8+) quality encoding.
pub const PHRED_OFFSET: u8 = 33;

/// Highest quality byte accepted by the Phred+33 encoding.
pub const PHRED_MAX_BYTE: u8 = 126;

/// Decodes one encoded quality byte into a Phred score.
///
/// This is the only place where quality bytes are interpreted. Bytes below the
/// encoding offset saturate to zero; callers that need to reject them use
/// [`validate_quality`].
#[inline]
#[must_use]
pub fn phred(byte: u8) -> u8 {
    byte.saturating_sub(PHRED_OFFSET)
}

/// Sums the Phred scores of `quality`.
#[inline]
#[must_use]
pub fn phred_sum(quality: &[u8]) -> u32 {
    quality.iter().map(|&b| u32::from(phred(b))).sum()
}

/// Verifies that every byte of `quality` is a valid Phred+33 value.
pub fn validate_quality(quality: &[u8], record: u64, mate: &'static str) -> Result<()> {
    for &byte in quality {
        if !(PHRED_OFFSET..=PHRED_MAX_BYTE).contains(&byte) {
            return Err(Error::InvalidQualityEncoding { record, mate, byte });
        }
    }
    Ok(())
}

/// A borrowed view of one mate of a record.
///
/// `quality` is `Some` only when the source CBQ file declares a quality
/// column. Records from quality-free files expose a synthetic quality buffer
/// through `binseq`, which must never be mistaken for real data.
#[derive(Debug, Clone, Copy)]
pub struct ReadView<'a> {
    pub sequence: &'a [u8],
    pub quality: Option<&'a [u8]>,
}

impl<'a> ReadView<'a> {
    /// Builds a view, validating the sequence/quality length invariant.
    pub fn new(
        sequence: &'a [u8],
        quality: Option<&'a [u8]>,
        record: u64,
        mate: &'static str,
    ) -> Result<Self> {
        if let Some(quality) = quality
            && quality.len() != sequence.len()
        {
            return Err(Error::SequenceQualityLengthMismatch {
                record,
                mate,
                sequence: sequence.len(),
                quality: quality.len(),
            });
        }
        Ok(Self { sequence, quality })
    }

    /// Builds a view without validation. Intended for tests and internal use
    /// where the invariant is already established.
    #[must_use]
    pub fn unchecked(sequence: &'a [u8], quality: Option<&'a [u8]>) -> Self {
        debug_assert!(quality.is_none_or(|q| q.len() == sequence.len()));
        Self { sequence, quality }
    }

    /// Number of bases in the underlying read.
    #[must_use]
    pub fn len(&self) -> usize {
        self.sequence.len()
    }

    /// Whether the underlying read has no bases. Required alongside [`Self::len`]
    /// by `clippy::len_without_is_empty`.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }
}

/// Half-open coordinate range retained from a read.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Span {
    pub start: usize,
    pub end: usize,
}

impl Span {
    /// The span covering an entire read of `length` bases.
    #[inline]
    #[must_use]
    pub fn full(length: usize) -> Self {
        Self {
            start: 0,
            end: length,
        }
    }

    /// Number of retained bases.
    #[inline]
    #[must_use]
    pub fn len(&self) -> usize {
        debug_assert!(self.start <= self.end);
        self.end - self.start
    }

    /// Whether the span retains no bases.
    #[inline]
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.start >= self.end
    }

    /// The retained sequence slice.
    #[inline]
    #[must_use]
    pub fn sequence<'a>(&self, read: ReadView<'a>) -> &'a [u8] {
        &read.sequence[self.start..self.end]
    }

    /// The retained quality slice, when the record carries qualities.
    #[inline]
    #[must_use]
    pub fn quality<'a>(&self, read: ReadView<'a>) -> Option<&'a [u8]> {
        read.quality.map(|q| &q[self.start..self.end])
    }

    /// Removes `bases` from the 5' end, clamping at the 3' boundary.
    #[inline]
    pub fn trim_front(&mut self, bases: usize) {
        self.start = self.start.saturating_add(bases).min(self.end);
    }

    /// Removes `bases` from the 3' end, clamping at the 5' boundary.
    #[inline]
    pub fn trim_back(&mut self, bases: usize) {
        self.end = self.end.saturating_sub(bases).max(self.start);
    }

    /// Retains at most `length` bases from the 5' end.
    #[inline]
    pub fn truncate_to(&mut self, length: usize) {
        if self.len() > length {
            self.end = self.start + length;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn phred_decoding_is_offset_33() {
        assert_eq!(phred(b'!'), 0);
        assert_eq!(phred(b'#'), 2);
        assert_eq!(phred(b'5'), 20);
        assert_eq!(phred(b'I'), 40);
        assert_eq!(phred(0), 0, "bytes below the offset saturate");
    }

    #[test]
    fn phred_sum_matches_manual_sum() {
        assert_eq!(phred_sum(b"IIII"), 160);
        assert_eq!(phred_sum(b""), 0);
    }

    #[test]
    fn quality_validation_rejects_out_of_range_bytes() {
        assert!(validate_quality(b"!~", 0, "R1").is_ok());
        let err = validate_quality(b"I\x7f", 7, "R2").unwrap_err();
        assert!(matches!(
            err,
            Error::InvalidQualityEncoding {
                record: 7,
                byte: 0x7f,
                ..
            }
        ));
        assert!(validate_quality(b"\x20", 0, "R1").is_err());
    }

    #[test]
    fn read_view_rejects_length_mismatch() {
        assert!(ReadView::new(b"ACGT", Some(b"IIII"), 0, "R1").is_ok());
        assert!(ReadView::new(b"ACGT", None, 0, "R1").is_ok());
        let err = ReadView::new(b"ACGT", Some(b"III"), 3, "R1").unwrap_err();
        assert!(matches!(
            err,
            Error::SequenceQualityLengthMismatch {
                record: 3,
                sequence: 4,
                quality: 3,
                ..
            }
        ));
    }

    #[test]
    fn span_slicing_is_consistent_across_sequence_and_quality() {
        let read = ReadView::unchecked(b"ACGTACGT", Some(b"IIIIFFFF"));
        let span = Span { start: 2, end: 6 };
        assert_eq!(span.sequence(read), b"GTAC");
        assert_eq!(span.quality(read), Some(b"IIFF".as_slice()));
        assert_eq!(span.len(), 4);
    }

    #[test]
    fn span_trimming_clamps_instead_of_wrapping() {
        let mut span = Span::full(10);
        span.trim_front(4);
        assert_eq!(span, Span { start: 4, end: 10 });
        span.trim_back(20);
        assert_eq!(
            span,
            Span { start: 4, end: 4 },
            "over-trim collapses to empty"
        );
        assert!(span.is_empty());

        let mut span = Span::full(10);
        span.trim_front(20);
        assert_eq!(span, Span { start: 10, end: 10 });

        let mut span = Span::full(10);
        span.truncate_to(4);
        assert_eq!(span, Span { start: 0, end: 4 });
        span.truncate_to(100);
        assert_eq!(span, Span { start: 0, end: 4 }, "truncation never extends");
    }
}
