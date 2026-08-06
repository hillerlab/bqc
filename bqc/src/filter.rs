// Copyright (c) 2026 The Hiller Lab at the Senckenberg Gesellschaft für Naturforschung
// Distributed under the terms of the GNU General Public License, Version 3.0.

//! Per-read acceptance decisions.
//!
//! Filters are evaluated against the *transformed* read (the retained span
//! after adapter removal and trimming). Every applicable predicate is
//! evaluated, so a rejected read reports all of its failure reasons rather
//! than only the first.

use serde::Serialize;

use crate::error::{Error, Result};
use crate::read::{phred, phred_sum};

bitflags::bitflags! {
    /// The set of filter predicates a read failed.
    #[derive(Debug, Clone, Copy, PartialEq, Eq, Default, serde::Serialize)]
    pub struct FilterReason: u32 {
        const TOO_SHORT         = 1 << 0;
        const TOO_LONG          = 1 << 1;
        const TOO_MANY_N        = 1 << 2;
        const TOO_MANY_LOW_QUAL = 1 << 3;
        const LOW_MEAN_QUAL     = 1 << 4;
        const LOW_COMPLEXITY    = 1 << 5;
        /// No linked adapter definition matched, under `--linked-unmatched fail`.
        const LINKED_UNMATCHED  = 1 << 6;
    }
}

/// Every reason in stable reporting order.
pub const REASONS: [(FilterReason, &str); 7] = [
    (FilterReason::TOO_SHORT, "TOO_SHORT"),
    (FilterReason::TOO_LONG, "TOO_LONG"),
    (FilterReason::TOO_MANY_N, "TOO_MANY_N"),
    (FilterReason::TOO_MANY_LOW_QUAL, "TOO_MANY_LOW_QUAL"),
    (FilterReason::LOW_MEAN_QUAL, "LOW_MEAN_QUAL"),
    (FilterReason::LOW_COMPLEXITY, "LOW_COMPLEXITY"),
    (FilterReason::LINKED_UNMATCHED, "LINKED_UNMATCHED"),
];

impl FilterReason {
    /// Whether the read passed every filter.
    #[must_use]
    pub fn passed(self) -> bool {
        self.is_empty()
    }

    /// Renders the reasons as `A/B/C`, or `PASS` when empty.
    #[must_use]
    pub fn label(self) -> String {
        let names: Vec<&str> = REASONS
            .into_iter()
            .filter(|&(reason, _)| self.contains(reason))
            .map(|(_, name)| name)
            .collect();
        if names.is_empty() {
            "PASS".to_string()
        } else {
            names.join("/")
        }
    }
}

/// Compiled filter configuration.
#[derive(Debug, Clone, Copy, Default, Serialize)]
pub struct FilterStage {
    pub min_length: Option<usize>,
    pub max_length: Option<usize>,
    pub max_n: Option<usize>,
    pub max_n_fraction: Option<f64>,
    pub qualified_quality: u8,
    pub max_unqualified_bases: Option<usize>,
    pub max_unqualified_fraction: Option<f64>,
    pub min_mean_quality: Option<u8>,
    pub min_complexity: Option<f64>,
}

impl FilterStage {
    /// Whether no predicate is configured.
    #[must_use]
    pub fn is_noop(&self) -> bool {
        self.min_length.is_none()
            && self.max_length.is_none()
            && self.max_n.is_none()
            && self.max_n_fraction.is_none()
            && self.max_unqualified_bases.is_none()
            && self.max_unqualified_fraction.is_none()
            && self.min_mean_quality.is_none()
            && self.min_complexity.is_none()
    }

    /// Whether any configured predicate reads quality values.
    #[must_use]
    pub fn needs_quality(&self) -> bool {
        self.max_unqualified_bases.is_some()
            || self.max_unqualified_fraction.is_some()
            || self.min_mean_quality.is_some()
    }

    /// Validates thresholds and rejects no-op or contradictory settings.
    pub fn validate(self) -> Result<Self> {
        if self.is_noop() {
            return Err(Error::config(
                "filtering requires at least one filter \
                 (--min-length, --length-limit, --max-n, --max-n-fraction, \
                 --max-unqualified-bases, --max-unqualified-fraction, \
                 --min-mean-quality or --min-complexity)",
            ));
        }
        check_fraction("--max-n-fraction", self.max_n_fraction)?;
        check_fraction("--max-unqualified-fraction", self.max_unqualified_fraction)?;
        check_fraction("--min-complexity", self.min_complexity)?;
        if let (Some(min), Some(max)) = (self.min_length, self.max_length)
            && min > max
        {
            return Err(Error::config(format!(
                "--min-length ({min}) exceeds --length-limit ({max}); no read can pass"
            )));
        }
        Ok(self)
    }

    /// Evaluates every configured predicate against one transformed mate.
    #[must_use]
    pub fn evaluate(&self, sequence: &[u8], quality: Option<&[u8]>) -> FilterReason {
        let mut reasons = FilterReason::empty();
        let length = sequence.len();

        if self.min_length.is_some_and(|min| length < min) {
            reasons |= FilterReason::TOO_SHORT;
        }
        if self.max_length.is_some_and(|max| length > max) {
            reasons |= FilterReason::TOO_LONG;
        }

        if self.max_n.is_some() || self.max_n_fraction.is_some() {
            // A scalar count is fast enough for short reads; revisit only if
            // profiling shows base counting rather than (de)compression is hot.
            #[allow(clippy::naive_bytecount)]
            let n_count = sequence.iter().filter(|&&b| b == b'N').count();
            let too_many = self.max_n.is_some_and(|max| n_count > max)
                || self
                    .max_n_fraction
                    .is_some_and(|max| length > 0 && n_count as f64 / length as f64 > max);
            if too_many {
                reasons |= FilterReason::TOO_MANY_N;
            }
        }

        if let Some(quality) = quality {
            if self.max_unqualified_bases.is_some() || self.max_unqualified_fraction.is_some() {
                let threshold = self.qualified_quality;
                let unqualified = quality.iter().filter(|&&b| phred(b) < threshold).count();
                let too_many = self
                    .max_unqualified_bases
                    .is_some_and(|max| unqualified > max)
                    || self
                        .max_unqualified_fraction
                        .is_some_and(|max| length > 0 && unqualified as f64 / length as f64 > max);
                if too_many {
                    reasons |= FilterReason::TOO_MANY_LOW_QUAL;
                }
            }
            // Integer comparison of `sum < threshold * len` rather than a
            // floating-point mean. A zero-length read has a zero-valued sum
            // and a zero-valued bound, so it is not rejected here; use
            // --min-length to reject empty reads.
            if let Some(minimum) = self.min_mean_quality
                && phred_sum(quality) < u32::from(minimum) * length as u32
            {
                reasons |= FilterReason::LOW_MEAN_QUAL;
            }
        }

        if let Some(minimum) = self.min_complexity
            && complexity(sequence) < minimum
        {
            reasons |= FilterReason::LOW_COMPLEXITY;
        }

        reasons
    }
}

fn check_fraction(flag: &str, value: Option<f64>) -> Result<()> {
    if let Some(value) = value
        && !(0.0..=1.0).contains(&value)
    {
        return Err(Error::config(format!(
            "{flag} must be within 0.0..=1.0 (got {value})"
        )));
    }
    Ok(())
}

/// Fraction of adjacent base pairs that differ.
///
/// This is the lightweight complexity metric popularized by fastp. It detects
/// homopolymer-like reads and is **not** an entropy measure. Reads shorter than
/// two bases have no adjacent pairs and are defined to have complexity `0.0`.
#[must_use]
pub fn complexity(sequence: &[u8]) -> f64 {
    if sequence.len() < 2 {
        return 0.0;
    }
    let changes = sequence.windows(2).filter(|w| w[0] != w[1]).count();
    changes as f64 / (sequence.len() - 1) as f64
}

#[cfg(test)]
#[allow(clippy::float_cmp)]
mod tests {
    use super::*;

    fn stage() -> FilterStage {
        FilterStage {
            qualified_quality: 15,
            ..FilterStage::default()
        }
    }

    #[test]
    fn reason_labels_are_stable_and_ordered() {
        assert_eq!(FilterReason::empty().label(), "PASS");
        assert_eq!(FilterReason::TOO_SHORT.label(), "TOO_SHORT");
        let combined =
            FilterReason::LOW_COMPLEXITY | FilterReason::TOO_MANY_N | FilterReason::TOO_SHORT;
        assert_eq!(combined.label(), "TOO_SHORT/TOO_MANY_N/LOW_COMPLEXITY");
    }

    #[test]
    fn min_length_uses_inclusive_boundary() {
        let filter = FilterStage {
            min_length: Some(4),
            ..stage()
        };
        assert!(filter.evaluate(b"ACGT", None).passed());
        assert_eq!(filter.evaluate(b"ACG", None), FilterReason::TOO_SHORT);
        assert_eq!(filter.evaluate(b"", None), FilterReason::TOO_SHORT);
    }

    #[test]
    fn max_length_rejects_only_longer_reads() {
        let filter = FilterStage {
            max_length: Some(4),
            ..stage()
        };
        assert!(filter.evaluate(b"ACGT", None).passed());
        assert_eq!(filter.evaluate(b"ACGTA", None), FilterReason::TOO_LONG);
    }

    #[test]
    fn max_n_count_and_fraction_are_independent() {
        let by_count = FilterStage {
            max_n: Some(2),
            ..stage()
        };
        assert!(by_count.evaluate(b"ANNGT", None).passed());
        assert_eq!(by_count.evaluate(b"ANNNT", None), FilterReason::TOO_MANY_N);

        // 2 of 5 bases == 0.4, which is not greater than 0.4.
        let by_fraction = FilterStage {
            max_n_fraction: Some(0.4),
            ..stage()
        };
        assert!(by_fraction.evaluate(b"ANNGT", None).passed());
        assert_eq!(
            by_fraction.evaluate(b"ANNNT", None),
            FilterReason::TOO_MANY_N
        );

        // A count limit alone tolerates a high fraction in short reads.
        let by_count = FilterStage {
            max_n: Some(2),
            ..stage()
        };
        assert!(by_count.evaluate(b"NN", None).passed());
        let by_fraction = FilterStage {
            max_n_fraction: Some(0.4),
            ..stage()
        };
        assert_eq!(by_fraction.evaluate(b"NN", None), FilterReason::TOO_MANY_N);
    }

    #[test]
    fn unqualified_base_limits_use_the_qualified_threshold() {
        // '0' == Q15, '/' == Q14.
        let filter = FilterStage {
            qualified_quality: 15,
            max_unqualified_bases: Some(1),
            ..stage()
        };
        assert!(filter.evaluate(b"ACGT", Some(b"000/")).passed());
        assert_eq!(
            filter.evaluate(b"ACGT", Some(b"00//")),
            FilterReason::TOO_MANY_LOW_QUAL
        );

        let filter = FilterStage {
            qualified_quality: 15,
            max_unqualified_fraction: Some(0.25),
            ..stage()
        };
        assert!(filter.evaluate(b"ACGT", Some(b"000/")).passed());
        assert_eq!(
            filter.evaluate(b"ACGT", Some(b"00//")),
            FilterReason::TOO_MANY_LOW_QUAL
        );
    }

    #[test]
    fn mean_quality_boundary_is_inclusive() {
        // '5' == Q20 exactly.
        let filter = FilterStage {
            min_mean_quality: Some(20),
            ..stage()
        };
        assert!(filter.evaluate(b"ACGT", Some(b"5555")).passed());
        assert_eq!(
            filter.evaluate(b"ACGT", Some(b"5554")),
            FilterReason::LOW_MEAN_QUAL
        );
    }

    #[test]
    fn complexity_definition_matches_documentation() {
        assert_eq!(complexity(b""), 0.0);
        assert_eq!(complexity(b"A"), 0.0, "reads shorter than 2 bases are 0.0");
        assert_eq!(complexity(b"AAAA"), 0.0);
        assert_eq!(complexity(b"ACGT"), 1.0);
        assert_eq!(complexity(b"ACAC"), 1.0);
        assert!((complexity(b"AACC") - 1.0 / 3.0).abs() < 1e-12);
    }

    #[test]
    fn complexity_filter_rejects_homopolymers() {
        let filter = FilterStage {
            min_complexity: Some(0.3),
            ..stage()
        };
        assert!(filter.evaluate(b"ACGTACGT", None).passed());
        assert_eq!(
            filter.evaluate(b"AAAAAAAA", None),
            FilterReason::LOW_COMPLEXITY
        );
        assert_eq!(filter.evaluate(b"A", None), FilterReason::LOW_COMPLEXITY);
    }

    #[test]
    fn every_applicable_reason_is_collected() {
        let filter = FilterStage {
            min_length: Some(10),
            max_n: Some(0),
            qualified_quality: 20,
            max_unqualified_bases: Some(0),
            min_mean_quality: Some(30),
            min_complexity: Some(0.5),
            ..stage()
        };
        let reasons = filter.evaluate(b"NNNN", Some(b"!!!!"));
        assert_eq!(
            reasons,
            FilterReason::TOO_SHORT
                | FilterReason::TOO_MANY_N
                | FilterReason::TOO_MANY_LOW_QUAL
                | FilterReason::LOW_MEAN_QUAL
                | FilterReason::LOW_COMPLEXITY
        );
        assert!(!reasons.passed());
    }

    #[test]
    fn quality_filters_are_skipped_when_quality_is_absent() {
        let filter = FilterStage {
            min_mean_quality: Some(40),
            max_unqualified_bases: Some(0),
            ..stage()
        };
        assert!(filter.evaluate(b"ACGT", None).passed());
        assert!(filter.needs_quality());
    }

    #[test]
    fn validation_rejects_empty_and_contradictory_configurations() {
        assert!(stage().validate().is_err());
        assert!(
            FilterStage {
                min_length: Some(1),
                ..stage()
            }
            .validate()
            .is_ok()
        );
        assert!(
            FilterStage {
                min_length: Some(50),
                max_length: Some(10),
                ..stage()
            }
            .validate()
            .is_err()
        );
        assert!(
            FilterStage {
                max_n_fraction: Some(1.5),
                ..stage()
            }
            .validate()
            .is_err()
        );
        assert!(
            FilterStage {
                min_complexity: Some(-0.1),
                ..stage()
            }
            .validate()
            .is_err()
        );
        assert!(
            FilterStage {
                min_complexity: Some(f64::NAN),
                ..stage()
            }
            .validate()
            .is_err()
        );
    }
}
