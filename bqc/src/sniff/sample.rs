// Copyright (c) 2026 The Hiller Lab at the Senckenberg Gesellschaft für Naturforschung
// Distributed under the terms of the GNU General Public License, Version 3.0.

//! Deterministic distributed sampling.
//!
//! Both `sniff` subcommands inspect a subset of the input rather than all of it,
//! and both must return the same answer on every run and at every thread count.
//!
//! Scanning the leading `N` records — what a FASTQ tool has to do — is biased by
//! everything that puts unusual data at the start of a file: concatenated lanes,
//! tile-ordered records, joined runs, calibration reads. CBQ is indexed and
//! block-addressable, so `bqc` spreads the sample evenly across the selected
//! range instead.
//!
//! The selection is arithmetic, not random: for `S` records drawn from `N`,
//!
//! ```text
//! index(i) = (2i + 1) * N / (2S)      for i in 0..S
//! ```
//!
//! which is `floor((i + 0.5) * N / S)` evaluated in integers. There is no
//! generator and no seed, so nothing has to be recorded to reproduce a run. The
//! `u128` intermediate keeps it exact for any record count a file can hold,
//! where the floating-point form would depend on rounding.
//!
//! The indices are strictly increasing (consecutive values differ by at least
//! `N/S >= 1`), all lie in `0..N`, and never repeat.

use std::sync::atomic::{AtomicUsize, Ordering};

use binseq::cbq::ColumnarBlock;
use binseq::prelude::BinseqRecord;
use zstd::zstd_safe;

use crate::error::Result;
use crate::io::CbqInput;

/// A deterministic, evenly distributed selection of record indices.
///
/// Indices are generated on demand rather than stored: a one-million-record
/// sample would otherwise hold 8 MB for numbers that are cheaper to recompute.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize)]
pub struct SamplePlan {
    /// First record index of the selected range.
    pub range_start: u64,
    /// One past the last record index of the selected range.
    pub range_end: u64,
    /// Records the caller asked for.
    pub requested: u64,
    /// Records actually selected: `min(requested, range length)`.
    pub selected: u64,
}

impl SamplePlan {
    /// Plans a sample of `requested` records from `range_start..range_end`.
    #[must_use]
    pub fn new(range_start: u64, range_end: u64, requested: u64) -> Self {
        let range_end = range_end.max(range_start);
        let available = range_end - range_start;
        Self {
            range_start,
            range_end,
            requested,
            selected: requested.min(available),
        }
    }

    /// Plans a sample over a whole input, honouring an optional `--span`.
    #[must_use]
    pub fn for_input(
        input: &CbqInput,
        span: Option<&std::ops::Range<u64>>,
        requested: u64,
    ) -> Self {
        let total = input.num_records();
        let (start, end) = match span {
            Some(span) => (span.start.min(total), span.end.min(total)),
            None => (0, total),
        };
        Self::new(start, end, requested)
    }

    /// Number of records the range contains.
    #[must_use]
    pub fn available(&self) -> u64 {
        self.range_end - self.range_start
    }

    /// The absolute record index of the `i`th sampled record.
    ///
    /// # Panics
    ///
    /// Panics when `i >= self.selected`.
    #[must_use]
    pub fn index(&self, i: u64) -> u64 {
        assert!(i < self.selected, "sample ordinal out of range");
        let n = u128::from(self.available());
        let s = u128::from(self.selected);
        let offset = ((2 * u128::from(i) + 1) * n) / (2 * s);
        self.range_start + offset as u64
    }

    /// Every selected record index, ascending.
    pub fn indices(&self) -> impl Iterator<Item = u64> + '_ {
        (0..self.selected).map(|i| self.index(i))
    }
}

/// The sampled records that live in one CBQ block.
///
/// `ordinals` indexes [`SamplePlan`], not the file. Because both the plan's
/// indices and the file's blocks are ascending, each block owns a contiguous
/// run of ordinals, so one work item describes it exactly.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BlockWork {
    /// Index into [`CbqInput::blocks`].
    pub block: usize,
    /// First plan ordinal in this block.
    pub first_ordinal: u64,
    /// One past the last plan ordinal in this block.
    pub end_ordinal: u64,
}

/// Groups a plan's records into the blocks that hold them.
///
/// Blocks containing no sampled record are absent, so a worker never
/// decompresses data the plan does not touch, and every listed block is
/// decompressed exactly once per pass.
#[must_use]
pub fn block_work(plan: &SamplePlan, input: &CbqInput) -> Vec<BlockWork> {
    let blocks = input.blocks();
    let mut work: Vec<BlockWork> = Vec::new();
    let mut cursor = 0usize;
    for ordinal in 0..plan.selected {
        let record = plan.index(ordinal);
        // Both sequences ascend, so the block cursor only moves forwards.
        while cursor < blocks.len() && blocks[cursor].end_record() <= record {
            cursor += 1;
        }
        if cursor >= blocks.len() {
            break;
        }
        match work.last_mut() {
            Some(last) if last.block == cursor => last.end_ordinal = ordinal + 1,
            _ => work.push(BlockWork {
                block: cursor,
                first_ordinal: ordinal,
                end_ordinal: ordinal + 1,
            }),
        }
    }
    work
}

/// Runs `observe` over every sampled record, in parallel over blocks.
///
/// Each worker owns its state; the caller merges. The merge must be
/// commutative — integer addition — because which worker claims which block
/// depends on scheduling.
pub(crate) fn for_each_sampled<S, M, F>(
    input: &CbqInput,
    work: &[BlockWork],
    plan: &SamplePlan,
    threads: usize,
    make: M,
    observe: F,
) -> Result<Vec<S>>
where
    S: Send,
    M: Fn() -> S + Sync,
    F: Fn(&mut S, &[u8], Option<&[u8]>) + Sync,
{
    let paired = input.schema().paired;
    let run_one = |state: &mut S, block: &ColumnarBlock, item: &BlockWork, range| {
        let mut ordinal = item.first_ordinal;
        for record in block.iter_records(range) {
            if ordinal >= item.end_ordinal {
                break;
            }
            // Both the record stream and the plan ascend, so one cursor walks
            // the intersection without searching.
            if record.index() != plan.index(ordinal) {
                continue;
            }
            ordinal += 1;
            // Both mates arrive together so a pair can be analysed as a pair —
            // the overlap between them is inferred once, not once per mate.
            observe(state, record.sseq(), paired.then(|| record.xseq()));
        }
    };

    if threads <= 1 || work.len() <= 1 {
        let mut state = make();
        let mut block = ColumnarBlock::new(input.header());
        let mut dctx = zstd_safe::DCtx::create();
        for item in work {
            let range = input.load(item.block, &mut block, &mut dctx)?;
            run_one(&mut state, &block, item, range);
        }
        return Ok(vec![state]);
    }

    let next = AtomicUsize::new(0);
    let mut states: Vec<Result<S>> = Vec::new();
    std::thread::scope(|scope| {
        let mut handles = Vec::with_capacity(threads);
        for _ in 0..threads.min(work.len()) {
            let next = &next;
            let make = &make;
            let run_one = &run_one;
            handles.push(scope.spawn(move || {
                let mut state = make();
                let mut block = ColumnarBlock::new(input.header());
                let mut dctx = zstd_safe::DCtx::create();
                loop {
                    let ordinal = next.fetch_add(1, Ordering::Relaxed);
                    let Some(item) = work.get(ordinal) else { break };
                    let range = input.load(item.block, &mut block, &mut dctx)?;
                    run_one(&mut state, &block, item, range);
                }
                Ok(state)
            }));
        }
        for handle in handles {
            states.push(
                handle
                    .join()
                    .unwrap_or(Err(crate::error::Error::WorkerPanic)),
            );
        }
    });
    states.into_iter().collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_sample_larger_than_the_input_selects_every_record() {
        let plan = SamplePlan::new(0, 10, 100);
        assert_eq!(plan.selected, 10);
        assert_eq!(
            plan.indices().collect::<Vec<_>>(),
            (0..10).collect::<Vec<_>>()
        );
    }

    #[test]
    fn a_sample_equal_to_the_input_selects_every_record() {
        let plan = SamplePlan::new(0, 7, 7);
        assert_eq!(
            plan.indices().collect::<Vec<_>>(),
            (0..7).collect::<Vec<_>>()
        );
    }

    #[test]
    fn an_empty_input_selects_nothing() {
        let plan = SamplePlan::new(0, 0, 100);
        assert_eq!(plan.selected, 0);
        assert_eq!(plan.indices().count(), 0);
    }

    #[test]
    fn one_record_is_selectable() {
        let plan = SamplePlan::new(0, 1, 1);
        assert_eq!(plan.indices().collect::<Vec<_>>(), vec![0]);
    }

    #[test]
    fn indices_are_sorted_distinct_and_in_range() {
        // Deliberately non-divisible combinations.
        for (n, s) in [
            (1000u64, 7u64),
            (7, 1000),
            (1_000_000, 262_144),
            (13, 5),
            (5, 13),
            (u64::from(u32::MAX), 1000),
        ] {
            let plan = SamplePlan::new(0, n, s);
            let indices: Vec<u64> = plan.indices().collect();
            assert_eq!(indices.len() as u64, plan.selected, "n={n} s={s}");
            for pair in indices.windows(2) {
                assert!(pair[0] < pair[1], "not strictly increasing: n={n} s={s}");
            }
            if let Some(&last) = indices.last() {
                assert!(last < n, "index past the end: n={n} s={s}");
            }
        }
    }

    #[test]
    fn the_sample_is_spread_across_the_whole_range() {
        // A leading-prefix sampler would put every index below 1000.
        let plan = SamplePlan::new(0, 1_000_000, 1000);
        let indices: Vec<u64> = plan.indices().collect();
        assert!(
            indices[0] < 1000,
            "first index {} is not near the start",
            indices[0]
        );
        assert!(
            *indices.last().unwrap() > 999_000,
            "last index {} is not near the end",
            indices.last().unwrap()
        );
        // Even coverage: every decile holds a tenth of the sample.
        for decile in 0..10u64 {
            let lo = decile * 100_000;
            let hi = lo + 100_000;
            let count = indices.iter().filter(|&&i| i >= lo && i < hi).count();
            assert_eq!(count, 100, "decile {decile} holds {count} of 1000");
        }
    }

    #[test]
    fn a_span_restricts_the_sample_to_that_range() {
        let plan = SamplePlan::new(500, 1500, 10);
        let indices: Vec<u64> = plan.indices().collect();
        assert_eq!(indices.len(), 10);
        assert!(indices.iter().all(|&i| (500..1500).contains(&i)));
        assert_eq!(indices[0], 550);
    }

    #[test]
    fn very_large_record_counts_stay_exact() {
        // The f64 form of this expression loses integer precision above 2^53.
        let n = 1u64 << 60;
        let plan = SamplePlan::new(0, n, 4);
        let indices: Vec<u64> = plan.indices().collect();
        assert_eq!(indices, vec![n / 8, 3 * (n / 8), 5 * (n / 8), 7 * (n / 8)]);
    }

    #[test]
    fn the_plan_does_not_depend_on_how_it_is_evaluated() {
        let plan = SamplePlan::new(3, 9999, 321);
        let once: Vec<u64> = plan.indices().collect();
        let again: Vec<u64> = (0..plan.selected).map(|i| plan.index(i)).collect();
        assert_eq!(once, again);
    }
}
