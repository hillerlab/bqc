// Copyright (c) 2026 The Hiller Lab at the Senckenberg Gesellschaft für Naturforschung
// Distributed under the terms of the GNU General Public License, Version 3.0.

//! The order-preserving processing engine.
//!
//! One CBQ block is one unit of work. Workers decode, process and re-encode
//! whole blocks independently; a single committer appends the finished,
//! already-compressed blocks in input order. Because the partitioning comes
//! from the input file rather than from the thread count, results — and in
//! practice the output bytes themselves — do not depend on `--threads`.
//!
//! ```text
//! input blocks 0..N          bounded channel        ordered commit
//!      |                          |                       |
//!      +--> worker 0 --+          |                       |
//!      +--> worker 1 --+--> completed chunks --> BTreeMap --> CBQ writer
//!      +--> worker k --+
//! ```

use std::collections::BTreeMap;
use std::ops::Range;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::mpsc;
use std::sync::{Condvar, Mutex};

use binseq::BinseqRecord;
use binseq::cbq::{ColumnarBlock, FileHeader};
use zstd::zstd_safe;

use crate::error::{Error, Result};
use crate::io::{
    CbqInput, CbqOutput, Fragment, MateOutput, Schema, TextOutput, push_record, write_reason_row,
};
use crate::process::{FailedMode, Mate, Workflow};
use crate::read::{ReadView, Span};
use crate::stats::Stats;

/// Runtime options that are independent of the biological configuration.
#[derive(Debug, Clone)]
pub struct RunOptions {
    pub threads: usize,
    /// Original record indices to process, or `None` for the whole file.
    pub span: Option<Range<u64>>,
    pub failed_mode: FailedMode,
}

/// The destinations a run writes to.
pub struct Outputs<'a> {
    pub accepted: &'a mut CbqOutput,
    pub failed: Option<&'a mut CbqOutput>,
    pub reasons: Option<&'a mut TextOutput>,
    /// Single-end output for pairs where only R1 passed (orphan policy).
    pub orphan_r1: Option<&'a mut CbqOutput>,
    /// Single-end output for pairs where only R2 passed (orphan policy).
    pub orphan_r2: Option<&'a mut CbqOutput>,
    /// Correction log, written in record order.
    pub corrections: Option<&'a mut TextOutput>,
    /// Segmentation provenance sidecar, written in source-record order.
    pub segments: Option<&'a mut TextOutput>,
}

/// Everything one worker produces for one input block.
struct ChunkOutput {
    ordinal: usize,
    /// Correction log rows for this chunk, committed in chunk order.
    corrections: Vec<u8>,
    /// Segmentation provenance rows for this chunk.
    segments: Vec<u8>,
    accepted: Fragment,
    failed: Option<Fragment>,
    orphan_r1: Option<Fragment>,
    orphan_r2: Option<Fragment>,
    reasons: Vec<u8>,
    stats: Stats,
}

/// Processes `input` and writes the results in input order.
pub fn run(
    input: &CbqInput,
    workflow: &Workflow,
    outputs: &mut Outputs<'_>,
    options: &RunOptions,
) -> Result<Stats> {
    let schema = input.schema();
    if workflow.needs_quality() && !schema.quality {
        return Err(Error::MissingQuality("quality trimming or filtering"));
    }

    let span = resolve_span(input, options.span.clone())?;
    let chunks: Vec<usize> = input
        .blocks()
        .iter()
        .enumerate()
        .filter(|(_, block)| block.first_record < span.end && block.end_record() > span.start)
        .map(|(index, _)| index)
        .collect();

    if let Some(reasons) = outputs.reasons.as_mut() {
        reasons.write(crate::io::REASON_HEADER.as_bytes())?;
    }
    if let Some(segments) = outputs.segments.as_mut() {
        segments.write(crate::io::SEGMENT_HEADER.as_bytes())?;
    }
    if let Some(corrections) = outputs.corrections.as_mut() {
        let detail = workflow
            .correction
            .map_or(crate::correct::LogDetail::Reads, |stage| stage.log_detail);
        corrections.write(crate::io::correction_log_header(detail).as_bytes())?;
    }

    // Both orphan destinations share the same single-end header; either one
    // can provide it to the worker fragments.
    let orphan_header = outputs
        .orphan_r1
        .as_deref()
        .or(outputs.orphan_r2.as_deref())
        .map(CbqOutput::header);
    let context = Context {
        input,
        workflow,
        schema,
        span,
        failed_mode: options.failed_mode,
        output_header: outputs.accepted.header(),
        orphan_header,
        want_failed: outputs.failed.is_some(),
        want_reasons: outputs.reasons.is_some(),
        want_corrections: outputs.corrections.is_some(),
        want_segments: outputs.segments.is_some(),
        correction_detail: workflow
            .correction
            .map_or(crate::correct::LogDetail::Reads, |stage| stage.log_detail),
        want_orphan_r1: outputs.orphan_r1.is_some(),
        want_orphan_r2: outputs.orphan_r2.is_some(),
    };

    if options.threads <= 1 {
        run_sequential(&context, &chunks, outputs)
    } else {
        run_parallel(&context, &chunks, outputs, options.threads)
    }
}

/// Worker-local buffers, reused across every record a worker touches.
#[derive(Default)]
struct Scratch {
    correction: crate::correct::CorrectionScratch,
    segment: crate::segment::SegmentScratch,
    fragments: Vec<crate::process::FragmentResult>,
    /// The header of the fragment currently being written.
    header: Vec<u8>,
}

/// Immutable state shared by every worker.
struct Context<'a> {
    input: &'a CbqInput,
    workflow: &'a Workflow,
    schema: Schema,
    span: Range<u64>,
    failed_mode: FailedMode,
    output_header: FileHeader,
    orphan_header: Option<FileHeader>,
    want_failed: bool,
    want_reasons: bool,
    want_corrections: bool,
    want_segments: bool,
    correction_detail: crate::correct::LogDetail,
    want_orphan_r1: bool,
    want_orphan_r2: bool,
}

fn resolve_span(input: &CbqInput, span: Option<Range<u64>>) -> Result<Range<u64>> {
    let total = input.num_records();
    match span {
        None => Ok(0..total),
        Some(span) => {
            if span.start > span.end {
                return Err(Error::config(format!(
                    "--span start ({}) is greater than end ({})",
                    span.start, span.end
                )));
            }
            if span.start > total {
                return Err(Error::config(format!(
                    "--span start ({}) is past the last record ({total} records)",
                    span.start
                )));
            }
            Ok(span.start..span.end.min(total))
        }
    }
}

/// Turns a panic inside chunk processing into a controlled command failure.
///
/// The committer then stops, no output is renamed into place, and every
/// temporary file is removed on drop. The scratch state the closure borrowed is
/// abandoned rather than reused, which is why `AssertUnwindSafe` is sound here.
fn guarded<F: FnOnce() -> Result<ChunkOutput>>(process: F) -> Result<ChunkOutput> {
    std::panic::catch_unwind(std::panic::AssertUnwindSafe(process))
        .unwrap_or(Err(Error::WorkerPanic))
}

#[derive(Default)]
struct WindowState {
    next: usize,
    committed: usize,
}

impl WindowState {
    /// Claims the next ordinal only while it is within `limit` of commit.
    fn claim(&mut self, total: usize, limit: usize) -> Option<usize> {
        if self.next >= total || self.next >= self.committed.saturating_add(limit) {
            return None;
        }
        let ordinal = self.next;
        self.next += 1;
        Some(ordinal)
    }
}

/// Sliding claim window tied to ordered commit progress.
///
/// A bounded result channel alone does not bound the reorder map: the committer
/// can keep draining later chunks while an early chunk is slow. This gate caps
/// all claimed work — processing, queued and pending combined — relative to the
/// number of chunks that have actually been committed.
struct CommitWindow {
    limit: usize,
    state: Mutex<WindowState>,
    advanced: Condvar,
}

impl CommitWindow {
    fn new(limit: usize) -> Self {
        Self {
            limit: limit.max(1),
            state: Mutex::new(WindowState::default()),
            advanced: Condvar::new(),
        }
    }

    fn claim(&self, total: usize, abort: &AtomicBool) -> Option<usize> {
        let mut state = self
            .state
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        loop {
            if abort.load(Ordering::Relaxed) || state.next >= total {
                return None;
            }
            if let Some(ordinal) = state.claim(total, self.limit) {
                return Some(ordinal);
            }
            state = self
                .advanced
                .wait(state)
                .unwrap_or_else(std::sync::PoisonError::into_inner);
        }
    }

    fn commit_through(&self, committed: usize) {
        let mut state = self
            .state
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        debug_assert!(committed >= state.committed);
        state.committed = committed;
        self.advanced.notify_all();
    }

    fn wake_all(&self) {
        let _guard = self
            .state
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        self.advanced.notify_all();
    }
}

fn new_stats(workflow: &Workflow) -> Stats {
    let (r1, r2) = match (workflow.adapter.as_ref(), workflow.segment.as_ref()) {
        (Some(stage), _) => (stage.r1.len(), stage.r2.len()),
        // Delimiters are counted in the same per-adapter table.
        (None, Some(stage)) => (stage.adapters.len(), 0),
        (None, None) => (0, 0),
    };
    let (linked_r1, linked_r2) = workflow
        .linked
        .as_ref()
        .map_or((0, 0), |stage| (stage.r1.len(), stage.r2.len()));
    Stats::new(r1, r2).with_linked(linked_r1, linked_r2)
}

fn run_sequential(
    context: &Context<'_>,
    chunks: &[usize],
    outputs: &mut Outputs<'_>,
) -> Result<Stats> {
    let mut stats = new_stats(context.workflow);
    let mut block = ColumnarBlock::new(context.input.header());
    let mut dctx = zstd_safe::DCtx::create();
    let mut scratch = Scratch::default();
    for (ordinal, &block_index) in chunks.iter().enumerate() {
        let mut chunk = guarded(|| {
            process_chunk(
                context,
                ordinal,
                block_index,
                &mut block,
                &mut dctx,
                &mut scratch,
            )
        })?;
        commit_chunk(&mut chunk, outputs, &mut stats)?;
    }
    Ok(stats)
}

fn run_parallel(
    context: &Context<'_>,
    chunks: &[usize],
    outputs: &mut Outputs<'_>,
    threads: usize,
) -> Result<Stats> {
    let abort = AtomicBool::new(false);
    // Each worker may be processing one chunk while up to two per worker wait
    // in the result path. Tying that total window to commit progress bounds the
    // reorder map even when the first chunk is arbitrarily slower than the rest.
    let in_flight = threads.saturating_mul(3).max(1);
    let window = CommitWindow::new(in_flight);
    let channel_capacity = threads.saturating_mul(2).max(1);
    let (sender, receiver) = mpsc::sync_channel::<Result<ChunkOutput>>(channel_capacity);
    let mut stats = new_stats(context.workflow);
    let mut failure: Option<Error> = None;

    std::thread::scope(|scope| {
        for _ in 0..threads {
            let sender = sender.clone();
            let abort = &abort;
            let window = &window;
            scope.spawn(move || {
                let mut block = ColumnarBlock::new(context.input.header());
                let mut dctx = zstd_safe::DCtx::create();
                let mut scratch = Scratch::default();
                while let Some(ordinal) = window.claim(chunks.len(), abort) {
                    let outcome = guarded(|| {
                        process_chunk(
                            context,
                            ordinal,
                            chunks[ordinal],
                            &mut block,
                            &mut dctx,
                            &mut scratch,
                        )
                    });
                    let failed = outcome.is_err();
                    if sender.send(outcome).is_err() || failed {
                        break;
                    }
                }
            });
        }
        // The committer holds the only remaining sender handle; dropping this
        // one lets `recv` terminate once every worker has finished.
        drop(sender);

        let mut pending: BTreeMap<usize, ChunkOutput> = BTreeMap::new();
        let mut expected = 0usize;
        while let Ok(message) = receiver.recv() {
            if failure.is_some() {
                continue; // keep draining so no worker blocks on send
            }
            match message {
                Err(error) => {
                    failure = Some(error);
                    abort.store(true, Ordering::Relaxed);
                    window.wake_all();
                }
                Ok(chunk) => {
                    pending.insert(chunk.ordinal, chunk);
                    while let Some(mut chunk) = pending.remove(&expected) {
                        if let Err(error) = commit_chunk(&mut chunk, outputs, &mut stats) {
                            failure = Some(error);
                            abort.store(true, Ordering::Relaxed);
                            window.wake_all();
                            break;
                        }
                        expected += 1;
                        window.commit_through(expected);
                    }
                }
            }
        }
    });

    match failure {
        Some(error) => Err(error),
        None => Ok(stats),
    }
}

/// Decodes one input block, processes its records and re-encodes the results.
fn process_chunk(
    context: &Context<'_>,
    ordinal: usize,
    block_index: usize,
    block: &mut ColumnarBlock,
    dctx: &mut zstd_safe::DCtx<'_>,
    scratch: &mut Scratch,
) -> Result<ChunkOutput> {
    let range = context.input.load(block_index, block, dctx)?;
    let orphan_fragment = |want: bool| -> Result<Option<Fragment>> {
        match (want, context.orphan_header) {
            (true, Some(header)) => Ok(Some(crate::io::fragment(header)?)),
            _ => Ok(None),
        }
    };
    let mut chunk = ChunkOutput {
        ordinal,
        accepted: crate::io::fragment(context.output_header)?,
        failed: if context.want_failed {
            Some(crate::io::fragment(context.output_header)?)
        } else {
            None
        },
        orphan_r1: orphan_fragment(context.want_orphan_r1)?,
        orphan_r2: orphan_fragment(context.want_orphan_r2)?,
        reasons: Vec::new(),
        corrections: Vec::new(),
        segments: Vec::new(),
        stats: new_stats(context.workflow),
    };
    let schema = context.schema;
    // Hoisted: the pipeline shape is fixed for the whole run, so the record loop
    // tests a local rather than reaching through the workflow every record.
    let segmenting = context.workflow.segment.is_some();

    for record in block.iter_records(range) {
        let index = record.index();
        if index < context.span.start || index >= context.span.end {
            continue;
        }

        let r1 = ReadView::new(
            record.sseq(),
            if schema.quality {
                Some(record.squal())
            } else {
                None
            },
            index,
            Mate::R1.name(),
        )?;
        let r2 = if schema.paired {
            Some(ReadView::new(
                record.xseq(),
                if schema.quality {
                    Some(record.xqual())
                } else {
                    None
                },
                index,
                Mate::R2.name(),
            )?)
        } else {
            None
        };

        if segmenting {
            segment_record(context, &mut chunk, &record, index, r1, scratch)?;
            continue;
        }

        let result = context
            .workflow
            .process_corrected(index, r1, r2, &mut scratch.correction)?;
        chunk.stats.record(&result);
        // Substitution values come from the plan, which holds the pre-mutation
        // bases by construction.
        for edit in &scratch.correction.edits {
            chunk
                .stats
                .record_substitution(edit.original_base, edit.corrected_base);
        }
        if context.want_corrections {
            write_correction_rows(
                context,
                &mut chunk.corrections,
                &record,
                &result,
                &scratch.correction,
            );
        }
        route_record(context, &mut chunk, &record, &result, &scratch.correction)?;
    }

    // Force every record into a completed block so the committer only has to
    // append already-compressed bytes in order.
    chunk.accepted.flush()?;
    if let Some(failed) = chunk.failed.as_mut() {
        failed.flush()?;
    }
    if let Some(orphan) = chunk.orphan_r1.as_mut() {
        orphan.flush()?;
    }
    if let Some(orphan) = chunk.orphan_r2.as_mut() {
        orphan.flush()?;
    }
    Ok(chunk)
}

/// Segments one source record and writes every fragment it produced.
///
/// One source record becomes zero, one or many output records. They are appended
/// in ascending coordinate order within the record, and records are visited in
/// index order, so the output is ordered by source index then segment index — the
/// same chunk-ordered commit that keeps the rest of the engine deterministic
/// needs no changes for this.
fn segment_record<R: BinseqRecord>(
    context: &Context<'_>,
    chunk: &mut ChunkOutput,
    record: &R,
    index: u64,
    read: ReadView<'_>,
    scratch: &mut Scratch,
) -> Result<()> {
    let outcome = context.workflow.process_segments(
        index,
        read,
        &mut scratch.segment,
        &mut scratch.fragments,
    )?;
    chunk.stats.record_segments(
        read.len(),
        outcome,
        scratch.segment.delimiters(),
        &scratch.fragments,
    );
    let schema = context.schema.unpaired();
    let names = context
        .workflow
        .segment
        .as_ref()
        .map_or(&[][..], |stage| stage.adapters.as_slice());
    for result in &scratch.fragments {
        if context.schema.headers {
            crate::segment::fragment_header(&mut scratch.header, record.sheader(), result.fragment);
        }
        if result.passed() {
            crate::io::push_segment(
                &mut chunk.accepted,
                schema,
                record,
                result.retained,
                &scratch.header,
            )?;
        } else if let Some(failed) = chunk.failed.as_mut() {
            // `original` cannot mean the whole source read here — that would put a
            // record of a different length in a per-fragment file — so the failed
            // output holds the fragment before trimming.
            let span = match context.failed_mode {
                FailedMode::Original => result.fragment.span,
                FailedMode::Processed => result.retained,
            };
            crate::io::push_segment(failed, schema, record, span, &scratch.header)?;
        }
        if context.want_segments {
            crate::io::write_segment_row(
                &mut chunk.segments,
                context.schema,
                record,
                result,
                |i| {
                    names
                        .get(i)
                        .map_or_else(|| ".".to_string(), |adapter| adapter.name.clone())
                },
            );
        }
        if context.want_reasons && !result.passed() {
            write_fragment_reason_row(&mut chunk.reasons, record.index(), result);
        }
    }
    Ok(())
}

/// Explains one rejected fragment in the failed-reason sidecar.
///
/// The `mate` column carries the segment index, because a fragment has no mate:
/// segmentation is single-end by definition.
fn write_fragment_reason_row(
    out: &mut Vec<u8>,
    index: u64,
    result: &crate::process::FragmentResult,
) {
    use std::io::Write as _;
    let _ = writeln!(
        out,
        "{index}\tsegment{}\tFAIL\t{}\t{}\t{}\t{}\t.\t.",
        result.fragment.index,
        result.reasons.label(),
        result.fragment.span.len(),
        result.fragment.span.len(),
        result.retained.len(),
    );
}

/// Appends this record's correction-log rows to the chunk-local buffer.
///
/// Chunk-local so workers never contend on a shared log; the committer appends
/// the buffers in chunk order, which makes the log ordered by record index, then
/// mate, then position — the order `plan` already sorted the edits into.
fn write_correction_rows<R: BinseqRecord>(
    context: &Context<'_>,
    out: &mut Vec<u8>,
    record: &R,
    result: &crate::process::PairResult,
    scratch: &crate::correct::CorrectionScratch,
) {
    use std::io::Write as _;

    let summary = result.correction;
    let Some(overlap) = result.overlap else {
        return;
    };
    if summary.corrected() == 0 {
        return; // only corrected pairs are logged
    }
    let disposition = result.disposition.name();
    let index = record.index();
    match context.correction_detail {
        crate::correct::LogDetail::Reads => {
            let _ = write!(out, "{index}\t");
            if context.schema.headers {
                crate::io::write_escaped(out, record.sheader());
                out.push(b'\t');
                crate::io::write_escaped(out, record.xheader());
            } else {
                out.extend_from_slice(b".\t.");
            }
            let _ = writeln!(
                out,
                "\t{}\t{}\t{}\t{}\t{}\t{}\t{disposition}",
                overlap.offset,
                overlap.overlap,
                summary.mismatches,
                summary.corrected_r1,
                summary.corrected_r2,
                summary.unresolved,
            );
        }
        crate::correct::LogDetail::Bases => {
            for edit in &scratch.edits {
                let _ = writeln!(
                    out,
                    "{index}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{disposition}",
                    edit.target.name(),
                    edit.target_position,
                    edit.donor().name(),
                    edit.donor_position,
                    char::from(edit.original_base),
                    char::from(edit.corrected_base),
                    crate::read::phred(edit.original_quality),
                    crate::read::phred(edit.corrected_quality),
                    overlap.offset,
                    overlap.overlap,
                );
            }
        }
    }
}

/// Routes one processed record to its destination fragments and the sidecar.
fn route_record<R: BinseqRecord>(
    context: &Context<'_>,
    chunk: &mut ChunkOutput,
    record: &R,
    result: &crate::process::PairResult,
    scratch: &crate::correct::CorrectionScratch,
) -> Result<()> {
    let schema = context.schema;
    let r2_span = result
        .r2
        .map_or(Span { start: 0, end: 0 }, |mate| mate.retained);
    // A corrected mate is written from worker scratch; everything else is
    // borrowed straight from the record.
    let corrected = |span: Span, mate: Mate| MateOutput::from_scratch(span, scratch, mate);
    match result.disposition {
        crate::process::PairDisposition::Accepted => {
            push_record(
                &mut chunk.accepted,
                schema,
                record,
                corrected(result.r1.retained, Mate::R1),
                corrected(r2_span, Mate::R2),
            )?;
            return Ok(());
        }
        crate::process::PairDisposition::Rejected => {
            if let Some(failed) = chunk.failed.as_mut() {
                // `original` means exactly that: the uncorrected pair as it was
                // read, so rejected data stays recoverable.
                let (r1_out, r2_out) = match context.failed_mode {
                    FailedMode::Original => (
                        MateOutput::borrowed(Span::full(result.r1.original_length)),
                        MateOutput::borrowed(Span::full(
                            result.r2.map_or(0, |mate| mate.original_length),
                        )),
                    ),
                    FailedMode::Processed => (
                        corrected(result.r1.retained, Mate::R1),
                        corrected(r2_span, Mate::R2),
                    ),
                };
                push_record(failed, schema, record, r1_out, r2_out)?;
            }
        }
        crate::process::PairDisposition::OrphanR1 => {
            // The surviving mate keeps its processed (trimmed) form.
            if let Some(orphan) = chunk.orphan_r1.as_mut() {
                crate::io::push_mate(
                    orphan,
                    schema.unpaired(),
                    record,
                    Mate::R1,
                    corrected(result.r1.retained, Mate::R1),
                )?;
            }
        }
        crate::process::PairDisposition::OrphanR2 => {
            if let Some(orphan) = chunk.orphan_r2.as_mut() {
                crate::io::push_mate(
                    orphan,
                    schema.unpaired(),
                    record,
                    Mate::R2,
                    corrected(r2_span, Mate::R2),
                )?;
            }
        }
    }
    // Every non-accepted record is explained in the sidecar, per mate.
    if context.want_reasons {
        for (mate, mate_result) in result.mates() {
            let adapter_name = mate_result.adapter_hit.and_then(|hit| {
                context
                    .workflow
                    .adapter
                    .as_ref()
                    .and_then(|stage| stage.adapters(mate).get(hit.adapter))
                    .map(|adapter| adapter.name.as_str())
            });
            write_reason_row(
                &mut chunk.reasons,
                record.index(),
                mate.name(),
                mate_result,
                adapter_name,
            );
        }
    }
    Ok(())
}

fn commit_chunk(
    chunk: &mut ChunkOutput,
    outputs: &mut Outputs<'_>,
    stats: &mut Stats,
) -> Result<()> {
    outputs.accepted.commit(&mut chunk.accepted)?;
    if let (Some(output), Some(fragment)) = (outputs.failed.as_mut(), chunk.failed.as_mut()) {
        output.commit(fragment)?;
    }
    if let (Some(output), Some(fragment)) = (outputs.orphan_r1.as_mut(), chunk.orphan_r1.as_mut()) {
        output.commit(fragment)?;
    }
    if let (Some(output), Some(fragment)) = (outputs.orphan_r2.as_mut(), chunk.orphan_r2.as_mut()) {
        output.commit(fragment)?;
    }
    if let Some(reasons) = outputs.reasons.as_mut() {
        reasons.write(&chunk.reasons)?;
    }
    if let Some(corrections) = outputs.corrections.as_mut() {
        corrections.write(&chunk.corrections)?;
    }
    if let Some(segments) = outputs.segments.as_mut() {
        segments.write(&chunk.segments)?;
    }
    stats.merge(&chunk.stats);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::WindowState;

    #[test]
    fn claim_window_advances_only_with_commit_progress() {
        let mut state = WindowState::default();
        let limit = 6;
        for ordinal in 0..limit {
            assert_eq!(state.claim(100, limit), Some(ordinal));
        }
        assert_eq!(state.claim(100, limit), None);

        state.committed = 2;
        assert_eq!(state.claim(100, limit), Some(6));
        assert_eq!(state.claim(100, limit), Some(7));
        assert_eq!(state.claim(100, limit), None);
        assert!(state.next - state.committed <= limit);
    }
}
