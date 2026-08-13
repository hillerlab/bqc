// Copyright (c) 2026 The Hiller Lab at the Senckenberg Gesellschaft für Naturforschung
// Distributed under the terms of the GNU General Public License, Version 3.0.

//! Exact deduplication with probabilistic acceleration.
//!
//! Two records are duplicates when their sequence payloads are byte-for-byte
//! equal; qualities, names and flags do not participate. For paired input the
//! ordered `(R1, R2)` pair is the key. The earliest occurrence survives and is
//! reproduced unchanged.
//!
//! The work is split into three order-preserving phases, each a parallel block
//! scan around a small serial critical section (the same pattern the processing
//! engine uses):
//!
//! ```text
//! pass 1   parallel decode/hash  ->  serial bitset update
//! pass 2a  parallel decode/hash  ->  serial exact classifier -> keep bitmaps
//! pass 2b  parallel decode/encode -> serial ordered commit
//! ```
//!
//! The exact classifier alone owns the candidate table, so first-occurrence
//! semantics and determinism need no locks. Bloom or hash collisions only cause
//! extra exact work; they can never drop a record, so the result is exact.

use std::collections::{BTreeMap, HashMap};
use std::path::Path;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::mpsc;

use binseq::BinseqRecord;
use binseq::cbq::ColumnarBlock;
use zstd::zstd_safe;

use crate::engine::CommitWindow;
use crate::error::{Error, Result};
use crate::io::{self, CbqInput, CbqOutput, Fragment, MateOutput, Schema};
use crate::read::Span;

/// Default memory budget in MiB for the Bloom filters and candidate arena.
pub const DEFAULT_MEMORY_MB: usize = 1024;

/// Deduplication options.
#[derive(Debug, Clone)]
pub struct DedupOptions {
    /// Memory budget in MiB. Split between the two Bloom filters; the exact
    /// candidate arena aborts rather than exceeding it.
    pub memory_mb: usize,
    /// Worker threads.
    pub threads: usize,
}

/// Counters for the dedup report.
#[derive(Debug, Clone, Default, serde::Serialize)]
pub struct DedupStats {
    pub records_seen: u64,
    pub records_kept: u64,
    pub duplicates_removed: u64,
    pub duplication_rate: f64,
    pub candidate_records: u64,
    pub candidate_families: u64,
    pub bloom_bytes: u64,
}

/// A one-bit-per-word membership filter. The fingerprint is never trusted for
/// equality, so a collision only causes extra exact work.
struct Bloom {
    bits: Vec<u64>,
    num_bits: usize,
}

impl Bloom {
    fn new(num_bits: usize) -> Self {
        Self {
            bits: vec![0u64; num_bits.div_ceil(64)],
            num_bits: num_bits.max(1),
        }
    }

    #[inline]
    fn contains(&self, h: u64) -> bool {
        let index = (h % self.num_bits as u64) as usize;
        self.bits[index / 64] & (1 << (index % 64)) != 0
    }

    #[inline]
    fn insert(&mut self, h: u64) {
        let index = (h % self.num_bits as u64) as usize;
        self.bits[index / 64] |= 1 << (index % 64);
    }

    fn bytes(&self) -> usize {
        self.bits.len() * 8
    }
}

/// A representative sequence stored in the candidate arena.
#[derive(Debug, Clone, Copy)]
struct CandidateRef {
    r1_start: usize,
    r1_len: usize,
    r2_start: usize,
    r2_len: usize,
}

/// Where a candidate record's bytes live in a worker's per-block sequence arena.
#[derive(Debug, Clone, Copy)]
struct CandidateSpan {
    r1_start: usize,
    r1_len: usize,
    r2_len: usize,
}

/// One record's fingerprint plus, for candidates, the location of its bytes.
#[derive(Debug, Clone, Copy)]
struct Fingerprinted {
    hash: u64,
    candidate: Option<CandidateSpan>,
}

/// What one decode worker sends to the exact classifier.
struct BlockHashes {
    records: Vec<Fingerprinted>,
    seqs: Vec<u8>,
}

/// The exact candidate table: a fingerprint-keyed map into a byte arena. Owned
/// by the serial classifier so first-occurrence order needs no lock.
struct Classifier {
    table: HashMap<u64, Vec<CandidateRef>>,
    arena: Vec<u8>,
    budget: usize,
}

impl Classifier {
    fn new(budget: usize) -> Self {
        Self {
            table: HashMap::new(),
            arena: Vec::new(),
            budget,
        }
    }

    /// Whether `(r1, r2)` already has a byte-identical representative.
    fn is_duplicate(&self, r1: &[u8], r2: Option<&[u8]>, h: u64) -> bool {
        let Some(bucket) = self.table.get(&h) else {
            return false;
        };
        bucket.iter().any(|candidate| {
            let r1_match =
                &self.arena[candidate.r1_start..candidate.r1_start + candidate.r1_len] == r1;
            let r2_match = match r2 {
                Some(r2) => {
                    &self.arena[candidate.r2_start..candidate.r2_start + candidate.r2_len] == r2
                }
                None => candidate.r2_len == 0,
            };
            r1_match && r2_match
        })
    }

    /// Appends a new representative and aborts when the arena exceeds the budget.
    fn store(&mut self, r1: &[u8], r2: Option<&[u8]>, h: u64) -> Result<()> {
        let r1_start = self.arena.len();
        self.arena.extend_from_slice(r1);
        let r2_start = self.arena.len();
        let r2_len = match r2 {
            Some(r2) => {
                self.arena.extend_from_slice(r2);
                r2.len()
            }
            None => 0,
        };
        self.table.entry(h).or_default().push(CandidateRef {
            r1_start,
            r1_len: r1.len(),
            r2_start,
            r2_len,
        });
        if self.arena.len() > self.budget {
            return Err(Error::config(
                "exact duplicate candidate table exceeded the configured memory budget; \
                 rerun with --memory-mb <larger>",
            ));
        }
        Ok(())
    }
}

/// Runs exact deduplication, writing the survivors to `output_path`.
pub fn run(
    input: &CbqInput,
    output_path: &Path,
    options: &DedupOptions,
    force: bool,
) -> Result<DedupStats> {
    let schema = input.schema();
    let header = input.header();
    let threads = options.threads.max(1);

    // Two ordinary bitsets. Seen and repeated split the budget evenly.
    let total_bits = options.memory_mb.saturating_mul(8 * 1024 * 1024).max(2);
    let mut seen = Bloom::new(total_bits / 2);
    let mut repeated = Bloom::new(total_bits - total_bits / 2);

    discover(input, schema, &mut seen, &mut repeated, threads)?;

    let mut stats = DedupStats {
        bloom_bytes: (seen.bytes() + repeated.bytes()) as u64,
        ..DedupStats::default()
    };
    // exact candidate families are RAM-backed. If real datasets
    // routinely hit this ceiling, spill hash partitions to disk; keep exact
    // equality and first-occurrence semantics unchanged.
    let mut classifier = Classifier::new(options.memory_mb.max(16).saturating_mul(1024 * 1024));

    let keep = classify(
        input,
        schema,
        &repeated,
        &mut classifier,
        threads,
        &mut stats,
    )?;
    encode(input, output_path, header, schema, &keep, force, threads)?;

    if stats.records_seen > 0 {
        stats.duplication_rate = stats.duplicates_removed as f64 / stats.records_seen as f64;
    }
    Ok(stats)
}

/// Pass 1: fingerprint every record in parallel, updating the two bitsets in
/// input order on the coordinator.
fn discover(
    input: &CbqInput,
    schema: Schema,
    seen: &mut Bloom,
    repeated: &mut Bloom,
    threads: usize,
) -> Result<()> {
    let blocks = input.blocks().len();
    if blocks == 0 {
        return Ok(());
    }
    let abort = AtomicBool::new(false);
    let window = CommitWindow::new(threads.saturating_mul(3).max(1));
    let (sender, receiver) = mpsc::sync_channel(threads.saturating_mul(2).max(1));
    let mut failure: Option<Error> = None;

    std::thread::scope(|scope| {
        for _ in 0..threads {
            let sender = sender.clone();
            let abort = &abort;
            let window = &window;
            scope.spawn(move || {
                let mut block = ColumnarBlock::new(input.header());
                let mut dctx = zstd_safe::DCtx::create();
                while let Some(ordinal) = window.claim(blocks, abort) {
                    let outcome = (|| {
                        let range = input.load(ordinal, &mut block, &mut dctx)?;
                        let mut hashes = Vec::new();
                        for record in block.iter_records(range) {
                            hashes.push(fingerprint(
                                record.sseq(),
                                schema.paired.then(|| record.xseq()),
                            ));
                        }
                        Ok(hashes)
                    })();
                    let failed = outcome.is_err();
                    if sender.send((ordinal, outcome)).is_err() || failed {
                        break;
                    }
                }
            });
        }
        drop(sender);

        let mut pending: BTreeMap<usize, Vec<u64>> = BTreeMap::new();
        let mut expected = 0usize;
        while let Ok((ordinal, outcome)) = receiver.recv() {
            if failure.is_some() {
                continue; // keep draining so no worker blocks on send
            }
            match outcome {
                Err(error) => {
                    failure = Some(error);
                    abort.store(true, Ordering::Relaxed);
                    window.wake_all();
                }
                Ok(hashes) => {
                    pending.insert(ordinal, hashes);
                    while let Some(hashes) = pending.remove(&expected) {
                        for h in hashes {
                            if seen.contains(h) {
                                repeated.insert(h);
                            }
                            seen.insert(h);
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
        None => Ok(()),
    }
}

/// Pass 2a: decode and fingerprint in parallel, classify in input order, and
/// return one keep bitmap per block.
fn classify(
    input: &CbqInput,
    schema: Schema,
    repeated: &Bloom,
    classifier: &mut Classifier,
    threads: usize,
    stats: &mut DedupStats,
) -> Result<Vec<Vec<bool>>> {
    let blocks = input.blocks().len();
    if blocks == 0 {
        return Ok(Vec::new());
    }
    let abort = AtomicBool::new(false);
    let window = CommitWindow::new(threads.saturating_mul(3).max(1));
    let (sender, receiver) = mpsc::sync_channel(threads.saturating_mul(2).max(1));
    let mut failure: Option<Error> = None;
    let mut keep_bitmaps: Vec<Vec<bool>> = Vec::with_capacity(blocks);

    std::thread::scope(|scope| {
        for _ in 0..threads {
            let sender = sender.clone();
            let abort = &abort;
            let window = &window;
            scope.spawn(move || {
                let mut block = ColumnarBlock::new(input.header());
                let mut dctx = zstd_safe::DCtx::create();
                while let Some(ordinal) = window.claim(blocks, abort) {
                    let outcome =
                        fingerprint_block(input, schema, repeated, ordinal, &mut block, &mut dctx);
                    let failed = outcome.is_err();
                    if sender.send((ordinal, outcome)).is_err() || failed {
                        break;
                    }
                }
            });
        }
        drop(sender);

        let mut pending: BTreeMap<usize, BlockHashes> = BTreeMap::new();
        let mut expected = 0usize;
        while let Ok((ordinal, outcome)) = receiver.recv() {
            if failure.is_some() {
                continue;
            }
            match outcome {
                Err(error) => {
                    failure = Some(error);
                    abort.store(true, Ordering::Relaxed);
                    window.wake_all();
                }
                Ok(block_hashes) => {
                    pending.insert(ordinal, block_hashes);
                    while let Some(block_hashes) = pending.remove(&expected) {
                        match classify_block(&block_hashes, schema, classifier, stats) {
                            Ok(keep) => keep_bitmaps.push(keep),
                            Err(error) => {
                                failure = Some(error);
                                abort.store(true, Ordering::Relaxed);
                                window.wake_all();
                                break;
                            }
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
        None => Ok(keep_bitmaps),
    }
}

/// Decodes one block, fingerprinting every record and copying the bytes of
/// candidate records into a per-block arena.
fn fingerprint_block(
    input: &CbqInput,
    schema: Schema,
    repeated: &Bloom,
    ordinal: usize,
    block: &mut ColumnarBlock,
    dctx: &mut zstd_safe::DCtx<'_>,
) -> Result<BlockHashes> {
    let range = input.load(ordinal, block, dctx)?;
    let mut records = Vec::new();
    let mut seqs = Vec::new();
    for record in block.iter_records(range) {
        let r1 = record.sseq();
        let r2 = schema.paired.then(|| record.xseq());
        let hash = fingerprint(r1, r2);
        let candidate = if repeated.contains(hash) {
            let r1_start = seqs.len();
            let r1_len = r1.len();
            seqs.extend_from_slice(r1);
            let r2_len = match r2 {
                Some(r2) => {
                    seqs.extend_from_slice(r2);
                    r2.len()
                }
                None => 0,
            };
            Some(CandidateSpan {
                r1_start,
                r1_len,
                r2_len,
            })
        } else {
            None
        };
        records.push(Fingerprinted { hash, candidate });
    }
    Ok(BlockHashes { records, seqs })
}

/// Classifies one block's fingerprints in record order, mutating the candidate
/// table and returning the keep bitmap.
fn classify_block(
    block_hashes: &BlockHashes,
    schema: Schema,
    classifier: &mut Classifier,
    stats: &mut DedupStats,
) -> Result<Vec<bool>> {
    let mut keep = Vec::with_capacity(block_hashes.records.len());
    for fingerprinted in &block_hashes.records {
        stats.records_seen += 1;
        let kept = match fingerprinted.candidate {
            None => true,
            Some(span) => {
                stats.candidate_records += 1;
                let r1 = &block_hashes.seqs[span.r1_start..span.r1_start + span.r1_len];
                let r2 = schema.paired.then(|| {
                    &block_hashes.seqs
                        [span.r1_start + span.r1_len..span.r1_start + span.r1_len + span.r2_len]
                });
                let duplicate = classifier.is_duplicate(r1, r2, fingerprinted.hash);
                if duplicate {
                    stats.duplicates_removed += 1;
                } else {
                    classifier.store(r1, r2, fingerprinted.hash)?;
                    stats.candidate_families += 1;
                }
                !duplicate
            }
        };
        if kept {
            stats.records_kept += 1;
        }
        keep.push(kept);
    }
    Ok(keep)
}

/// Pass 2b: re-decode and encode each block in parallel, committing the
/// already-compressed blocks in input order.
fn encode(
    input: &CbqInput,
    output_path: &Path,
    header: binseq::cbq::FileHeader,
    schema: Schema,
    keep: &[Vec<bool>],
    force: bool,
    threads: usize,
) -> Result<()> {
    let mut output = CbqOutput::create(output_path, input.path(), header, force)?;
    let blocks = input.blocks().len();
    if blocks == 0 {
        return output.finish();
    }
    let abort = AtomicBool::new(false);
    let window = CommitWindow::new(threads.saturating_mul(3).max(1));
    let (sender, receiver) = mpsc::sync_channel(threads.saturating_mul(2).max(1));
    let mut failure: Option<Error> = None;

    std::thread::scope(|scope| {
        for _ in 0..threads {
            let sender = sender.clone();
            let abort = &abort;
            let window = &window;
            scope.spawn(move || {
                let mut block = ColumnarBlock::new(header);
                let mut dctx = zstd_safe::DCtx::create();
                while let Some(ordinal) = window.claim(blocks, abort) {
                    let outcome = encode_block(
                        input,
                        schema,
                        header,
                        ordinal,
                        &keep[ordinal],
                        &mut block,
                        &mut dctx,
                    );
                    let failed = outcome.is_err();
                    if sender.send((ordinal, outcome)).is_err() || failed {
                        break;
                    }
                }
            });
        }
        drop(sender);

        let mut pending: BTreeMap<usize, Fragment> = BTreeMap::new();
        let mut expected = 0usize;
        while let Ok((ordinal, outcome)) = receiver.recv() {
            if failure.is_some() {
                continue;
            }
            match outcome {
                Err(error) => {
                    failure = Some(error);
                    abort.store(true, Ordering::Relaxed);
                    window.wake_all();
                }
                Ok(fragment) => {
                    pending.insert(ordinal, fragment);
                    while let Some(mut fragment) = pending.remove(&expected) {
                        if let Err(error) = output.commit(&mut fragment) {
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
        None => output.finish(),
    }
}

/// Re-decodes one block and writes its surviving records.
fn encode_block(
    input: &CbqInput,
    schema: Schema,
    header: binseq::cbq::FileHeader,
    ordinal: usize,
    keep: &[bool],
    block: &mut ColumnarBlock,
    dctx: &mut zstd_safe::DCtx<'_>,
) -> Result<Fragment> {
    let range = input.load(ordinal, block, dctx)?;
    let mut fragment = io::fragment(header)?;
    for (record, &kept) in block.iter_records(range).zip(keep) {
        if kept {
            push_record(&mut fragment, schema, &record)?;
        }
    }
    fragment.flush()?;
    Ok(fragment)
}

/// FNV-1a 64-bit. Deterministic and dependency-free; the fingerprint is never
/// trusted for equality, so a collision is harmless.
fn fnv1a(bytes: &[u8]) -> u64 {
    let mut hash = 0xcbf2_9ce4_8422_2325u64;
    for &byte in bytes {
        hash ^= u64::from(byte);
        hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
    }
    hash
}

/// Fingerprint of one record. The R1 length is mixed in so `(AC, CGT)` cannot
/// alias `(ACC, GT)`.
fn fingerprint(r1: &[u8], r2: Option<&[u8]>) -> u64 {
    let mut hash = fnv1a(r1);
    if let Some(r2) = r2 {
        hash ^= r1.len() as u64;
        hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
        hash ^= fnv1a(r2);
    }
    hash
}

/// Writes one record unchanged, propagating sequence, quality, header and flag.
fn push_record<R: BinseqRecord>(fragment: &mut Fragment, schema: Schema, record: &R) -> Result<()> {
    io::push_record(
        fragment,
        schema,
        record,
        MateOutput::borrowed(Span::full(record.sseq().len())),
        MateOutput::borrowed(Span::full(if schema.paired {
            record.xseq().len()
        } else {
            0
        })),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fingerprints_do_not_alias_across_mate_boundaries() {
        // (R1=AC, R2=CGT) must not fingerprint like (R1=ACC, R2=GT).
        let a = fingerprint(b"AC", Some(b"CGT"));
        let b = fingerprint(b"ACC", Some(b"GT"));
        assert_ne!(a, b);
    }

    #[test]
    fn exact_equality_distinguishes_colliding_fingerprints() {
        let mut classifier = Classifier::new(1024 * 1024);
        // Force a collision by using the same key for two distinct sequences.
        classifier.store(b"ACGT", None, 7).unwrap();
        assert!(classifier.is_duplicate(b"ACGT", None, 7));
        assert!(!classifier.is_duplicate(b"ACGA", None, 7));
    }
}
