// Copyright (c) 2026 The Hiller Lab at the Senckenberg Gesellschaft für Naturforschung
// Distributed under the terms of the GNU General Public License, Version 3.0.

//! Composite adapter discovery.
//!
//! One detector, not a menu of algorithms. A user asking "what adapter is in
//! this file" should get the best answer the evidence supports without having to
//! know which internal method produced it.
//!
//! ```text
//! pass 1   known-library scan (seed-indexed)      -> support per library entry
//!          exact 3' 10-mer counting               -> enriched seeds
//!          paired-overlap overhangs               -> adapter past the insert
//! pass 2   consensus extension around each seed   -> assembled candidates
//!          deduplicate against the known entries
//! pass 3   verify every candidate with the trimming matcher
//!            -> support, positional evidence, error rates
//!          gates    -> high | medium | low
//!          decision -> confident | mixed | inconclusive
//! ```
//!
//! Every pass reads the *same* deterministic sample, so each measures exactly
//! the reads the previous one drew its conclusions from.
//!
//! **Support is measured once.** A read can carry evidence for a known entry, an
//! enriched k-mer and an overlap overhang at the same time; adding those tallies
//! would count it three times. The verification pass is therefore the single
//! source of every support number, and the evidence sources a candidate came
//! from are recorded as provenance only.
//!
//! Determinism comes from the merge being integer addition over per-worker
//! tallies: which worker sees which block does not change a sum, so the result
//! is identical at any thread count. Nothing here retains reads or allocates per
//! record.

use serde::Serialize;

use crate::adapter::{Adapter, AdapterHit, AdapterParams, verify_at};
use crate::detect::{KNOWN_ADAPTERS, KNOWN_ADAPTERS_VERSION, SeedIndex};
use crate::error::Result;
use crate::io::CbqInput;
use crate::process::Mate;
use crate::sniff::sample::{BlockWork, SamplePlan, block_work, for_each_sampled};
use crate::sniff::{Confidence, Decision, EvidenceSource, Gates, fraction, related};

/// Window at the 3' end treated as the tail for positional evidence.
///
/// Also the k-mer collection window in the de novo stage, so tail statistics
/// mean the same thing to every evidence source.
pub const TAIL_WINDOW: usize = 60;

/// Library entries carried from the scan into verification, per mate.
///
/// Verification costs one matcher pass per candidate per sampled read, so the
/// shortlist is small and fixed. Entries below it cannot reach `medium` anyway.
const SHORTLIST: usize = 8;

/// Length of the de novo enrichment k-mer.
///
/// Ten bases is short enough that a real adapter's seed survives the sequencing
/// errors in a sample, and long enough that a specific 10-mer is rare by chance:
/// `4^10` is a million, against tens of thousands of sampled tails.
pub const KMER_K: usize = 10;

/// Number of distinct k-mers, and so the length of the dense counter.
///
/// An exact `Vec<u32>` over this key space is 4 MiB — smaller and faster than a
/// hash table at this size, with no counting collisions to reason about. A
/// probabilistic sketch would buy nothing here.
const KMER_SPACE: usize = 1 << (2 * KMER_K);

/// Mask keeping the low `KMER_K` bases of a rolled 2-bit k-mer.
const KMER_MASK: u32 = (KMER_SPACE - 1) as u32;

/// Enriched seeds carried into consensus extension, per mate.
const SEED_SHORTLIST: usize = 32;

/// Longest sequence a consensus may reach, including its seed.
const MAX_CONSENSUS: usize = 60;

/// Most k-mers a tail window can hold, for the per-read stack buffer.
const MAX_TAIL_KMERS: usize = TAIL_WINDOW - KMER_K + 1;

/// Largest match coordinate tracked exactly in the position histograms.
const MAX_POSITION: usize = 1024;

/// Share of a tail window that one base must occupy to be called an artifact.
const ARTIFACT_SHARE: f64 = 0.9;

/// Bases permitted after the median candidate end while still calling it
/// terminal. This absorbs one- or two-base alignment shadows without treating a
/// genuine internal motif as read-through.
const END_REACH_SLACK: u64 = 3;

/// Adapter-sniffing thresholds.
#[derive(Debug, Clone, Copy, Serialize)]
pub struct Params {
    /// Records sampled across the selected range.
    pub sample_size: u64,
    /// Candidates reported per mate.
    pub top: usize,
    /// Matcher thresholds used for verification.
    pub matcher: AdapterParams,
    /// Classification gates.
    pub gates: Gates,
}

impl Default for Params {
    fn default() -> Self {
        Self {
            sample_size: 262_144,
            top: 5,
            matcher: AdapterParams::default(),
            gates: Gates::default(),
        }
    }
}

impl Params {
    /// Validates sniffing thresholds.
    pub fn validate(self) -> Result<Self> {
        if self.sample_size == 0 {
            return Err(crate::error::Error::config(
                "--sample-size must be at least 1",
            ));
        }
        if self.top == 0 {
            return Err(crate::error::Error::config("--top must be at least 1"));
        }
        Ok(self)
    }
}

/// One reported adapter candidate.
#[derive(Debug, Clone, Serialize)]
pub struct Candidate {
    #[serde(serialize_with = "crate::report::serialize_bytes")]
    pub sequence: Vec<u8>,
    pub known_name: Option<String>,
    pub known_category: Option<String>,
    pub confidence: Confidence,
    pub evidence_sources: Vec<EvidenceSource>,
    pub supporting_reads: u64,
    pub support_fraction: f64,
    pub tail_connected_fraction: f64,
    /// Normalized tail/body enrichment. `None` means tail-only support: the body
    /// rate is zero while the tail rate is nonzero.
    pub tail_enrichment: Option<f64>,
    pub median_start: u64,
    pub median_distance_to_end: u64,
    pub exact_matches: u64,
    pub substitution_matches: u64,
    pub indel_matches: u64,
    pub mean_error_rate: f64,
    pub paired_overlap_support: u64,
    /// Declaration order in the bundled library; orders otherwise-tied
    /// candidates. Not part of the user-facing result.
    #[serde(skip)]
    library_index: usize,
}

/// The result for one mate.
#[derive(Debug, Clone, Serialize)]
pub struct MateResult {
    pub mate: &'static str,
    pub decision: Decision,
    pub recommended_sequence: Option<String>,
    pub recommended_name: Option<String>,
    pub candidates: Vec<Candidate>,
    pub sampled_reads: u64,
    pub informative_reads: u64,
    /// Sampled reads whose 3' window is a poly-A run. Reported as an artifact
    /// signal, never recommended as an adapter.
    pub poly_a_signal: u64,
    pub poly_g_signal: u64,
}

/// The whole `sniff adapters` result.
#[derive(Debug, Clone, Serialize)]
pub struct AdapterSniff {
    pub database_version: u32,
    pub r1: MateResult,
    pub r2: Option<MateResult>,
}

impl AdapterSniff {
    /// The strictest decision across the mates present.
    ///
    /// A file is only confidently characterised when every mate is, so a
    /// pipeline gating on the whole file cannot be satisfied by one good mate.
    #[must_use]
    pub fn decision(&self) -> Decision {
        let mates = std::iter::once(&self.r1).chain(self.r2.as_ref());
        mates
            .map(|mate| mate.decision)
            .max()
            .unwrap_or(Decision::Confident)
    }

    /// Mates whose recommendation can be written into a configuration file.
    #[must_use]
    pub fn recommendation(&self) -> Option<(String, Option<String>)> {
        if !self.decision().is_confident() {
            return None;
        }
        let r1 = self.r1.recommended_sequence.clone()?;
        let r2 = match self.r2.as_ref() {
            Some(mate) => Some(mate.recommended_sequence.clone()?),
            None => None,
        };
        Some((r1, r2))
    }
}

/// Pass-1 tallies for one mate. Merged by addition.
struct Scan {
    sampled: u64,
    informative: u64,
    poly_a: u64,
    poly_g: u64,
    /// Searchable start coordinates in the tail window, summed over reads.
    tail_positions: u64,
    /// Searchable start coordinates outside it.
    body_positions: u64,
    hits: Vec<u64>,
    /// Reads containing each tail k-mer. Empty when de novo discovery is off.
    kmers: Vec<u32>,
    /// Bases past the inferred insert boundary, voted column by column.
    overhang: Consensus,
}

impl Scan {
    fn new(entries: usize, de_novo: bool) -> Self {
        Self {
            sampled: 0,
            informative: 0,
            poly_a: 0,
            poly_g: 0,
            tail_positions: 0,
            body_positions: 0,
            hits: vec![0; entries],
            kmers: if de_novo {
                vec![0; KMER_SPACE]
            } else {
                Vec::new()
            },
            overhang: Consensus::new(),
        }
    }

    fn merge(&mut self, other: &Self) {
        self.sampled += other.sampled;
        self.informative += other.informative;
        self.poly_a += other.poly_a;
        self.poly_g += other.poly_g;
        self.tail_positions += other.tail_positions;
        self.body_positions += other.body_positions;
        for (slot, add) in self.hits.iter_mut().zip(&other.hits) {
            *slot += add;
        }
        for (slot, add) in self.kmers.iter_mut().zip(&other.kmers) {
            *slot += add;
        }
        self.overhang.merge(&other.overhang);
    }

    /// Counts every distinct tail k-mer of one read once.
    ///
    /// Read support, not occurrence count: a tandem repeat would otherwise
    /// out-vote a genuine adapter that appears once per read. Deduplication uses
    /// a fixed stack buffer — at most `TAIL_WINDOW - KMER_K + 1` codes — rather
    /// than a hash set per read.
    fn count_tail_kmers(&mut self, sequence: &[u8]) {
        if self.kmers.is_empty() {
            return;
        }
        let window = &sequence[sequence.len().saturating_sub(TAIL_WINDOW)..];
        let mut codes = [0u32; MAX_TAIL_KMERS];
        let mut found = 0usize;
        let mut code = 0u32;
        let mut valid = 0usize;
        for &base in window {
            let Some(slot) = crate::detect::base_slot(base) else {
                // An `N` invalidates every window it belongs to.
                valid = 0;
                continue;
            };
            code = ((code << 2) | slot as u32) & KMER_MASK;
            valid += 1;
            if valid >= KMER_K && found < codes.len() {
                codes[found] = code;
                found += 1;
            }
        }
        let codes = &mut codes[..found];
        codes.sort_unstable();
        let mut previous = None;
        for &code in codes.iter() {
            if previous != Some(code) {
                self.kmers[code as usize] += 1;
                previous = Some(code);
            }
        }
    }

    /// The library entries worth verifying, best supported first.
    fn shortlist(&self) -> Vec<usize> {
        let mut order: Vec<usize> = (0..self.hits.len()).filter(|&i| self.hits[i] > 0).collect();
        // Most hits first; ties to declaration order, as everywhere else.
        order.sort_by_key(|&i| (std::cmp::Reverse(self.hits[i]), i));
        order.truncate(SHORTLIST);
        order
    }

    /// The enriched, informative k-mer seeds worth extending into a consensus.
    ///
    /// Ties break to the smaller encoding so the shortlist never depends on scan
    /// order. Homopolymers and low-complexity runs are excluded here rather than
    /// after extension: they are the most abundant tail k-mers in real data and
    /// would otherwise fill the whole shortlist.
    fn seeds(&self, gates: Gates) -> Vec<(u32, u64)> {
        if self.kmers.is_empty() {
            return Vec::new();
        }
        let floor = gates.min_supporting_reads.max(1);
        let mut best: Vec<(u32, u64)> = Vec::with_capacity(SEED_SHORTLIST + 1);
        for (code, &count) in self.kmers.iter().enumerate() {
            let count = u64::from(count);
            if count < floor {
                continue;
            }
            let code = code as u32;
            if !informative_seed(code, gates) {
                continue;
            }
            let position = best.partition_point(|&(other_code, other)| {
                other > count || (other == count && other_code < code)
            });
            if position < SEED_SHORTLIST {
                best.insert(position, (code, count));
                best.truncate(SEED_SHORTLIST);
            }
        }
        best
    }
}

/// Longest repeat period treated as a low-complexity artifact.
const MAX_ARTIFACT_PERIOD: usize = 3;

/// Whether a k-mer can be adapter evidence at all.
///
/// Poly-A and poly-G tails are the most abundant 3' k-mers in real data and are
/// instrument or biology artifacts, never adapters. They are counted and
/// reported separately rather than silently dropped.
///
/// Three tests, because no one of them is sufficient. A dominant base catches
/// homopolymers. Adjacent-base complexity — fastp's metric, and the one the
/// filter stage already uses — catches runs. Neither catches a short-period
/// repeat: `ATATATATAT` has two bases in equal share and changes at *every*
/// position, so it scores a perfect 1.0. Such a tract would then extend into a
/// long, perfectly stable consensus and be recommended as an adapter, which is
/// precisely the failure the seed exclusion exists to prevent.
fn informative_seed(code: u32, gates: Gates) -> bool {
    let bases = crate::detect::decode_kmer(code, KMER_K);
    let mut counts = [0u64; 4];
    for &base in &bases {
        if let Some(slot) = crate::detect::base_slot(base) {
            counts[slot] += 1;
        }
    }
    let dominant = counts.iter().max().copied().unwrap_or(0);
    if dominant as f64 > crate::detect::MAX_SEED_BASE_SHARE * KMER_K as f64 {
        return false;
    }
    if periodic(&bases, MAX_ARTIFACT_PERIOD) {
        return false;
    }
    crate::filter::complexity(&bases) >= gates.min_complexity
}

/// Whether `bases` is an exact repeat of a unit no longer than `max_period`.
fn periodic(bases: &[u8], max_period: usize) -> bool {
    (1..=max_period.min(bases.len().saturating_sub(1)))
        .any(|period| bases.iter().skip(period).zip(bases).all(|(a, b)| a == b))
}

/// Pass-2 tallies for one candidate. Merged by addition.
#[derive(Clone)]
struct Verified {
    supporting: u64,
    tail_connected: u64,
    tail_hits: u64,
    body_hits: u64,
    exact: u64,
    substitution: u64,
    indel: u64,
    error_sum: u64,
    overlap_sum: u64,
    starts: Vec<u32>,
    distances: Vec<u32>,
}

impl Verified {
    fn new() -> Self {
        Self {
            supporting: 0,
            tail_connected: 0,
            tail_hits: 0,
            body_hits: 0,
            exact: 0,
            substitution: 0,
            indel: 0,
            error_sum: 0,
            overlap_sum: 0,
            starts: vec![0; MAX_POSITION + 1],
            distances: vec![0; MAX_POSITION + 1],
        }
    }

    fn merge(&mut self, other: &Self) {
        self.supporting += other.supporting;
        self.tail_connected += other.tail_connected;
        self.tail_hits += other.tail_hits;
        self.body_hits += other.body_hits;
        self.exact += other.exact;
        self.substitution += other.substitution;
        self.indel += other.indel;
        self.error_sum += other.error_sum;
        self.overlap_sum += other.overlap_sum;
        for (slot, add) in self.starts.iter_mut().zip(&other.starts) {
            *slot += add;
        }
        for (slot, add) in self.distances.iter_mut().zip(&other.distances) {
            *slot += add;
        }
    }

    fn record(&mut self, hit: &AdapterHit, read_len: usize, adapter_len: usize) {
        self.supporting += 1;
        if hit.start + adapter_len >= read_len {
            self.tail_connected += 1;
        }
        if hit.start + TAIL_WINDOW >= read_len {
            self.tail_hits += 1;
        } else {
            self.body_hits += 1;
        }
        match (hit.errors, hit.consumed == hit.overlap) {
            (0, _) => self.exact += 1,
            (_, true) => self.substitution += 1,
            (_, false) => self.indel += 1,
        }
        self.error_sum += hit.errors as u64;
        self.overlap_sum += hit.overlap as u64;
        self.starts[hit.start.min(MAX_POSITION)] += 1;
        self.distances[read_len.saturating_sub(hit.start).min(MAX_POSITION)] += 1;
    }
}

/// The value at or below which half the observations fall.
fn median(histogram: &[u32], total: u64) -> u64 {
    if total == 0 {
        return 0;
    }
    let half = total.div_ceil(2);
    let mut seen = 0u64;
    for (value, &count) in histogram.iter().enumerate() {
        seen += u64::from(count);
        if seen >= half {
            return value as u64;
        }
    }
    0
}

/// Searchable start coordinates in a read, split at the tail window.
fn positions(read_len: usize, min_overlap: usize) -> (u64, u64) {
    if read_len < min_overlap {
        return (0, 0);
    }
    let total = (read_len - min_overlap + 1) as u64;
    let tail = TAIL_WINDOW.min(read_len) as u64;
    let tail = tail.min(total);
    (tail, total - tail)
}

/// Whether a read's 3' window is a homopolymer artifact.
fn artifact(sequence: &[u8], base: u8) -> bool {
    let window = &sequence[sequence.len().saturating_sub(TAIL_WINDOW)..];
    if window.len() < TAIL_WINDOW {
        return false;
    }
    let matching = window
        .iter()
        .fold(0usize, |n, &b| n + usize::from(b == base));
    matching as f64 >= ARTIFACT_SHARE * window.len() as f64
}

/// Column votes around one enriched seed, accumulated streaming.
///
/// The reads themselves are never retained: each sampled read votes on the
/// columns around the seed as it goes past, and the consensus is read off the
/// merged columns afterwards. Retaining the sample instead would cost one heap
/// allocation per read and hold the whole sample in memory at once.
#[derive(Clone)]
struct Consensus {
    /// Reads that contained the seed and voted.
    voting: u64,
    /// Columns before the seed, nearest first.
    left: Vec<[u64; 4]>,
    /// Columns from the seed's first base onwards.
    right: Vec<[u64; 4]>,
}

impl Consensus {
    fn new() -> Self {
        Self {
            voting: 0,
            left: Vec::new(),
            right: Vec::new(),
        }
    }

    fn merge(&mut self, other: &Self) {
        self.voting += other.voting;
        for (columns, add) in [
            (&mut self.left, &other.left),
            (&mut self.right, &other.right),
        ] {
            if columns.len() < add.len() {
                columns.resize(add.len(), [0; 4]);
            }
            for (column, extra) in columns.iter_mut().zip(add) {
                for (slot, count) in column.iter_mut().zip(extra) {
                    *slot += count;
                }
            }
        }
    }

    /// Votes a sequence known to begin at the adapter's first base.
    ///
    /// An overlap overhang needs no leftward walk: the insert boundary already
    /// names where the adapter starts, so there is nothing before it to assemble.
    fn observe_prefix(&mut self, sequence: &[u8]) {
        self.voting += 1;
        for (column, &base) in sequence.iter().take(MAX_CONSENSUS).enumerate() {
            crate::detect::vote(&mut self.right, column, base);
        }
    }

    /// Votes the bases around one occurrence of the seed.
    fn observe(&mut self, sequence: &[u8], start: usize) {
        self.voting += 1;
        for (column, &base) in sequence[start..].iter().take(MAX_CONSENSUS).enumerate() {
            crate::detect::vote(&mut self.right, column, base);
        }
        for distance in 1..=start.min(MAX_CONSENSUS) {
            crate::detect::vote(&mut self.left, distance - 1, sequence[start - distance]);
        }
    }

    /// Reads the consensus off the columns.
    ///
    /// Extension stops at the first column without enough coverage or without a
    /// stable majority — the same rule in both directions, so a candidate never
    /// grows past the evidence supporting it. The leftward walk matters: an
    /// enriched seed is usually an *internal* k-mer of the adapter, not its 5'
    /// start, so the sequence a user needs begins to the left of the seed.
    fn sequence(&self) -> Vec<u8> {
        let mut left: Vec<u8> = Vec::new();
        for column in &self.left {
            match crate::detect::majority_base(column, self.voting) {
                Some(base) => left.push(base),
                None => break,
            }
        }
        left.reverse();
        let mut consensus = left;
        for column in &self.right {
            if consensus.len() >= MAX_CONSENSUS {
                break;
            }
            match crate::detect::majority_base(column, self.voting) {
                Some(base) => consensus.push(base),
                None => break,
            }
        }
        if consensus.len() > MAX_CONSENSUS {
            consensus.truncate(MAX_CONSENSUS);
        }
        consensus
    }
}

/// The mates present in one sampled record.
fn mates<'a>(r1: &'a [u8], r2: Option<&'a [u8]>) -> impl Iterator<Item = (Mate, &'a [u8])> {
    std::iter::once((Mate::R1, r1)).chain(r2.map(|sequence| (Mate::R2, sequence)))
}

/// Per-mate worker state for the scan pass.
struct ScanState {
    r1: Scan,
    r2: Scan,
}

impl ScanState {
    fn mate(&mut self, mate: Mate) -> &mut Scan {
        match mate {
            Mate::R1 => &mut self.r1,
            Mate::R2 => &mut self.r2,
        }
    }
}

/// Per-mate worker state for the consensus pass.
struct ConsensusState {
    r1: Vec<Consensus>,
    r2: Vec<Consensus>,
}

/// Votes the columns around every shortlisted seed, over the same sample.
///
/// One rolling pass per read serves every seed: the tail window's k-mers are
/// rolled once and each is looked up in the sorted shortlist, so the cost does
/// not grow with the number of seeds.
fn consensus_pass(
    input: &CbqInput,
    work: &[BlockWork],
    plan: &SamplePlan,
    threads: usize,
    seeds: &[Vec<(u32, u64)>; 2],
) -> Result<[Vec<Consensus>; 2]> {
    // `(code, shortlist position)`, sorted by code for bisection.
    let lookup: [Vec<(u32, usize)>; 2] = std::array::from_fn(|mate| {
        let mut table: Vec<(u32, usize)> = seeds[mate]
            .iter()
            .enumerate()
            .map(|(position, &(code, _))| (code, position))
            .collect();
        table.sort_unstable();
        table
    });

    let states = for_each_sampled(
        input,
        work,
        plan,
        threads,
        || ConsensusState {
            r1: vec![Consensus::new(); seeds[0].len()],
            r2: vec![Consensus::new(); seeds[1].len()],
        },
        |state, r1, r2| {
            for (mate, sequence) in mates(r1, r2) {
                observe_seeds(state, &lookup[mate.index()], mate, sequence);
            }
        },
    )?;

    let mut r1 = vec![Consensus::new(); seeds[0].len()];
    let mut r2 = vec![Consensus::new(); seeds[1].len()];
    for state in &states {
        for (slot, add) in r1.iter_mut().zip(&state.r1) {
            slot.merge(add);
        }
        for (slot, add) in r2.iter_mut().zip(&state.r2) {
            slot.merge(add);
        }
    }
    Ok([r1, r2])
}

/// Votes every shortlisted seed that occurs in one read.
///
/// One rolling pass serves every seed: the tail window's k-mers are rolled once
/// and each looked up in the sorted shortlist, so the cost does not grow with
/// the number of seeds.
fn observe_seeds(state: &mut ConsensusState, table: &[(u32, usize)], mate: Mate, sequence: &[u8]) {
    {
        if table.is_empty() {
            return;
        }
        let columns = match mate {
            Mate::R1 => &mut state.r1,
            Mate::R2 => &mut state.r2,
        };
        let offset = sequence.len().saturating_sub(TAIL_WINDOW);
        let window = &sequence[offset..];
        let mut code = 0u32;
        let mut valid = 0usize;
        let mut seen = [false; SEED_SHORTLIST];
        for (position, &base) in window.iter().enumerate() {
            let Some(slot) = crate::detect::base_slot(base) else {
                valid = 0;
                continue;
            };
            code = ((code << 2) | slot as u32) & KMER_MASK;
            valid += 1;
            if valid < KMER_K {
                continue;
            }
            let index = table.partition_point(|&(other, _)| other < code);
            let Some(&(found, target)) = table.get(index) else {
                continue;
            };
            // Only the first occurrence in a read votes, so a repeat cannot
            // vote twice and skew its own consensus.
            if found != code || seen[target] {
                continue;
            }
            seen[target] = true;
            columns[target].observe(sequence, offset + position + 1 - KMER_K);
        }
    }
}

/// Per-mate worker state for the verification pass.
struct VerifyState {
    r1: Vec<Verified>,
    r2: Vec<Verified>,
}

/// Runs adapter discovery over a deterministic distributed sample.
pub fn sniff(
    input: &CbqInput,
    params: Params,
    span: Option<&std::ops::Range<u64>>,
    threads: usize,
) -> Result<(AdapterSniff, SamplePlan)> {
    let paired = input.schema().paired;
    let plan = SamplePlan::for_input(input, span, params.sample_size);
    let work = block_work(&plan, input);

    let library = [
        SeedIndex::new(Mate::R1, params.matcher),
        SeedIndex::new(Mate::R2, params.matcher),
    ];
    let scans = scan_pass(input, &work, &plan, threads, params, &library)?;

    // De novo extension: enriched tail seeds become consensus sequences.
    let seeds: [Vec<(u32, u64)>; 2] = std::array::from_fn(|mate| {
        if mate == 0 || paired {
            scans[mate].seeds(params.gates)
        } else {
            Vec::new()
        }
    });
    let consensus = consensus_pass(input, &work, &plan, threads, &seeds)?;

    // Known and de novo evidence meet here, and are deduplicated before the
    // expensive verification rather than after it.
    let proposals: [Vec<Proposed>; 2] =
        std::array::from_fn(|mate| propose(&scans[mate], &library[mate], &consensus[mate]));
    let adapters: [Vec<Adapter>; 2] = std::array::from_fn(|mate| {
        proposals[mate]
            .iter()
            .map(|proposal| proposal.adapter.clone())
            .collect()
    });
    let tallies = verify_pass(input, &work, &plan, threads, params, &adapters)?;

    let r1 = assemble(Mate::R1, &scans[0], &proposals[0], &tallies[0], params);
    let r2 = paired.then(|| assemble(Mate::R2, &scans[1], &proposals[1], &tallies[1], params));

    Ok((
        AdapterSniff {
            database_version: KNOWN_ADAPTERS_VERSION,
            r1,
            r2,
        },
        plan,
    ))
}

/// Pass 1: rough support per library entry, plus the positional denominators.
fn scan_pass(
    input: &CbqInput,
    work: &[BlockWork],
    plan: &SamplePlan,
    threads: usize,
    params: Params,
    library: &[SeedIndex; 2],
) -> Result<[Scan; 2]> {
    let entries = library[0].adapters().len();
    let min_overlap = params.matcher.min_overlap;
    let paired = input.schema().paired;
    // The existing inference, with its existing thresholds. Sniffing does not
    // introduce a second overlap algorithm or a second set of thresholds.
    let overlap_params = crate::overlap::OverlapParams::default();
    let states = for_each_sampled(
        input,
        work,
        plan,
        threads,
        || ScanState {
            r1: Scan::new(entries, true),
            r2: Scan::new(entries, paired),
        },
        |state, r1, r2| {
            for (mate, sequence) in mates(r1, r2) {
                let counts = state.mate(mate);
                counts.sampled += 1;
                let (tail, body) = positions(sequence.len(), min_overlap);
                counts.tail_positions += tail;
                counts.body_positions += body;
                if artifact(sequence, b'A') {
                    counts.poly_a += 1;
                }
                if artifact(sequence, b'G') {
                    counts.poly_g += 1;
                }
                counts.count_tail_kmers(sequence);
                if let Some((entry, _)) = library[mate.index()].find(sequence, params.matcher) {
                    counts.hits[entry] += 1;
                    counts.informative += 1;
                }
            }

            // One overlap analysis per pair, shared by both mates: the inferred
            // insert boundary names where read-through begins, and everything
            // past it is adapter by construction rather than by matching.
            if let Some(second) = r2
                && let Some(overlap) = crate::overlap::find_overlap(r1, second, overlap_params)
            {
                for (mate, sequence) in mates(r1, r2) {
                    let boundary = crate::overlap::insert_boundary(overlap, sequence.len());
                    if sequence.len() >= boundary + KMER_K {
                        state
                            .mate(mate)
                            .overhang
                            .observe_prefix(&sequence[boundary..]);
                    }
                }
            }
        },
    )?;

    let mut r1 = Scan::new(entries, true);
    let mut r2 = Scan::new(entries, paired);
    for state in &states {
        r1.merge(&state.r1);
        r2.merge(&state.r2);
    }
    Ok([r1, r2])
}

/// Pass 2: the single source of every support number.
fn verify_pass(
    input: &CbqInput,
    work: &[BlockWork],
    plan: &SamplePlan,
    threads: usize,
    params: Params,
    candidates: &[Vec<Adapter>; 2],
) -> Result<[Vec<Verified>; 2]> {
    let states = for_each_sampled(
        input,
        work,
        plan,
        threads,
        || VerifyState {
            r1: vec![Verified::new(); candidates[0].len()],
            r2: vec![Verified::new(); candidates[1].len()],
        },
        |state, r1, r2| {
            for (mate, sequence) in mates(r1, r2) {
                let tallies = match mate {
                    Mate::R1 => &mut state.r1,
                    Mate::R2 => &mut state.r2,
                };
                for (tally, adapter) in tallies.iter_mut().zip(&candidates[mate.index()]) {
                    if let Some(hit) = best_match(adapter, params.matcher, sequence) {
                        tally.record(&hit, sequence.len(), adapter.len());
                    }
                }
            }
        },
    )?;

    let mut r1 = vec![Verified::new(); candidates[0].len()];
    let mut r2 = vec![Verified::new(); candidates[1].len()];
    for state in &states {
        for (slot, add) in r1.iter_mut().zip(&state.r1) {
            slot.merge(add);
        }
        for (slot, add) in r2.iter_mut().zip(&state.r2) {
            slot.merge(add);
        }
    }
    Ok([r1, r2])
}

/// A candidate sequence awaiting verification.
struct Proposed {
    adapter: Adapter,
    known_name: Option<String>,
    known_category: Option<String>,
    sources: Vec<EvidenceSource>,
    /// Declaration order for a library entry; `usize::MAX` for a de novo one, so
    /// a named entry wins an otherwise exact tie.
    library_index: usize,
    /// Pairs whose overlap overhang supported this sequence. Provenance, not
    /// support: the verification pass measures support.
    overlap_support: u64,
}

/// Collects the candidates worth verifying, deduplicated.
///
/// A read carrying a known adapter also carries its k-mers, so the same sequence
/// usually arrives from both evidence sources. Merging them here — rather than
/// reporting two candidates, or summing their support — is what keeps a single
/// library from looking like a mixed one, and is why support is measured once
/// afterwards by verification.
fn propose(scan: &Scan, library: &SeedIndex, consensus: &[Consensus]) -> Vec<Proposed> {
    let mut proposals: Vec<Proposed> = scan
        .shortlist()
        .into_iter()
        .map(|entry| Proposed {
            adapter: library.adapters()[entry].clone(),
            known_name: Some(library.adapters()[entry].name.clone()),
            known_category: Some(KNOWN_ADAPTERS[entry].category().to_string()),
            sources: vec![EvidenceSource::KnownDatabase],
            library_index: entry,
            overlap_support: 0,
        })
        .collect();

    for columns in consensus {
        add_evidence(
            &mut proposals,
            &columns.sequence(),
            EvidenceSource::KmerConsensus,
            0,
        );
    }
    // Overlap overhangs are the strongest evidence available on paired input:
    // the insert boundary is inferred from the mates themselves, so everything
    // past it is adapter without any sequence having to be matched.
    add_evidence(
        &mut proposals,
        &scan.overhang.sequence(),
        EvidenceSource::PairedOverlap,
        scan.overhang.voting,
    );
    proposals
}

/// Folds one assembled sequence into the proposal set.
///
/// A sequence equivalent to something already proposed records the extra
/// evidence source and keeps the existing entry, whose name is more useful to a
/// reader; support is never summed across sources.
fn add_evidence(
    proposals: &mut Vec<Proposed>,
    sequence: &[u8],
    source: EvidenceSource,
    overlap_support: u64,
) {
    if sequence.len() < KMER_K {
        return;
    }
    if let Some(existing) = proposals
        .iter_mut()
        .find(|proposal| related(&proposal.adapter.sequence, sequence))
    {
        if !existing.sources.contains(&source) {
            existing.sources.push(source);
        }
        existing.overlap_support = existing.overlap_support.max(overlap_support);
        return;
    }
    let Ok(adapter) = Adapter::new("de novo consensus", sequence) else {
        return;
    };
    proposals.push(Proposed {
        adapter,
        known_name: None,
        known_category: None,
        sources: vec![source],
        library_index: usize::MAX,
        overlap_support,
    });
}

/// The earliest accepted match of one adapter in a read.
///
/// The same acceptance rule the trimming matcher uses, restricted to one
/// adapter, so a candidate's measured support is the support it would have as a
/// configured adapter.
fn best_match(adapter: &Adapter, params: AdapterParams, sequence: &[u8]) -> Option<AdapterHit> {
    crate::adapter::find_three_prime(std::slice::from_ref(adapter), params, sequence)
        .or_else(|| verify_at(adapter, params, sequence, 0))
}

/// Deterministic candidate ordering.
///
/// Every key is derived from the data, so the order never depends on how the
/// sample was traversed or on how many threads ran.
fn order(a: &Candidate, b: &Candidate) -> std::cmp::Ordering {
    b.confidence
        .cmp(&a.confidence)
        .then(b.supporting_reads.cmp(&a.supporting_reads))
        .then_with(|| {
            b.tail_enrichment
                .unwrap_or(f64::INFINITY)
                .partial_cmp(&a.tail_enrichment.unwrap_or(f64::INFINITY))
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .then(
            b.tail_connected_fraction
                .partial_cmp(&a.tail_connected_fraction)
                .unwrap_or(std::cmp::Ordering::Equal),
        )
        .then(b.sequence.len().cmp(&a.sequence.len()))
        .then(
            a.mean_error_rate
                .partial_cmp(&b.mean_error_rate)
                .unwrap_or(std::cmp::Ordering::Equal),
        )
        .then(a.library_index.cmp(&b.library_index))
        .then(a.sequence.cmp(&b.sequence))
}

/// Turns tallies into classified, ordered candidates and a decision.
fn assemble(
    mate: Mate,
    scan: &Scan,
    proposals: &[Proposed],
    tallies: &[Verified],
    params: Params,
) -> MateResult {
    let gates = params.gates;
    let mut candidates: Vec<Candidate> = proposals
        .iter()
        .zip(tallies)
        .map(|(proposal, tally)| build(proposal, tally, scan, gates))
        .collect();

    candidates.sort_by(order);
    let decision = decide(&mut candidates, gates);
    // Deciding demotes the losers of a family, so re-order: the recommended
    // candidate must be the one a reader sees first.
    candidates.sort_by(order);
    let recommended = decision
        .is_confident()
        .then(|| candidates.iter().find(|c| c.confidence == Confidence::High))
        .flatten();
    let recommended_sequence =
        recommended.map(|c| String::from_utf8_lossy(&c.sequence).into_owned());
    let recommended_name = recommended.and_then(|c| c.known_name.clone());

    candidates.truncate(params.top);
    MateResult {
        mate: mate.name(),
        decision,
        recommended_sequence,
        recommended_name,
        candidates,
        sampled_reads: scan.sampled,
        informative_reads: scan.informative,
        poly_a_signal: scan.poly_a,
        poly_g_signal: scan.poly_g,
    }
}

fn build(proposal: &Proposed, tally: &Verified, scan: &Scan, gates: Gates) -> Candidate {
    let adapter = &proposal.adapter;
    let supporting = tally.supporting;
    let support_fraction = fraction(supporting, scan.sampled);
    let tail_connected_fraction = fraction(tally.tail_connected, supporting);
    // Normalized by searchable positions, so an abundant sequence that happens
    // to be everywhere does not look 3'-specific just because reads are long.
    let tail_rate = fraction(tally.tail_hits, scan.tail_positions);
    let body_rate = fraction(tally.body_hits, scan.body_positions);
    let tail_enrichment = if body_rate > 0.0 {
        Some(tail_rate / body_rate)
    } else if tail_rate > 0.0 {
        None
    } else {
        Some(0.0)
    };
    let mean_error_rate = if tally.overlap_sum == 0 {
        0.0
    } else {
        tally.error_sum as f64 / tally.overlap_sum as f64
    };

    let mut candidate = Candidate {
        sequence: adapter.sequence.clone(),
        known_name: proposal.known_name.clone(),
        known_category: proposal.known_category.clone(),
        confidence: Confidence::Low,
        evidence_sources: proposal.sources.clone(),
        supporting_reads: supporting,
        support_fraction,
        tail_connected_fraction,
        tail_enrichment,
        median_start: median(&tally.starts, supporting),
        median_distance_to_end: median(&tally.distances, supporting),
        exact_matches: tally.exact,
        substitution_matches: tally.substitution,
        indel_matches: tally.indel,
        mean_error_rate,
        paired_overlap_support: proposal.overlap_support,
        library_index: proposal.library_index,
    };
    candidate.confidence = classify(&candidate, scan, gates);
    candidate
}

/// Applies the classification gates.
///
/// A database name proves sequence identity, not that trimming from that
/// coordinate is safe. Known candidates therefore need corroboration: reaching
/// the read end, enough inferred pair boundaries, 3' localization, or a stable
/// consensus assembled independently from read tails. Unnamed candidates retain
/// the consensus-stability and tail-enrichment gates that establish provenance.
fn classify(candidate: &Candidate, scan: &Scan, gates: Gates) -> Confidence {
    let enough_sample = scan.sampled >= gates.min_sample;
    let enough_reads = candidate.supporting_reads >= gates.min_supporting_reads;
    let long_enough = candidate.sequence.len() >= gates.min_candidate_length;
    if !(enough_sample && enough_reads && long_enough) {
        return Confidence::Low;
    }
    let complex = complex_enough(&candidate.sequence, gates);
    let clean = candidate.mean_error_rate <= gates.max_mean_error_rate;
    if !(complex && clean) {
        return Confidence::Low;
    }
    let localized = candidate.tail_connected_fraction >= gates.min_tail_connected_fraction;
    let enriched = candidate
        .tail_enrichment
        .is_none_or(|value| value >= gates.min_tail_enrichment);
    let end_gap = candidate
        .median_distance_to_end
        .saturating_sub(candidate.sequence.len() as u64);
    let end_reaching = end_gap <= END_REACH_SLACK;
    let paired = candidate.paired_overlap_support >= gates.min_supporting_reads;
    let known = candidate.known_name.is_some();
    let assembled = candidate
        .evidence_sources
        .iter()
        .any(|source| matches!(source, EvidenceSource::KmerConsensus))
        || (paired
            && candidate
                .evidence_sources
                .contains(&EvidenceSource::PairedOverlap));
    let extended = assembled && candidate.sequence.len() >= gates.min_consensus_length;
    let trimming_safe = if known {
        paired || end_reaching || localized || extended
    } else {
        extended || (localized && enriched)
    };
    if candidate.support_fraction >= gates.min_support_fraction && trimming_safe {
        Confidence::High
    } else if trimming_safe || localized || enriched || known {
        Confidence::Medium
    } else {
        Confidence::Low
    }
}

/// Sequence complexity gate: no dominant base, and enough adjacent variation.
fn complex_enough(sequence: &[u8], gates: Gates) -> bool {
    let mut counts = [0usize; 256];
    for &base in sequence {
        counts[base as usize] += 1;
    }
    let dominant = counts.iter().max().copied().unwrap_or(0);
    if dominant as f64 > gates.max_base_share * sequence.len() as f64 {
        return false;
    }
    crate::filter::complexity(sequence) >= gates.min_complexity
}

/// Resolves the mate-level decision, demoting candidates that lose.
///
/// Two unrelated high-confidence candidates mean a mixed library, and the plan
/// forbids resolving that silently. A much weaker unrelated candidate is not a
/// second library, though: below `competitor_share` of the leader's support it
/// is reported as evidence but demoted, so ordinary trace contamination does not
/// make every file `mixed`.
fn decide(candidates: &mut [Candidate], gates: Gates) -> Decision {
    let Some(first) = candidates
        .iter()
        .position(|c| c.confidence == Confidence::High)
    else {
        return Decision::Inconclusive;
    };

    // Which spelling of the leading family to report. Support alone is the wrong
    // key: the library holds the same adapter at several frame shifts, and a
    // sequence one base shorter matches marginally *more* reads while leaving
    // that base behind on every trimmed read. Measured on a 200 000-pair overlap
    // fixture, the truncated spelling led the correct one 70 933 to 70 638 —
    // 0.4% — while the correct one carried three independent evidence sources to
    // the other's one. Corroboration decides it, which is the whole premise of a
    // composite detector. This applies only *within* a family, so a weakly
    // supported candidate can never displace an unrelated stronger one.
    let head = candidates[first].sequence.clone();
    let leader = candidates
        .iter()
        .enumerate()
        .filter(|(_, candidate)| {
            candidate.confidence == Confidence::High && related(&head, &candidate.sequence)
        })
        .min_by_key(|(index, candidate)| {
            (
                std::cmp::Reverse(candidate.evidence_sources.len()),
                std::cmp::Reverse(candidate.sequence.len()),
                std::cmp::Reverse(candidate.supporting_reads),
                *index,
            )
        })
        .map_or(first, |(index, _)| index);

    let leader_support = candidates[leader].supporting_reads;
    let leader_sequence = candidates[leader].sequence.clone();
    let mut competitors = 0usize;
    for (index, candidate) in candidates.iter_mut().enumerate() {
        if index == leader || candidate.confidence != Confidence::High {
            continue;
        }
        if related(&leader_sequence, &candidate.sequence) {
            // Same family: a prefix variant of the winner, not a rival.
            candidate.confidence = Confidence::Medium;
            continue;
        }
        let share = candidate.supporting_reads as f64 / leader_support.max(1) as f64;
        if share >= gates.competitor_share {
            competitors += 1;
        } else {
            candidate.confidence = Confidence::Medium;
        }
    }
    if competitors > 0 {
        Decision::Mixed
    } else {
        Decision::Confident
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn candidate(sequence: &[u8], supporting: u64, known: bool) -> Candidate {
        Candidate {
            sequence: sequence.to_vec(),
            known_name: known.then(|| "test".to_string()),
            known_category: known.then(|| "adapter".to_string()),
            confidence: Confidence::High,
            evidence_sources: vec![EvidenceSource::KnownDatabase],
            supporting_reads: supporting,
            support_fraction: 0.5,
            tail_connected_fraction: 1.0,
            tail_enrichment: Some(100.0),
            median_start: 0,
            median_distance_to_end: 0,
            exact_matches: supporting,
            substitution_matches: 0,
            indel_matches: 0,
            mean_error_rate: 0.0,
            paired_overlap_support: 0,
            library_index: 0,
        }
    }

    #[test]
    fn one_high_candidate_is_confident() {
        let mut candidates = vec![candidate(b"AGATCGGAAGAGCACACG", 1000, true)];
        assert_eq!(
            decide(&mut candidates, Gates::default()),
            Decision::Confident
        );
    }

    #[test]
    fn two_unrelated_high_candidates_are_mixed() {
        let mut candidates = vec![
            candidate(b"AGATCGGAAGAGCACACG", 1000, true),
            candidate(b"CTGTCTCTTATACACATCT", 900, true),
        ];
        assert_eq!(decide(&mut candidates, Gates::default()), Decision::Mixed);
    }

    #[test]
    fn prefix_variants_are_one_family_not_a_mixed_library() {
        let mut candidates = vec![
            candidate(b"AGATCGGAAGAGCACACG", 1000, true),
            candidate(b"AGATCGGAAGAGCACACGTCTGAACT", 950, true),
        ];
        assert_eq!(
            decide(&mut candidates, Gates::default()),
            Decision::Confident
        );
        // Exactly one survives, and it is the longer spelling even though the
        // truncated one matched fifty more reads: trimming with the short one
        // would leave eight adapter bases behind on every read.
        let high: Vec<&Candidate> = candidates
            .iter()
            .filter(|c| c.confidence == Confidence::High)
            .collect();
        assert_eq!(high.len(), 1);
        assert_eq!(high[0].sequence, b"AGATCGGAAGAGCACACGTCTGAACT");
    }

    #[test]
    fn corroboration_decides_between_two_spellings_of_one_adapter() {
        // The real case from the overlap fixture: a frame-shifted spelling with
        // marginally more support, against the true adapter with three sources.
        let mut truncated = candidate(b"GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT", 70_933, true);
        truncated.evidence_sources = vec![EvidenceSource::KnownDatabase];
        let mut correct = candidate(b"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT", 70_638, true);
        correct.evidence_sources = vec![
            EvidenceSource::KnownDatabase,
            EvidenceSource::KmerConsensus,
            EvidenceSource::PairedOverlap,
        ];
        let mut candidates = vec![truncated, correct];
        assert_eq!(
            decide(&mut candidates, Gates::default()),
            Decision::Confident
        );
        let leader = candidates
            .iter()
            .find(|c| c.confidence == Confidence::High)
            .unwrap();
        assert_eq!(leader.sequence, b"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT");
    }

    #[test]
    fn a_weak_unrelated_candidate_does_not_make_a_library_mixed() {
        let mut candidates = vec![
            candidate(b"AGATCGGAAGAGCACACG", 1000, true),
            candidate(b"CTGTCTCTTATACACATCT", 10, true),
        ];
        assert_eq!(
            decide(&mut candidates, Gates::default()),
            Decision::Confident
        );
        assert_eq!(candidates[1].confidence, Confidence::Medium);
    }

    #[test]
    fn no_high_candidate_is_inconclusive() {
        let mut candidates = vec![candidate(b"AGATCGGAAGAGCACACG", 1000, true)];
        candidates[0].confidence = Confidence::Medium;
        assert_eq!(
            decide(&mut candidates, Gates::default()),
            Decision::Inconclusive
        );
    }

    #[test]
    fn a_small_sample_cannot_produce_a_recommendation() {
        let scan = Scan {
            sampled: 10,
            ..Scan::new(0, false)
        };
        let candidate = candidate(b"AGATCGGAAGAGCACACG", 5, true);
        assert_eq!(
            classify(&candidate, &scan, Gates::default()),
            Confidence::Low
        );
    }

    #[test]
    fn a_terminal_known_adapter_without_enrichment_is_still_high() {
        // A full adapter that reaches the read end is safe even when the
        // normalized enrichment statistic contributes no evidence.
        let scan = Scan {
            sampled: 10_000,
            ..Scan::new(0, false)
        };
        let mut candidate = candidate(b"AGATCGGAAGAGCACACG", 5_000, true);
        candidate.median_start = 30;
        candidate.median_distance_to_end = candidate.sequence.len() as u64;
        candidate.tail_connected_fraction = 0.0;
        candidate.tail_enrichment = Some(0.0);
        assert_eq!(
            classify(&candidate, &scan, Gates::default()),
            Confidence::High
        );
    }

    #[test]
    fn a_start_only_known_sequence_is_not_high_confidence() {
        let scan = Scan {
            sampled: 10_000,
            ..Scan::new(0, false)
        };
        let mut candidate = candidate(b"AATGATACGGCGACCACCGACAGGTTCAGAGT", 5_000, true);
        candidate.median_start = 0;
        candidate.median_distance_to_end = candidate.sequence.len() as u64 + 100;
        candidate.tail_connected_fraction = 0.0;
        candidate.tail_enrichment = Some(0.0);
        assert_eq!(
            classify(&candidate, &scan, Gates::default()),
            Confidence::Medium
        );
    }

    #[test]
    fn paired_overlap_needs_the_absolute_support_floor() {
        let scan = Scan {
            sampled: 10_000,
            ..Scan::new(0, false)
        };
        let gates = Gates::default();
        let mut candidate = candidate(b"AGATCGGAAGAGCACACG", 5_000, true);
        candidate
            .evidence_sources
            .push(EvidenceSource::PairedOverlap);
        candidate.median_start = 30;
        candidate.median_distance_to_end = candidate.sequence.len() as u64 + 100;
        candidate.tail_connected_fraction = 0.0;
        candidate.tail_enrichment = Some(0.0);

        candidate.paired_overlap_support = 1;
        assert_eq!(classify(&candidate, &scan, gates), Confidence::Medium);

        candidate.paired_overlap_support = gates.min_supporting_reads;
        assert_eq!(classify(&candidate, &scan, gates), Confidence::High);
    }

    #[test]
    fn a_homopolymer_is_never_recommended() {
        let scan = Scan {
            sampled: 10_000,
            ..Scan::new(0, false)
        };
        let candidate = candidate(b"AAAAAAAAAAAAAAAAAA", 5_000, false);
        assert_eq!(
            classify(&candidate, &scan, Gates::default()),
            Confidence::Low
        );
    }

    #[test]
    fn medians_come_from_the_histogram() {
        let mut tally = Verified::new();
        for start in [10usize, 20, 20, 30] {
            tally.starts[start] += 1;
        }
        assert_eq!(median(&tally.starts, 4), 20);
        assert_eq!(median(&tally.starts, 0), 0);
    }

    #[test]
    fn searchable_positions_split_at_the_tail_window() {
        // 150 bases, min overlap 8: 143 start coordinates, 60 of them in the tail.
        assert_eq!(positions(150, 8), (60, 83));
        // A read shorter than the window is entirely tail.
        assert_eq!(positions(40, 8), (33, 0));
        // A read shorter than the minimum overlap has none.
        assert_eq!(positions(4, 8), (0, 0));
    }

    #[test]
    fn a_kmer_counts_once_per_read_however_often_it_repeats() {
        let mut scan = Scan::new(0, true);
        // A tandem repeat of one 10-mer, filling the tail window.
        let repeat = b"ACGTACGGTT".repeat(6);
        scan.count_tail_kmers(&repeat);
        let code = crate::detect::encode_kmer(b"ACGTACGGTT").unwrap();
        assert_eq!(
            scan.kmers[code as usize], 1,
            "a repeat must not out-vote a single occurrence"
        );

        // A second read carrying it once adds exactly one more.
        let mut once = vec![b'C'; 50];
        once.extend_from_slice(b"ACGTACGGTT");
        scan.count_tail_kmers(&once);
        assert_eq!(scan.kmers[code as usize], 2);
    }

    #[test]
    fn ambiguous_bases_invalidate_the_windows_they_touch() {
        let mut clean = Scan::new(0, true);
        clean.count_tail_kmers(b"ACGTACGGTTCCAA");
        let mut dirty = Scan::new(0, true);
        dirty.count_tail_kmers(b"ACGTACGNTTCCAA");
        let code = crate::detect::encode_kmer(b"ACGTACGGTT").unwrap();
        assert_eq!(clean.kmers[code as usize], 1);
        assert_eq!(dirty.kmers[code as usize], 0);
        // Windows clear of the `N` still count.
        let after = crate::detect::encode_kmer(b"GNTTCCAA").is_none();
        assert!(after, "a k-mer spanning the N cannot be encoded");
    }

    #[test]
    fn only_the_tail_window_is_counted() {
        let mut scan = Scan::new(0, true);
        let mut read = b"ACGTACGGTT".to_vec();
        read.extend(std::iter::repeat_n(b'C', TAIL_WINDOW));
        scan.count_tail_kmers(&read);
        let code = crate::detect::encode_kmer(b"ACGTACGGTT").unwrap();
        assert_eq!(
            scan.kmers[code as usize], 0,
            "a k-mer outside the tail window is not adapter evidence"
        );
    }

    #[test]
    fn artifact_and_low_complexity_seeds_are_never_extended() {
        let gates = Gates::default();
        let poly_a = crate::detect::encode_kmer(b"AAAAAAAAAA").unwrap();
        let poly_g = crate::detect::encode_kmer(b"GGGGGGGGGG").unwrap();
        let dinucleotide = crate::detect::encode_kmer(b"ATATATATAT").unwrap();
        let trinucleotide = crate::detect::encode_kmer(b"ACGACGACGA").unwrap();
        let real = crate::detect::encode_kmer(b"AGATCGGAAG").unwrap();
        assert!(!informative_seed(poly_a, gates));
        assert!(!informative_seed(poly_g, gates));
        // Scores a perfect 1.0 on adjacent-base complexity, so only the
        // periodicity test rejects it.
        assert!(!informative_seed(dinucleotide, gates));
        assert!(!informative_seed(trinucleotide, gates));
        assert!(informative_seed(real, gates));
    }

    #[test]
    fn periodicity_only_rejects_exact_repeats() {
        assert!(periodic(b"AAAAAAAAAA", 3));
        assert!(periodic(b"ATATATATAT", 3));
        assert!(periodic(b"ACGACGACGA", 3));
        // One base out of phase is no longer an exact repeat.
        assert!(!periodic(b"ATATATATAG", 3));
        assert!(!periodic(b"AGATCGGAAG", 3));
        // Period four is left to the other tests.
        assert!(!periodic(b"ACGTACGTAC", 3));
    }

    #[test]
    fn the_seed_shortlist_is_ordered_and_bounded() {
        let mut scan = Scan::new(0, true);
        let gates = Gates::default();
        // Three informative seeds with distinct support.
        for (kmer, count) in [
            (b"AGATCGGAAG", 500u32),
            (b"CCTGAACTCC", 900),
            (b"GTCGATCGTA", 700),
        ] {
            scan.kmers[crate::detect::encode_kmer(kmer).unwrap() as usize] = count;
        }
        // Below the floor, so excluded however informative.
        scan.kmers[crate::detect::encode_kmer(b"TTCAGACGTG").unwrap() as usize] = 3;
        let seeds = scan.seeds(gates);
        let sequences: Vec<Vec<u8>> = seeds
            .iter()
            .map(|&(code, _)| crate::detect::decode_kmer(code, KMER_K))
            .collect();
        assert_eq!(sequences.len(), 3);
        assert_eq!(sequences[0], b"CCTGAACTCC");
        assert_eq!(sequences[1], b"GTCGATCGTA");
        assert_eq!(sequences[2], b"AGATCGGAAG");
    }

    #[test]
    fn a_consensus_extends_through_agreement_and_stops_at_disagreement() {
        let adapter = b"GTCGATCGTACGGCATCCGATCGTACGATCGG";
        let mut consensus = Consensus::new();
        // Ten reads whose inserts differ but whose tails agree.
        for insert in 0..10usize {
            let mut read = vec![b'A'; 12 + insert];
            read.extend_from_slice(adapter);
            // The seed sits at the adapter's start.
            consensus.observe(&read, 12 + insert);
        }
        let sequence = consensus.sequence();
        assert!(
            sequence.ends_with(adapter),
            "extension stopped inside the adapter: {}",
            String::from_utf8_lossy(&sequence)
        );
        assert!(sequence.len() <= MAX_CONSENSUS);
    }

    #[test]
    fn a_consensus_does_not_grow_past_its_evidence() {
        // The seed is shared but what follows differs in every read, as for a
        // genomic repeat. Extension must stop almost immediately.
        let mut consensus = Consensus::new();
        for (index, tail) in [
            b"ACGTACGGTTAAAAAAAA".as_slice(),
            b"ACGTACGGTTCCCCCCCC".as_slice(),
            b"ACGTACGGTTGGGGGGGG".as_slice(),
            b"ACGTACGGTTTTTTTTTT".as_slice(),
        ]
        .into_iter()
        .enumerate()
        {
            let mut read = vec![b'C'; 5 + index];
            read.extend_from_slice(tail);
            consensus.observe(&read, 5 + index);
        }
        let sequence = consensus.sequence();
        assert!(
            sequence.len() < Gates::default().min_consensus_length,
            "a repeat with varied flanks assembled {} bases: {}",
            sequence.len(),
            String::from_utf8_lossy(&sequence)
        );
    }

    #[test]
    fn merging_consensus_columns_is_order_independent() {
        let adapter = b"GTCGATCGTACGGCATCCGATCGTACGATCGG";
        let build = |inserts: &[usize]| {
            let mut consensus = Consensus::new();
            for &insert in inserts {
                let mut read = vec![b'A'; insert];
                read.extend_from_slice(adapter);
                consensus.observe(&read, insert);
            }
            consensus
        };
        let mut forward = build(&[10, 11, 12]);
        forward.merge(&build(&[13, 14]));
        let mut backward = build(&[13, 14]);
        backward.merge(&build(&[10, 11, 12]));
        assert_eq!(forward.voting, backward.voting);
        assert_eq!(forward.sequence(), backward.sequence());
    }

    #[test]
    fn a_de_novo_candidate_needs_a_stable_extension_to_be_recommended() {
        let scan = Scan {
            sampled: 10_000,
            ..Scan::new(0, false)
        };
        let gates = Gates::default();
        let mut candidate = candidate(b"GTCGATCGTACGGCATCCGATCGTACGATCGG", 5_000, false);
        candidate.evidence_sources = vec![EvidenceSource::KmerConsensus];
        // Read-through covers body and tail alike, so neither positional test
        // fires; the stable extension is what carries the unnamed candidate.
        candidate.tail_connected_fraction = 0.27;
        candidate.tail_enrichment = Some(1.0);
        assert_eq!(classify(&candidate, &scan, gates), Confidence::High);

        // The same evidence behind a short consensus is not enough.
        let mut short = candidate.clone();
        short.sequence = b"GTCGATCGTACGGC".to_vec();
        assert_eq!(classify(&short, &scan, gates), Confidence::Low);
    }

    #[test]
    fn stable_tail_assembly_corroborates_a_known_adapter_with_long_readthrough() {
        let scan = Scan {
            sampled: 10_000,
            ..Scan::new(0, false)
        };
        let gates = Gates::default();
        let mut candidate = candidate(b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA", 5_000, true);
        candidate.evidence_sources =
            vec![EvidenceSource::KnownDatabase, EvidenceSource::KmerConsensus];
        candidate.tail_connected_fraction = 0.27;
        candidate.tail_enrichment = Some(1.0);
        candidate.median_start = 80;
        candidate.median_distance_to_end = 70;

        assert_eq!(classify(&candidate, &scan, gates), Confidence::High);
    }

    #[test]
    fn an_internal_known_sequence_is_not_recommended_without_end_evidence() {
        let scan = Scan {
            sampled: 4_000,
            ..Scan::new(0, false)
        };
        let mut candidate = candidate(b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA", 4_000, true);
        candidate.median_start = 30;
        candidate.median_distance_to_end = 70;
        candidate.tail_connected_fraction = 0.0;
        candidate.tail_enrichment = Some(1.0);
        assert_eq!(
            classify(&candidate, &scan, Gates::default()),
            Confidence::Medium
        );
    }

    #[test]
    fn artifact_windows_are_recognised() {
        let poly_a = vec![b'A'; 150];
        assert!(artifact(&poly_a, b'A'));
        assert!(!artifact(&poly_a, b'G'));
        let mixed = b"ACGT".repeat(40);
        assert!(!artifact(&mixed, b'A'));
        // Too short to judge.
        assert!(!artifact(b"AAAA", b'A'));
    }
}
