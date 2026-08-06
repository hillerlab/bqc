# bqc

CBQ-native adapter removal, read trimming and per-read filtering.

`bqc` reads a [CBQ](https://docs.rs/binseq) file, transforms records in a
single pass and writes CBQ back out. There is no FASTQ conversion, no
intermediate file between stages, and no reordering of records.

```bash
bqc adapter   reads.cbq -o out.cbq --adapter-r1 AGATCGG...
bqc adapter   reads.cbq -o out.cbq --auto-detect          # infer the adapter
bqc trim      reads.cbq -o out.cbq --quality-tail 20 --trim-terminal-n
bqc filter    reads.cbq -o out.cbq --min-length 30 --failed rejected.cbq
bqc correct   pairs.cbq -o out.cbq --correction-log corrections.tsv
bqc workflow  reads.cbq -o out.cbq --config illumina.toml -T 8
bqc sniff adapters reads.cbq                              # inspect, never modify
bqc sniff strand   reads.cbq --index salmon-index         # RNA-seq orientation
```

## Scope

| Capability                                    | Where it lives     |
| --------------------------------------------- | ------------------ |
| FASTQ/BQ/VBQ ↔ CBQ conversion, concatenation, splitting, inspection | `bqtools` |
| Dataset-level QC reports                      | `bqtools qc`       |
| Adapter removal, with indels, overlap inference or auto-detection | `bqc adapter` |
| Positional, quality, N and poly-X trimming    | `bqc trim`      |
| Per-read and per-pair acceptance decisions    | `bqc filter`    |
| Paired-overlap base correction                | `bqc correct`   |
| Linked adapters: retain the insert between two flanks | `bqc adapter --linked-*` |
| Split reads at internal adapters into separate records | `bqc segment` |
| All of the above fused into one pass          | `bqc workflow`  |
| Non-destructive inspection: which adapters are in this file? | `bqc sniff adapters` |
| RNA-seq library strandedness, against a transcriptome | `bqc sniff strand` |

`bqtools qc` describes a dataset; `bqc filter` decides whether an individual
read stays. They are complementary:

```bash
bqtools qc raw.cbq -o qc-before
bqc workflow raw.cbq -o clean.cbq --config illumina.toml
bqtools qc clean.cbq -o qc-after
```

CBQ is the only supported input and output format. BQ, VBQ and FASTQ inputs are
rejected with a pointer to `bqtools encode`.

## Guarantees

* **Order preservation.** Output records appear in input order. One CBQ block is
  one unit of work, so the partitioning comes from the input file rather than
  from `--threads`. In practice the output is byte-for-byte identical for
  `-T 1`, `-T 8` and everything in between; identical *decoded* records are the
  hard guarantee.
* **Schema preservation.** Pairing, qualities, headers and flags are read from
  the input file header and reproduced exactly. `binseq` synthesizes record
  headers and quality values for files that do not store them; `bqc` never
  writes those synthetic values out, so a header-free file stays header-free.
* **Determinism.** Adapter tie-breaking, failure-reason ordering and stage order
  are fixed. The same input and configuration always produce the same output.
* **Atomic, collision-safe outputs.** Every file is written to an
  invocation-unique hidden temporary path in the destination directory and
  published only after the run succeeds. Without `--force`, publication is an
  atomic no-clobber operation, so a destination created concurrently is never
  overwritten. A failed run leaves no partial output.
* **Bounded memory.** Input is memory mapped; at most `3 × threads` blocks are
  in flight. Peak RSS is roughly the mapped file plus a few MB per thread.

## Processing order

The stage order is part of the output contract, not an implementation detail:

```text
correct → adapter → trim → filter
```

Correction comes first so that corrected bases *and* corrected qualities are
visible to every later decision — a base rescued from its mate can lift a read
back over a mean-quality threshold. Within `trim`:

```text
fixed front/tail
  → quality front → quality right → quality tail
  → terminal N
  → poly-G
  → poly-X
  → maximum-length truncation
```

This differs in detail from fastp's internal ordering; `bqc` does not claim
fastp parity. Changing the order would change output, so it is covered by
golden and property tests.

## Adapter matching

`bqc adapter` removes a 3' adapter. For every candidate coordinate `p` in the
read, the adapter prefix is compared against `read[p .. p + overlap]` where
`overlap = min(adapter_length, read_length - p)`. A candidate matches when

```text
overlap >= --min-overlap            (default 8)
errors  <= floor(overlap * --max-error-rate)   (default 0.10)
errors  <= --max-errors             (when supplied)
```

Among all matches the winner is chosen by: **earliest coordinate**, then fewest
errors, then longest overlap, then adapter declaration order. Trimming to the
earliest coordinate removes the longest contaminated suffix.

Notes:

* By default only substitutions are considered; see **Indel-aware matching**
  below for `--allow-indels`.
* `N` never matches — in the read or in the adapter — and counts as one error.
  This is deliberately conservative: a poly-N tail would otherwise match any
  adapter. Use a larger `--max-error-rate` if you want N-tolerant matching.
* A partial adapter only matches at the 3' end. A truncated adapter copy in the
  middle of a read is compared against the whole adapter, and therefore does not
  match.
* An adapter shorter than `--min-overlap` is rejected at startup rather than
  silently ignored.
* `--adapter-fasta` adapters apply to both mates; `--adapter-r1` and
  `--adapter-r2` are mate-specific. Mate-specific adapters come first in
  declaration order, then FASTA records in file order.

### Indel-aware matching

`--allow-indels` replaces the substitution scan with a banded edit-distance
alignment: insertions and deletions each cost one error and draw on the same
budget. The reported overlap is the number of aligned adapter bases, which with
indels can differ from the number of read bases consumed. Every
substitution-only alignment is also an edit-distance alignment of equal cost, so
nothing that matched before stops matching.

Selection still prefers the earliest coordinate, with one necessary refinement.
Edit distance casts *shadows*: an occurrence at `p` is usually also reachable at
`p - d` by consuming the `d` extra read bases as insertions, which costs at least
`d` errors. Taking the earliest match blindly therefore shaves a base or two off
otherwise clean matches. Since a shadow at distance `d` costs at least `d`
errors, the true occurrence lies within `errors` bases of the earliest accepted
match; that window is scanned and the candidate with the most matched adapter
bases wins, ties going to the earliest coordinate. Matches beyond the window are
separate occurrences, so a read carrying two adapter copies is still trimmed at
the first.

The alignment runs per candidate coordinate, which costs roughly an order of
magnitude more than substitution matching (see **Performance**). It is worth it
for chemistries that really do produce adapter indels; otherwise leave it off.

### Paired-overlap inference

`--paired-overlap` infers the insert length instead of matching a sequence. R1
and R2 are read from opposite ends of one insert, so the reverse complement of R2
must align to R1 wherever the mates genuinely overlap; a successful alignment
reveals the insert length, and everything past it in either mate is adapter
read-through. This trims adapters whose sequence you do not know.

Every placement is scored by Hamming distance over the aligned region (`N` never
matches) and accepted when `overlap >= --paired-overlap-min-overlap` (default 30)
and `mismatches <= floor(overlap * --max-error-rate)`. The winner is the longest
overlap, then fewest mismatches, then smallest offset: length wins first because
a short perfect alignment is cheap to fake, while a long nearly-perfect one is
almost always the true shared insert.

When an alignment succeeds it fixes the boundary for both mates and sequence
matching is skipped for that pair — the boundary is already known, and a match
inside the insert would be a false positive. Pairs with no acceptable alignment
fall back to explicit adapter matching.

### Auto-detection

`--auto-detect` infers the adapters from the data before the main pass. It is
the same detector `bqc sniff adapters` reports — there is one — invoked here
for its recommendation instead of for a report, so the two can never disagree
about what is in a file. See [Sniffing](#sniffing) for how the evidence is
gathered and weighed.

`--detect-sample-size` (default 262144) and `--detect-min-support` (default 0.01)
tune the sample and the support gate; `--span` restricts both.

Detection **refuses to trim** rather than guess, in two distinct situations:

* **Nothing clears the gates.** The run aborts with the evidence, asking for
  explicit sequences — unless `--paired-overlap` was also requested, which needs
  no adapter sequence and proceeds.
* **Two unrelated adapters clear them.** A mixed library is a fact about the
  data, not a tie to break; trimming half the reads with the wrong sequence is
  worse than stopping. Run `bqc sniff adapters` to see both candidates and
  choose explicitly.

Note that the mates are inferred independently. A mate with no adapter evidence
gets no adapter: assuming both mates share a chemistry would configure a
sequence nothing in the data supports.

Detected adapters are appended to any explicit ones and appear in the resolved
configuration, so the report always states what was actually trimmed. The full
sniff result — every candidate, its evidence sources and its statistics — is
embedded under `adapter.detection`.

### Linked adapters

A linked adapter is a *pair* of flanks with the sequence of interest between
them, which is how amplicon and many small-RNA libraries are built. Instead of
removing a 3' end, `--linked-*` retains the insert:

```text
read     [5' flank] INSERT [3' flank] read-through
retained            INSERT
```

```bash
bqc adapter amplicon.cbq -o inserts.cbq \
  --linked-5p-r1 AGGTCAGTCTAC --linked-3p-r1 CTTACGGATCCA \
  --linked-min-insert-length 20 --linked-unmatched fail
```

```text
    --linked-5p-r1 <SEQ>            5' flank on R1 (with --linked-3p-r1)
    --linked-3p-r1 <SEQ>            3' flank on R1
    --linked-5p-r2 <SEQ>            5' flank on R2 (paired input only)
    --linked-3p-r2 <SEQ>            3' flank on R2
    --linked-require <REQUIRE>      both (default) | either
    --linked-max-5p-offset <INT>    bases the 5' flank may start within [3]
    --linked-max-3p-overhang <INT>  3' flank bases that may hang past the read end
    --linked-min-insert-length <INT>  shortest insert retained [1]
    --linked-unmatched <POLICY>     continue (default) | keep | fail
```

The 5' flank is anchored: it must begin within `--linked-max-5p-offset` bases of
the read's 5' end and be consumed whole, unless the read's own start cuts it
short. The 3' flank is found by the ordinary 3' matcher, and must begin at or
after the end of the 5' flank — a 3' flank inside the 5' flank is not a linked
structure. Both flanks share the matcher thresholds (`--min-overlap`,
`--max-error-rate`, `--max-errors`, `--allow-indels`).

Candidate definitions are compared by: both flanks before one flank, then fewer
errors, then lower error rate, then longer flank overlap, then a longer retained
insert, then declaration order. A definition is rejected when the flanks appear
out of order or the insert falls below the minimum.

`--linked-unmatched` decides the rest of the library:

| Policy | Effect on a read with no linked match |
| --- | --- |
| `continue` (default) | ordinary adapter matching still runs on it |
| `keep` | written unchanged; adapter matching is skipped |
| `fail` | rejected with the reason `LINKED_UNMATCHED` |

A successful linked match fixes both boundaries, so ordinary adapter matching is
skipped for that read — the same precedent paired-overlap inference sets.

`--linked-5p-r1` and `--linked-3p-r1` declare one definition per mate; a
configuration file can declare several, each with its own name and thresholds:

```toml
[[adapter.linked_r1]]
name = "amplicon-a"
five_prime = "AGGTCAGTCTAC"
three_prime = "CTTACGGATCCA"
minimum_insert_length = 20

[[adapter.linked_r1]]
name = "amplicon-b"
five_prime = "TTCGCAGTCAGT"
three_prime = "GGATCCTAAGCC"
```

## Sniffing

`bqc sniff` inspects a file without changing it. It never trims, filters,
reorders or rewrites anything; the input is opened read-only and is
byte-identical afterwards.

```bash
bqc sniff adapters reads.cbq                       # human-readable summary
bqc sniff adapters reads.cbq --format json -o adapters.json
bqc sniff adapters reads.cbq --require-confident \
    --emit-config adapters.toml                       # for a pipeline
```

### Distributed sampling, not the first N reads

Sniffing inspects a subset, and *which* subset matters. Scanning the leading
records — what a FASTQ tool has to do — is biased by everything that puts
unusual data at the start of a file: concatenated lanes, tile-ordered records,
joined runs, calibration reads. CBQ is indexed and block-addressable, so the
sample is spread evenly across the file instead. For `S` records out of `N`:

```text
index(i) = (2i + 1) * N / (2S)      for i in 0..S
```

evaluated in integers. There is no generator and no seed, so nothing has to be
recorded to reproduce a run, and the answer is identical at every thread count.
`--sample-size` sets `S` (default 262144) and `--span START..END` restricts the
range.

### What counts as evidence

Three sources, combined:

| Source | What it is |
| --- | --- |
| `known_database` | A versioned library of 234 published adapter, primer and PhiX sequences. A fixed 5' 12-mer handles the common case; a lossless q-gram index covers accepted errors and terminal partials, and unusually permissive settings fall back to exhaustive matching. Every candidate is verified with the same matcher that trims. |
| `kmer_consensus` | Exact 10-mers of each read's final 60 bases, counted **once per read** so a tandem repeat cannot out-vote a real adapter. Enriched seeds are extended column by column into a consensus. |
| `paired_overlap` | On paired input, the insert boundary inferred from the mates themselves. Everything past it is adapter by construction, with no sequence matching involved. |

A read usually carries evidence for several of these at once, so **support is
measured once**, by a final verification pass over the same sample. Summing the
sources would count the same read three times. The sources a candidate came from
are reported as provenance.

### Confident, mixed, inconclusive

Each candidate is classified `high`, `medium` or `low` against explicit gates —
sample size, supporting reads, support fraction, length, complexity, matcher
error rate, and positional or consensus evidence — rather than an opaque score.
Every threshold is serialized into the report, so a result explains itself.

```text
confident      exactly one unrelated candidate is high-confidence
mixed          two or more are: pooled libraries, or concatenated runs
inconclusive   nothing clears the gates
```

A `mixed` result is never resolved automatically. Two spellings of the *same*
adapter are not competitors, though: the library carries frame-shifted variants
of one chemistry, so candidates are compared by sliding them past each other,
and a family reports its best-corroborated member.

**Why abundant is not the same as adapter.** A genomic repeat can be as frequent
in read tails as real read-through. What separates them is the consensus: adapter
read-through is the same sequence in every read, so it extends; a repeat is
followed by different sequence in every read, so it stops within a few bases of
its seed. Positional tests alone do not work — on a library with 40%
read-through a correctly assembled adapter scores only 27% 3'-connected and 1.0x
tail enrichment, because read-through covers the read body as densely as its
tail.

**Poly-A and poly-G are reported, not recommended.** They are the most abundant
3' k-mers in real data and are instrument or biology artifacts. They are counted
and surfaced as an artifact signal; `--poly-g` and `--poly-x` remove them.

### Output

`--format text` (default), `json` or `tsv`. JSON is the stable pipeline
interface and carries a `schema_version`; TSV writes one row per candidate for
cohort aggregation. `-o` writes to a file, atomically; otherwise the report goes
to stdout.

Adapter TSV rows contain `input, mate, decision, sequence, known_name,
known_category, confidence, evidence_sources, supporting_reads,
support_fraction, tail_connected_fraction, tail_enrichment, median_start,
median_distance_to_end, exact_matches, substitution_matches, indel_matches,
mean_error_rate, paired_overlap_support`. `known_category` is `adapter`,
`primer` or `control` for catalogued sequences and `.` for de novo candidates.

When a candidate occurs only in read tails, its body rate is zero and its
enrichment is mathematically unbounded. JSON represents this as `null`, TSV as
`tail_only`, and text as `tail only`; no format emits a non-finite number.

```text
0  the analysis completed
2  a command line, configuration or runtime error
3  the result was not confident, and --require-confident was given
```

Without `--require-confident` an inconclusive result is a *successful analysis*
and exits 0 — the command answered the question, and the answer was "not
determinable from this data".

`--emit-config PATH` writes a minimal configuration fragment, and only for a
uniquely confident result:

```toml
# Written by `bqc sniff adapters`.
[adapter]
r1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
r2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
```

which feeds the next command directly:

```bash
bqc sniff adapters sample.cbq --require-confident --emit-config sample.toml \
    --format json -o sample.adapters.json
bqc workflow sample.cbq --config sample.toml -o sample.clean.cbq
```

### Strandedness

```bash
bqc sniff strand reads.cbq --index salmon-index --format json -o strand.json
```

Requires the `sniff-strand` feature, which is off by default:

```bash
cargo install bqc --features sniff-strand
```

Strandedness is not a property of the reads on their own — it is how they relate
to *oriented transcripts* — so a reference is required and there is no
composition-based guess. A deterministic sample is mapped through Salmon's
selective-alignment mapper, using its Rust crates directly against borrowed CBQ
slices: no FASTQ conversion, no temporary files, no subprocess.

Two related results are reported. Salmon's canonical library type preserves the
full mapping orientation:

```text
ISF ISR IU   OSF OSR OU   MSF MSR MU        paired
SF  SR  U                                   single-end
```

and the pipeline classification is what most workflows actually consume:

```text
forward       forward fraction >= --stranded-threshold   (default 0.80)
reverse       reverse fraction >= --stranded-threshold
unstranded    |forward - reverse| < --unstranded-threshold (default 0.10)
undetermined  none of the above, or not enough evidence
```

with `featurecounts_strand` and `htseq_stranded` alongside them — `null` when
undetermined, because a downstream parameter must never be manufactured from an
answer nobody established.

**Evidence comes before inference.** Salmon's own inference function returns
*unstranded* from an empty count array, which is right for a quantifier and wrong
for a detector: it is indistinguishable from a genuinely unstranded library. So
the gates — `--min-informative` (5000), `--min-informative-fraction` (0.05) —
are applied first, and below them the answer is `undetermined` with
`insufficient_mapping_evidence`, not a library type.

Pair orientation is reported as measured. An outward or matching library is
surfaced with a warning rather than collapsed into inward, because that is
exactly the finding worth seeing.

**The index must be a Salmon 2.x index.** An index built by salmon 1.x uses the
pufferfish format and cannot be read; `bqc` says so and gives the rebuild
command. Its reference count, k, decoy state and content hashes are recorded with
the result, so a report can be matched back to the index that produced it.

Sampling stops at a batch boundary once `--target-informative` (50000)
observations have accumulated, so a high mapping rate does not cost a full pass —
and, because batches are fixed and counts merge by addition, the answer is
identical at every thread count.

## Splitting reads at internal adapters

`bqc segment` treats adapters as *delimiters* rather than as ends to trim
from, and emits each piece between them as its own record:

```text
source   PREFIX [ADAPTER A] SEGMENT 1 [ADAPTER B] SEGMENT 2 [ADAPTER C] SUFFIX
output   PREFIX, SEGMENT 1, SEGMENT 2, SUFFIX
```

This is concatemer resolution: linked-read, single-molecule and some amplicon
protocols produce reads containing several inserts joined by known sequences.

```bash
bqc segment concatemers.cbq -o fragments.cbq \
  --adapter-fasta delimiters.fa --segments provenance.tsv \
  --min-segment-length 20 --terminal-fragments discard
```

```text
    --adapter-r1 <SEQ>              delimiter sequence
    --adapter-fasta <PATH>          FASTA of delimiter sequences
    --terminal-fragments <MODE>     keep (default) | discard
    --min-segment-length <INT>      shortest fragment emitted [1]
    --max-segments-per-read <INT>   safety limit per source read [64]
    --segments <PATH>               provenance sidecar, TSV
```

It is a separate command because it is the one operation whose output cardinality
differs from its input's: one record becomes zero, one or many. Nothing that
assumes one-record-in-one-record-out — correction, linked adapters, 3' trimming,
paired input — composes with it. `--front`, `--quality-tail`, `--min-length` and
the rest of `trim` and `filter` do, and apply to each fragment individually.

**Single-end only.** Splitting a read into fragments has no defined effect on its
mate, so a paired input is rejected with that message rather than guessed at.

**Candidate selection.** Every coordinate at which any delimiter matches becomes
a candidate. Candidates are ordered by start coordinate, then fewer errors, then
lower error rate, then greater overlap, then declaration order; they are accepted
greedily from the left, and a candidate overlapping an already accepted delimiter
is suppressed as a description of the same adapter copy. The report counts
suppressions, so an ambiguous delimiter set is visible rather than silent.

**Fragments.** A fragment is *internal* when a delimiter bounds it on both sides
and *terminal* when a read end does. Empty fragments — from a delimiter at a read
end, or two adjacent delimiters — are always discarded. `--terminal-fragments
discard` keeps only internal fragments, which also drops reads containing no
delimiter at all, since their single fragment is bounded by two read ends.
Fragments below `--min-segment-length` are discarded, and `--max-segments-per-read`
caps a pathological read; both are counted separately in the report.

**Provenance.** Each fragment's header gains a suffix, and `--segments` writes the
same facts as columns:

```text
@SRR1.42 1:N:0|segment=1|span=44-74
```

```text
source_record_index  segment_index  source_mate  start  end  length  left_adapter  right_adapter  original_header  status  filter_reasons
1                    0              R1           0      30   30      .             left           read1            PASS    PASS
1                    1              R1           44     74   30      left          .              read1            PASS    PASS
```

`start` and `end` are coordinates in the source read, after any trimming, so a
fragment is always traceable to the bases it came from. A header-free input stays
header-free — `bqc` never invents headers — so `--segments` is *required*
there: the sidecar is then the only surviving provenance. Every emitted fragment
gets a row, accepted or rejected. Flags describe the source molecule and are
copied to every fragment of it.

Output order is source record index ascending, then segment index ascending, at
any thread count.

## Base correction

Where R1 and R2 genuinely overlap they sequence the same molecule twice, from
opposite ends, so a disagreement between them is a sequencing error in one of the
two reads. `bqc correct` (or `--correction` in a workflow) replaces the
doubtful base with the confident one:

```text
R1 >= --donor-quality and R2 <= --recipient-quality  →  correct R2 from R1
R2 >= --donor-quality and R1 <= --recipient-quality  →  correct R1 from R2
otherwise                                            →  leave both unchanged
```

```text
--donor-quality <PHRED>                  default 30
--recipient-quality <PHRED>              default 14
--correction-log <PATH>
--correction-log-detail <reads|bases>    default reads
```

Thresholds are inclusive, and `donor > recipient` is required — which is what
makes the two rules mutually exclusive, so no base is ever both donor and
recipient. Details:

* The donor base is written in the recipient's own orientation, so it is
  complemented on the way across, and the donor's **raw quality byte** is copied
  exactly. Bases that already agree are never touched and their qualities are
  never raised.
* A low-quality `N` can be corrected, but an `N` — or any non-canonical base — is
  never used as donor evidence. Those mismatches are reported separately from ones
  the quality rule declined.
* Only ungapped overlaps are corrected. `bqc`'s overlap inference is ungapped
  by construction, so every accepted overlap qualifies and the "gapped overlaps
  skipped" counter is structurally zero.
* Reads are never merged, and lengths never change.

Correction requires **paired input with stored qualities**; single-end or
quality-free input is refused before any processing starts. It reuses the same
overlap analysis and the same `--paired-overlap-min-overlap` / `--max-error-rate`
thresholds as `--paired-overlap`, and when both are enabled the alignment is
computed **once** per pair.

### Zero-copy

Every edit is planned against the original pair before anything is mutated, and
only a corrected mate is copied into reused worker buffers. A pair with no
corrections copies nothing, and if only R2 changes then R1 stays borrowed from the
memory-mapped block. `--failed-mode original` therefore still writes the genuinely
uncorrected pair.

### Correction log

`--correction-log` writes a TSV in record order. `--correction-log-detail reads`
gives one row per corrected pair:

```text
record_index  r1_header  r2_header  overlap_offset  overlap_length
overlap_mismatches  corrected_r1_bases  corrected_r2_bases
unresolved_mismatches  final_disposition
```

`bases` gives one row per corrected base, ordered by record, then mate, then
position:

```text
record_index  mate  read_position  donor_mate  donor_position
original_base  corrected_base  original_phred  corrected_phred
overlap_offset  overlap_length  final_disposition
```

Records are identified by their original CBQ index; headers are never modified and
CBQ flags are never repurposed. Tabs, newlines, carriage returns and backslashes
in headers are escaped, so a hostile header cannot spill into another row or
column.

Correction counts describe corrections **applied to the read**, before trimming: a
corrected base that a later stage trims away is still counted, and
`corrected_pairs_by_disposition` records where those pairs ended up. The
per-pair histogram and the substitution matrix appear in the JSON report only —
the console summary keeps the aggregate lines.

## Trimming

```text
--front / --front-r1 / --front-r2      fixed 5' cut
--tail  / --tail-r1  / --tail-r2       fixed 3' cut

--quality-front <PHRED> [--quality-front-window <INT>]   default window 4
--quality-tail  <PHRED> [--quality-tail-window  <INT>]
--quality-right <PHRED> [--quality-right-window <INT>]

--trim-terminal-n
--poly-g [--poly-g-min-length <INT>] [--poly-g-max-mismatch-rate <FLOAT>]
--poly-x [--poly-x-min-length <INT>] [--poly-x-max-mismatch-rate <FLOAT>]

--max-length / --max-length-r1 / --max-length-r2
```

* Quality thresholds are compared as integer sums: a window qualifies when
  `sum(phred) >= threshold * window_length`. Qualities are Phred+33; bytes
  outside `33..=126` are rejected as an encoding error.
* `--quality-front` and `--quality-tail` slide the window one base at a time
  until a qualifying window is found. When fewer bases remain than the window
  width, the remaining bases are used as a shorter window.
* `--quality-right` scans left to right and truncates the read at the **start of
  the first failing window**, which can cut earlier than the visible quality
  drop. It is mutually exclusive with `--quality-tail`.
* `--trim-terminal-n` removes contiguous `N` runs from both current boundaries
  and never touches internal `N`s.
* Poly tails remove the **longest** suffix of length `i` for which
  `i >= min_length` and `mismatches <= floor(i * max_mismatch_rate)`. `--poly-x`
  evaluates every canonical base and prefers the longest tail, breaking ties in
  `A, C, G, T` order. When both are enabled, poly-G runs first and poly-X then
  examines what remains, so no suffix is counted twice.
* `--max-length` **truncates** (a transformation). To *reject* long reads use
  `--length-limit` in `filter`.

## Filtering

```text
--min-length <INT>
--length-limit <INT>              reject reads longer than INT
--max-n <INT> / --max-n-fraction <FLOAT>
--qualified-quality <PHRED>       default 15
--max-unqualified-bases <INT> / --max-unqualified-fraction <FLOAT>
--min-mean-quality <PHRED>
--min-complexity <FLOAT>
```

Every applicable predicate is evaluated, so a rejected read reports **all** of
its reasons: `TOO_SHORT`, `TOO_LONG`, `TOO_MANY_N`, `TOO_MANY_LOW_QUAL`,
`LOW_MEAN_QUAL`, `LOW_COMPLEXITY`.

Complexity is the fraction of adjacent positions whose bases differ:

```text
complexity = |{ i : seq[i] != seq[i-1] }| / (len - 1)
```

This is fastp's lightweight metric. It detects homopolymer-like reads and is
**not** an entropy measure — it is deliberately not called entropy. Reads
shorter than two bases have no adjacent pairs and are defined to have complexity
`0.0`.

Naming note: `filter` uses `--length-limit` where the design document proposed
`--max-length`, because in `workflow` that name is already taken by the trim-time
truncation and the two mean opposite things. The TOML keys keep the documented
names (`[trim] max_length` and `[filter] max_length`) since the tables
disambiguate them.

## Paired-end behaviour

CBQ pairing is a file-wide property, so pair retention has explicit semantics.
The default is strict retention:

```text
R1 passes and R2 passes → accepted output
R1 fails or  R2 fails   → failed output (both mates together)
```

The reason sidecar records each mate independently, so a rejected pair can show
`FAIL` for one mate and `PASS` for the other.

`--pair-policy orphan --orphan-prefix <PREFIX>` keeps the surviving mate of a
broken pair instead of discarding it:

```text
clean.cbq          pairs where both mates pass        (paired schema)
PREFIX.R1.cbq      R1 passed, R2 failed               (single-end schema)
PREFIX.R2.cbq      R2 passed, R1 failed               (single-end schema)
rejected.cbq       pairs where neither mate passes    (paired schema)
```

Orphan files are single-end; the rest of the schema — qualities, headers, flags —
is preserved, and the surviving mate keeps its own header and the record's flag.
Every input record still lands in exactly one destination, and the reason sidecar
explains both mates of every non-accepted record. The orphan policy requires a
paired input and a filter stage, since nothing can be orphaned without one.

There is deliberately no "keep the pair if either mate passes" policy: it would
leave a known-failing mate in a paired file, which is exactly the ambiguity these
two policies avoid.

## Outputs

```text
-o, --output <PATH>          accepted records (required)
    --failed <PATH>          rejected records, as CBQ
    --failed-mode <MODE>     original (default) | processed
    --failed-reasons <PATH>  per-mate reason sidecar, TSV
    --pair-policy <POLICY>   strict (default) | orphan
    --orphan-prefix <PREFIX> writes <PREFIX>.R1.cbq and <PREFIX>.R2.cbq
    --report <PATH>          structured report
    --report-format <FMT>    json (default) | tsv
    --force                  allow overwriting existing files
```

`--failed-mode original` writes the untransformed record, so rejected data stay
recoverable and can be retried with different thresholds; the sidecar still
records what the processed state would have been. `--failed-mode processed`
writes the transformed record instead.

The sidecar identifies records by **original CBQ record index** and never
modifies read headers:

```text
record_index  mate  status  reasons               original_length  adapter_trimmed_length  final_length  adapter_name  adapter_start
10291         R1    FAIL    TOO_SHORT/TOO_MANY_N  151              40                      23            illumina      40
10291         R2    PASS    PASS                  151              151                     147           .             .
```

After filtering a header-free file, the generated indices in the output are
contiguous again and no longer correspond to input positions — that is exactly
why the sidecar uses original indices.

The JSON report contains the tool and `binseq` versions, input metadata and
schema, the **fully resolved** configuration (including defaults the user never
typed and any auto-detected adapters), stage order, thread count, record and base
counts, per-adapter hits, bases removed per operation, detection evidence, filter
reason counts and combinations, accepted/rejected/orphan counts, output paths and
throughput. A linked or segmenting run adds its own section: `linked` counts
reads by which flanks matched, and `segment` counts source records and fragments
separately — source records seen, records split, delimiters accepted, candidates
suppressed, fragments emitted by kind, fragments discarded by cause, and the
distribution of fragments per source record. Aggregate distributions belong to `bqtools qc` and are deliberately
absent.

## Common options

```text
-T, --threads <INT>          0 uses every available core
    --span <START..END>      restrict to original record indices
    --compression-level <INT>  default: inherited from the input
    --block-size <SIZE>        default: inherited from the input (accepts K/M/G)
-q, --quiet                  suppress the stderr summary
```

`--span` is interpreted in original record indices, so `--span 10..20` processes
exactly the eleventh through twentieth records of the input and writes ten
records.

## Configuration files

`bqc workflow --config <PATH>` reads TOML. Command line arguments override
file values field by field, and unknown keys are rejected.

```toml
threads = 8
steps = ["adapter", "trim", "filter"]   # optional; --no-* refines the default

[adapter]
r1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
r2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
min_overlap = 8
max_error_rate = 0.10

[trim.quality_tail]
minimum_phred = 20
window = 4

[trim]
terminal_n = true

[filter]
min_length = 30
max_n = 5
qualified_quality = 15
max_unqualified_fraction = 0.40

[output]
failed = "rejected.cbq"
report = "clean.bqc.json"
```

Presence of a subtable enables the operation: `[trim.quality_tail]` turns on
quality-tail trimming, `[trim.poly_g]` turns on poly-G trimming.

The optional features have keys too:

```toml
[adapter]
allow_indels = true
paired_overlap = true
paired_overlap_min_overlap = 30
auto_detect = true
detect_sample_size = 262144
detect_min_support = 0.01

[correction]
enabled = true
donor_quality = 30
recipient_quality = 14
log = "corrections.tsv"
log_detail = "reads"

[output]
pair_policy = "orphan"
orphan_prefix = "surviving"
```

`--config` belongs to `workflow`; the single-stage commands and `segment` take
their options on the command line. A `[segment]` table in a workflow
configuration is rejected rather than ignored, because `workflow` cannot segment:
segmentation is its own command.

## Safe defaults

Nothing biological happens unless it is requested:

* `adapter` requires an adapter source: explicit sequences, a FASTA file,
  `--paired-overlap` or `--auto-detect`.
* `trim` requires at least one trimming operation.
* `segment` requires a delimiter source: `--adapter-r1` or `--adapter-fasta`.
* `filter` requires at least one predicate.
* `workflow` errors when no stage is effectively configured.
* `--failed` and `--failed-reasons` require a filter stage; nothing can be
  rejected without one. `--pair-policy orphan` additionally requires a paired
  input and an `--orphan-prefix`.
* Quality-based operations on a CBQ file without a quality column are a hard
  error, never silently treated as high quality.
* Auto-detection refuses to trim on weak evidence rather than guessing, and
  never invents an adapter shorter than `--min-overlap`.
* Base correction requires paired input with stored qualities, and refuses a
  donor threshold that is not above the recipient threshold.
* Options that only qualify another option are rejected on their own:
  `--detect-sample-size` without `--auto-detect`, `--poly-g-min-length` without
  `--poly-g`, `--quality-tail-window` without `--quality-tail`,
  `--paired-overlap-min-overlap` without `--paired-overlap`, `--failed-mode`
  without `--failed`, and mate-2 options on a single-end file.

Window sizes, mismatch rates and overlaps do have defaults, but only once their
operation is enabled.

## Using it as a library

The core algorithms operate on byte slices and never touch CBQ, `clap` or the
filesystem, so they can be unit-tested and reused directly:

```rust
use bqc::adapter::{Adapter, AdapterParams, AdapterStage};
use bqc::process::Workflow;
use bqc::read::ReadView;

let stage = AdapterStage::new(
    vec![Adapter::new("illumina", b"AGATCGGAAGAGC")?],
    Vec::new(),
    AdapterParams::default(),
    None,
)?;
let workflow = Workflow::new(Some(stage), None, None, None, None)?;

let sequence = b"ACGTACGTACGTACGTAGATCGGAAGAGC";
let read = ReadView::new(sequence, None, 0, "R1")?;
let result = workflow.process(0, read, None)?;
assert_eq!(result.r1.final_length(), 16);
# Ok::<(), bqc::error::Error>(())
```

A transformation is a pair of coordinates (`read::Span`) into the borrowed
record: adapter removal, fixed and quality trimming, terminal-N, poly-G/poly-X
and truncation all update that span, and only the writer ever materializes the
result.

## Performance

500,000 single-end 150 bp records (41 MB CBQ) through the fused workflow with
adapter removal, quality-tail trimming, terminal-N, poly-G and four filters, on
an 8-core Ryzen 7 5700X:

| Threads | Wall time | Reads/s | Bases/s | Peak RSS |
| ------- | --------- | ------- | ------- | -------- |
| 1       | 1.96 s    | 255 k   | 38 M    | 51 MB    |
| 2       | 0.99 s    | 505 k   | 76 M    | 57 MB    |
| 4       | 0.53 s    | 943 k   | 141 M   | 66 MB    |
| 8       | 0.31 s    | 1.61 M  | 242 M   | 93 MB    |

Peak RSS includes the memory-mapped input. Output is byte-identical at every
thread count.

Two figures are worth knowing when tuning: decoding this file costs 0.12 s and
re-encoding the output costs about 1.2 s, so **compression is over half of a
typical run** and `--compression-level` is the biggest single knob available. A
full account of where the time goes, what was optimised and what was measured and
rejected is in [PROFILE.md](PROFILE.md).

`--allow-indels` runs a banded alignment per candidate coordinate and costs
roughly 20× substitution matching; `--paired-overlap` scores every placement of
R2 against R1. Both are opt-in and scale with threads. `--auto-detect` adds a
bounded sampling pass of about 2 s regardless of file size.

Base correction is nearly free once the overlap is known: on 200 000 genuinely
overlapping pairs, overlap inference alone costs 4.90 s at `-T 1`, correction
alone 4.96 s, and **both together 4.91 s** — the alignment is computed once.
Logging adds under 2%. Correction scales with threads like everything else
(4.98 s at `-T 1`, 0.49 s at `-T 16`), and costs nothing when disabled.

The optional stages cost more, on the same dataset and machine:

| Stage                              | T=1     | T=8     | Notes |
| ---------------------------------- | ------- | ------- | ----- |
| adapter, substitutions             | 1.9 s   | 0.5 s   | baseline |
| adapter, `--allow-indels`          | 89 s    | 13.5 s  | a banded DP per candidate coordinate |
| `--paired-overlap`                 | 12.5 s  | 1.8 s   | every placement of R2 against R1 |
| `--auto-detect` / `sniff adapters`  | +1.6 s  | +0.6 s  | 262 k records sampled, three passes |

`--allow-indels` is the expensive one: it aligns the adapter at every candidate
start instead of comparing bases, so it costs roughly 20× substitution matching.
Both stay deterministic and scale with threads; a seed index would cut the indel
cost substantially and is the obvious next optimization if profiling on real data
justifies it. Detection samples on one thread before the main pass, so its cost
is a fixed addition rather than a per-record one.

## Not implemented

Deliberately out of scope:

* FASTQ input/output, interleaved FASTQ, quality-encoding conversion.
* Read merging, UMI handling, deduplication, output splitting.
* Correction of gapped overlaps, and correction as a consensus (bases are only
  ever copied from one mate to the other, never merged).
* HTML reports, per-position quality tables, GC distributions.
* Hand-written SIMD, and a seed index for indel-aware matching.
* Five-prime and "anywhere" adapter modes; only 3' trimming exists.

## Known upstream limitation

A CBQ file with zero records — which `bqc` will produce if every record is
filtered out — cannot be reopened with `binseq`'s `MmapReader`
(`IndexCastingError`: the empty index cannot be cast). `binseq`'s streaming
reader and `bqc` itself handle such files correctly. This affects any
zero-record CBQ file, including ones written by `binseq`'s own writer.

## Development

```bash
cargo fmt --check
cargo clippy --all-targets
cargo test

# Deterministic benchmark fixtures, and the decode-only floor of a file.
cargo run --release --example make_fixtures -- /tmp/bqc-fixtures
cargo run --release --example decode_floor -- /tmp/bqc-fixtures/se150.cbq 8
```

`.github/workflows/ci.yml` runs those checks plus release-mode integration tests,
an `aarch64` build, a thread-determinism job that compares real multi-block output
across `-T 1/3/8`, and Miri over the pure algorithms. Miri cannot run the
integration tests: they memory-map files and call zstd through FFI.

The test suite covers the full CBQ schema matrix (16 combinations of pairing,
qualities, headers and flags), order and thread-count equivalence, standalone
versus fused equivalence, staged versus fused equivalence, the reason sidecar,
report rendering, orphan routing, indel-aware matching, paired-overlap inference,
auto-detection including its refusal path, base correction (both directions,
threshold boundaries, ambiguous donors and recipients, negative overlap offsets,
log escaping and ordering), and property tests for the span, accounting, metadata
and correction invariants.
