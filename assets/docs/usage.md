<p align="center">

<p align="center">
  <picture>
    <source
      media="(prefers-color-scheme: dark)"
      srcset="../figures/bqc-dark.png"
    >
    <source
      media="(prefers-color-scheme: light)"
      srcset="../figures/bqc-light.png"
    >
    <img
      width="200"
      alt="Hiller Lab"
      src="../figures/bqc-light.png"
    >
  </picture>
</p>

<p align="center">
  <picture>
    <source
      media="(prefers-color-scheme: dark)"
      srcset="../figures/hillerlab-dark.png"
    >
    <source
      media="(prefers-color-scheme: light)"
      srcset="../figures/hillerlab-light.png"
    >
    <img
      width="200"
      alt="Hiller Lab"
      src="../figures/hillerlab-light.png"
    >
  </picture>
</p>

  <span>
    <h1 align="center">
        bqc
    </h1>
  </span>

  <span>
    <h2 align="center">
        USER GUIDE
    </h2>
  </span>

  <p align="center">
    <a href="https://github.com/hillerlab/bqc" reference="_blank">
      <img alt="GitHub License" src="https://img.shields.io/github/license/hillerlab/bqc?color=blue">
    </a>
  </p>

  <p align="center">
    <samp>
        <span> CBQ-native all-in-one quality control tool </span>
        <br>
        <span> The Hiller Lab at the Senckenberg Research Institute </span>
        <br>
        <br>
        <a href="https://github.com/ArcInstitute/binseq">binseq</a> .
        <a href="https://github.com/hillerlab/bqc/blob/master./docs/usage.md">usage</a> .
        <a href="https://hillerlab.com/">us</a> 
    </samp>
  </p>

</p>

---

`bqc` removes adapter sequences, trims and filters sequencing reads stored
in [CBQ](https://docs.rs/binseq) files. It works natively on CBQ: no FASTQ
conversion, no intermediate file, and every stage runs in one pass over the
input.

```bash
bqc adapter   reads.cbq -o out.cbq --adapter-r1 AGATCGG...
bqc adapter   reads.cbq -o out.cbq --auto-detect          # infer the adapter
bqc trim      reads.cbq -o out.cbq --quality-tail 20 --trim-terminal-n
bqc filter    reads.cbq -o out.cbq --min-length 30 --failed rejected.cbq
bqc correct   pairs.cbq -o out.cbq --correction-log corrections.tsv
bqc workflow  reads.cbq -o out.cbq --config illumina.toml -T 8
bqc umi       reads.cbq -o out.cbq --umi-location read1 --umi-length 8
bqc dedup     reads.cbq -o unique.cbq -T 8
bqc sniff adapters reads.cbq                              # inspect, never modify
bqc sniff strand   reads.cbq --index salmon-index         # RNA-seq orientation
```

Use `bqtools` for conversion to/from FASTQ, concatenation, splitting and
dataset-level QC reports; `bqc` itself only reads and writes CBQ.

## 1. Installation

```bash
cargo install bqc
```

`bqc sniff strand` is optional and off by default — its mapping engine
pulls in a large dependency tree and needs a C compiler:

```bash
cargo install bqc --features sniff-strand
```

Everything in this guide was verified against bqc 0.0.1; run `bqc
--help` (or `<command> --help`) for the live option list.

## 2. Quick start

Clean a read file end to end in one pass:

```bash
bqc workflow reads.cbq -o clean.cbq \
  --adapter-r1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
  --quality-tail 20 --trim-terminal-n --poly-g \
  --min-length 30 --max-n 5 -T 8
```

Don't know the adapter? Ask the file first (`bqc sniff adapters
reads.cbq`, or `bqc sniff strand reads.cbq --index salmon-index` for
RNA-seq orientation).

## 3. How it works

### One pass, nothing in between

Every command decodes the input once and writes CBQ back out. The stages of a
run always execute in this fixed order:

```text
umi → correct → adapter → trim → filter
```

`workflow` runs any combination of these in that single pass. The order is
part of the output contract: UMI removal runs first so a read-head UMI never
leaks into correction or trimming, then correction rescues bases from the
mate so they are visible to every later decision.

### Guarantees

* **Deterministic.** The same input and configuration produce byte-for-byte
  identical output at any thread count and on any machine.
* **Schema-preserving.** Pairing, qualities, headers and flags are read from
  the input file header and reproduced exactly; a header-free file stays
  header-free.
* **Atomic.** Every output is written to a hidden temporary path and renamed
  only after the run succeeds; a failed run leaves no partial output.
* **Bounded memory.** The input is memory-mapped; peak memory is roughly the
  mapped file plus a few MB per thread, never the whole dataset decoded.

### Pairing

CBQ pairing is a file-wide property. The default is strict retention: a pair
is accepted only if **both** mates pass, and both mates go to the failed
output if either one fails. There is deliberately no "keep the pair if one
mate passes" policy; use `--pair-policy orphan` to rescue surviving mates
instead (see [Outputs](#7-outputs)).

## 4. Options shared by every command

```text
-o, --output <PATH>             accepted output CBQ file (required)
-T, --threads <INT>             worker threads; 0 uses every available core
    --span <START..END>         restrict processing to original record indices
    --compression-level <INT>   output zstd level (default: inherited from input)
    --block-size <SIZE>         output block size, e.g. 1M (default: inherited)
    --report <PATH>             write a structured report
    --report-format <FMT>       json (default) | tsv
    --failed <PATH>             rejected records, as CBQ (requires filtering)
    --failed-reasons <PATH>     per-mate reason sidecar, TSV
    --failed-mode <MODE>        original (default) | processed
    --pair-policy <POLICY>      strict (default) | orphan
    --orphan-prefix <PREFIX>    writes <PREFIX>.R1.cbq and <PREFIX>.R2.cbq
    --force                     overwrite existing output files
-q, --quiet                     suppress the stderr summary
```

`--span` uses 0-based original record indices: `--span 10..20` processes
exactly the eleventh through twentieth records; open ends are allowed (`10..`,
`..50`, `..`).

## 5. Commands

### 5.1 `bqc adapter`

Removes adapter sequences from the 3' end of each read.

```text
    --adapter-r1 <SEQ>          adapter trimmed from R1
    --adapter-r2 <SEQ>          adapter trimmed from R2
    --adapter-fasta <PATH>      FASTA of adapters applied to both mates
    --min-overlap <INT>         minimum adapter/read overlap [8]
    --max-error-rate <FLOAT>    max mismatch fraction over the overlap [0.10]
    --max-errors <INT>          absolute cap on mismatches (in addition)
    --allow-indels              count insertions/deletions while matching
    --paired-overlap            infer insert length from the R1/R2 overlap
    --paired-overlap-min-overlap <INT>   min accepted overlap [30]
    --auto-detect               infer adapters from the data
    --detect-sample-size <INT>  sample for auto-detection [262144]
    --detect-min-support <FLOAT> min support fraction [0.01]
```

**The matching rule.** For every candidate coordinate `p` in the read, the
adapter prefix is compared against the read from `p` to the end. A match
requires `overlap >= --min-overlap` (default 8), `errors <= floor(overlap *
--max-error-rate)`, and, when supplied, `errors <= --max-errors`. Among
matches the winner is the **earliest coordinate**, then fewest errors, then
longest overlap — trimming at the earliest coordinate removes the longest
contaminated suffix.

Notes:

* `N` never matches — in the read or the adapter — and counts as one error,
  deliberately conservative: a poly-N tail would otherwise match any adapter.
  A partial adapter only matches at the 3' end; an adapter copy in the middle
  of a read is not trimmed.
* `--adapter-fasta` adapters apply to both mates; `--adapter-r1`/`--adapter-r2`
  are mate-specific.

**`--allow-indels`** replaces the substitution scan with a banded
edit-distance alignment; insertions and deletions each cost one error from
the same budget. It is roughly 20x slower than substitution matching, so
reserve it for chemistries that really produce adapter indels.

**`--paired-overlap`** needs no adapter sequence at all (paired input only).
R1 and R2 are read from opposite ends of one insert, so where they genuinely
overlap, the reverse complement of R2 aligns to R1. A successful alignment
reveals the insert length; everything past it in either mate is adapter
read-through, and sequence matching is skipped for that pair.

**`--auto-detect`** infers the adapter from the data before trimming, using
the same detector that `bqc sniff adapters` reports (there is exactly
one). The mates are inferred independently: a mate with no evidence gets no
adapter. Detection **refuses to trim rather than guess** when nothing clears
the evidence gates (the run aborts and asks for explicit sequences) or when
two unrelated adapters clear them — a mixed library is a fact about the data,
not a tie to break; run `bqc sniff adapters` to see both candidates and
choose.

#### Linked adapters

A linked adapter is a *pair* of flanks with the sequence of interest between
them — how amplicon and many small-RNA libraries are built. Instead of
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
    --linked-max-3p-overhang <INT>  3' flank bases that may hang past the end
    --linked-min-insert-length <INT> shortest insert retained [1]
    --linked-unmatched <POLICY>     continue (default) | keep | fail
```

The 5' flank is anchored near the read start; the 3' flank is found by the
ordinary 3' matcher and must begin after the 5' flank ends.
`--linked-unmatched` decides the rest of the library: `continue` (default)
runs ordinary adapter matching on unmatched reads, `keep` writes them
unchanged, `fail` rejects them with the reason `LINKED_UNMATCHED`.

### 5.2 `bqc trim`

Shortens reads by position, quality, terminal Ns or homopolymer tails.

```bash
bqc trim reads.cbq -o trimmed.cbq --quality-tail 20 --trim-terminal-n --poly-g
```

```text
--front / --front-r1 / --front-r2     fixed 5' cut
--tail  / --tail-r1  / --tail-r2      fixed 3' cut
--quality-front <PHRED> [--quality-front-window <INT>]    default window 4
--quality-tail  <PHRED> [--quality-tail-window  <INT>]    default window 4
--quality-right <PHRED> [--quality-right-window <INT>]    default window 4
--trim-terminal-n
--poly-g [--poly-g-min-length <INT>] [--poly-g-max-mismatch-rate <FLOAT>]
--poly-x [--poly-x-min-length <INT>] [--poly-x-max-mismatch-rate <FLOAT>]
--max-length / --max-length-r1 / --max-length-r2
```

* **Quality windows** slide one base at a time and keep the first window whose
  Phred sum is at or above `threshold * window_length`. `--quality-tail` cuts
  from the 3' end, `--quality-front` from the 5' end, `--quality-right`
  truncates at the start of the first failing window (mutually exclusive with
  `--quality-tail`). Quality operations on a quality-free CBQ are a hard
  error.
* **`--trim-terminal-n`** removes contiguous `N` runs from both ends; internal
  Ns are never touched.
* **Poly tails** remove the longest qualifying suffix of length >= the minimum
  (defaults: min length 10, mismatch rate 0.10). `--poly-x` considers every
  base and prefers the longest tail; poly-G runs first when both are enabled.
* **`--max-length` truncates** (a transformation). To *reject* long reads, use
  `--length-limit` in `filter`.

### 5.3 `bqc filter`

Accepts or rejects reads (and pairs) against per-read predicates.

```text
--min-length <INT>                    reject reads shorter than INT
--length-limit <INT>                  reject reads longer than INT
--max-n <INT> / --max-n-fraction <FLOAT>       reject too many Ns
--qualified-quality <PHRED>           Phred at or above which a base counts [15]
--max-unqualified-bases <INT> / --max-unqualified-fraction <FLOAT>
--min-mean-quality <PHRED>            reject reads below the mean quality
--min-complexity <FLOAT>              reject low-complexity reads
```

Every applicable predicate is evaluated, so a rejected read reports **all** of
its reasons: `TOO_SHORT`, `TOO_LONG`, `TOO_MANY_N`, `TOO_MANY_LOW_QUAL`,
`LOW_MEAN_QUAL`, `LOW_COMPLEXITY`. Complexity is the fraction of adjacent
positions whose bases differ — fastp's lightweight metric, not entropy
(`complexity = |{ i : seq[i] != seq[i-1] }| / (len - 1)`).

Naming note: `filter` uses `--length-limit` where you might expect
`--max-length` — in `workflow` that name is already taken by trim-time
truncation, and the two mean opposite things.

### 5.4 `bqc correct`

Corrects low-quality bases from the other mate where the pair overlaps.
Requires **paired input with stored qualities**; single-end or quality-free
input is refused before processing starts.

```bash
bqc correct pairs.cbq -o corrected.cbq --correction-log corrections.tsv
```

```text
--donor-quality <PHRED>              lowest Phred accepted as donor [30]
--recipient-quality <PHRED>          highest Phred that may be overwritten [14]
--correction-log <PATH>              write a correction log, TSV
--correction-log-detail <DETAIL>     reads (default) | bases
--paired-overlap-min-overlap <INT>   minimum accepted overlap [30]
--max-error-rate <FLOAT>             max mismatch fraction [0.10]
```

Where R1 and R2 overlap they sequence the same molecule from opposite ends, so
a disagreement between them is a sequencing error in one mate. Correction
applies per base: `R1 >= donor and R2 <= recipient` writes R1's base into R2,
`R2 >= donor and R1 <= recipient` writes R2's base into R1, and any other
combination leaves both unchanged.

The donor base is written in the recipient's own orientation, and the donor's
raw quality byte is copied exactly. Reads are never merged, lengths never
change, and bases that already agree are never touched. The overlap alignment
is computed once per pair even when `--paired-overlap` is also enabled.

### 5.5 `bqc segment`

Splits reads at internal adapter occurrences, treating adapters as
*delimiters* rather than ends to trim from:
`PREFIX [A] SEGMENT 1 [B] SEGMENT 2 [C] SUFFIX` becomes
`PREFIX, SEGMENT 1, SEGMENT 2, SUFFIX`. One read becomes zero, one or many
records — which is why it is its own command and cannot be part of a
`workflow`. **Single-end only**; a paired input is rejected.

```text
    --adapter-r1 <SEQ>              delimiter sequence
    --adapter-fasta <PATH>          FASTA of delimiter sequences
    --terminal-fragments <MODE>     keep (default) | discard
    --min-segment-length <INT>      shortest fragment emitted [1]
    --max-segments-per-read <INT>   safety limit per source read [64]
    --segments <PATH>               provenance sidecar, TSV
```

A fragment is *internal* when a delimiter bounds it on both sides, *terminal*
when a read end does. Empty fragments are always discarded;
`--terminal-fragments discard` keeps only internal fragments, which also
drops reads containing no delimiter at all. `--segments` writes one row per
emitted fragment — on a header-free input it is **required**, the sidecar
being the only surviving provenance. `trim` and `filter` options apply to
each fragment individually, so the `--front`/`--quality-tail`/`--min-length`
family is accepted here too.

### 5.6 `bqc workflow`

Runs any combination of the stages in one pass.

```bash
bqc workflow reads.cbq -o clean.cbq --config illumina.toml -T 8
```

```text
--steps <LIST>       stages to run, e.g. umi,correct,adapter,trim,filter
--no-adapter         skip the adapter stage
--no-trim            skip the trim stage
--no-filter          skip the filter stage
--correction         enable paired-overlap base correction
--config <PATH>      TOML configuration file (see below)
```

`--steps` and the `--no-*` flags conflict with each other. Every stage not
explicitly skipped is requested, but a stage only runs if it is configured —
so a workflow with just `--min-length 30` runs only the filter stage, and a
workflow with no stage effectively configured is an error.

### 5.7 `bqc umi`

Extracts a unique molecular identifier and relocates it into the read name,
matching fastp's `--umi`. This is extraction/relocation only — **not** UMI
family clustering, error correction or consensus calling.

```bash
bqc umi reads.cbq -o umi.cbq --umi-location read1 --umi-length 8
```

```text
--umi-location <LOCATION>    read1|read2|index1|index2|per_index|per_read
--umi-length <INT>           UMI length, required for read-derived UMIs
--umi-skip <INT>             bases to skip after the UMI [0]
--umi-prefix <STR>           prefix before the UMI in the name [empty]
--umi-delimiter <STR>        delimiter before the UMI [":"]
```

| location | UMI source | sequence modification |
|---|---|---|
| `index1` | first index | none |
| `index2` | second index | none |
| `read1` | R1 prefix | remove UMI + skip from R1 |
| `read2` | R2 prefix | remove UMI + skip from R2 |
| `per_index` | `index1_index2` | none |
| `per_read` | `r1umi_r2umi` | remove prefixes from both |

The tag is inserted before the first space in the header (appended when there
is no space) and applied to both mate names. A read shorter than
`length + skip` is an explicit error — `bqc` never silently truncates a UMI,
unlike fastp. UMI processing requires stored read headers; `read2`, `index2`,
`per_index` and `per_read` require paired input.

UMI removal runs **before** correction (it is the first stage of a workflow),
so a read-head UMI never participates in overlap inference, adapter matching,
trimming or filtering, and is never credited to adapter removal in the report.

### 5.8 `bqc dedup`

Removes exact duplicate reads across the whole dataset, keeping the earliest
occurrence unchanged. It is a separate command, not a workflow stage, because
deduplication needs global cross-block state that would dismantle the
workflow's single-pass, independent-block property.

```bash
bqc dedup reads.cbq -o unique.cbq -T 8
```

```text
--memory-mb <INT>   memory budget in MiB for the Bloom filters and the exact
                    candidate arena [1024]
```

Two records are duplicates when their sequence payloads are byte-for-byte
equal; qualities, names and flags do not participate. For paired input the
ordered `(R1, R2)` pair is the key, so `(AC, CGT)` never aliases `(ACC, GT)`.

The pipeline is two passes: a Bloom discovery pass marks fingerprints that
repeat, then an exact classifier re-reads in input order and verifies each
candidate against the exact bytes, so hash or Bloom collisions can only add
work — never drop a record. The result is **exact** (zero false deletions)
and deterministic. If the candidate arena exceeds the memory budget the run
aborts rather than silently switching to approximate deletion.

Deduplication operates on the sequences physically present in the input —
sequence deduplication, not UMI-aware molecular deduplication (that needs
alignment coordinates and belongs elsewhere).

### 5.9 `bqc sniff adapters`

Infers which adapter sequences contaminate the reads. **Non-destructive**: it
never trims, filters, reorders or rewrites anything — the input is opened
read-only and is byte-identical afterwards.

```bash
bqc sniff adapters reads.cbq                          # human-readable summary
bqc sniff adapters reads.cbq --format json -o adapters.json
```

```text
--sample-size <INT>      records sampled, spread evenly [262144]
--span <START..END>      restrict sampling to original record indices
-T, --threads <INT>      worker threads; 0 uses every available core
-o, --output <PATH>      write the report here instead of stdout
--format <FORMAT>        text (default) | json | tsv
--require-confident      exit with status 3 when not confident
--force                  overwrite an existing output file
--top <INT>              candidates reported per mate [5]
--emit-config <PATH>     write a TOML fragment for uniquely confident results
```

**Sampling is distributed, not "the first N reads".** CBQ is indexed and
block-addressable, so the sample is spread evenly across the file — a
leading-prefix sample would be biased by concatenated lanes and joined runs.
There is no random generator and no seed, so every run gives the same answer,
at every thread count; with fewer records than the sample size, all are
inspected, with a warning.

**Evidence comes from three sources:**

| Source | What it is |
| --- | --- |
| `known_database` | A versioned library of 234 published adapter, primer and PhiX sequences, matched with the same matcher that trims. |
| `kmer_consensus` | Exact 10-mers of each read's final 60 bases, counted once per read, extended column by column into a consensus. |
| `paired_overlap` | On paired input, the insert boundary inferred from the mates themselves; everything past it is adapter by construction. |

A read usually carries evidence from several of these at once, so **support
is measured once** by a final verification pass over the same sample; the
sources are provenance only. Each candidate is classified `high`, `medium` or
`low` against explicit gates (sample size, supporting reads, support
fraction, length, complexity, matcher error rate, positional evidence) —
never an opaque score. The mate-level (and file-level) decision:

```text
confident      exactly one unrelated candidate is high-confidence
mixed          two or more are: pooled libraries, or concatenated runs
inconclusive   nothing clears the gates
```

A `mixed` result is never resolved automatically. Two spellings of the *same*
adapter are not competitors: candidates are slid past each other, so
frame-shifted variants of one chemistry count as one family, and a much
weaker unrelated hit is demoted to `medium` rather than making every file
`mixed`. Poly-A and poly-G are the most abundant 3' k-mers in real data and
are instrument or biology artifacts; they are surfaced as an "Artifact
signal" in the text report, to be removed with `--poly-g` / `--poly-x` —
reported, never recommended.

Example output (real run, paired 150 bp reads):

```text
Input: reads.cbq
Sampling: deterministic-distributed

R1 adapter result: confident
  Sequence: AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
  Known as: illumina-truseq
  Confidence: high
  Supporting reads: 47,056 (18.82%)
```

**`--emit-config`** writes a minimal valid `bqc` TOML fragment, and only
for a uniquely confident result — never for `mixed` or `inconclusive` —
feeding the next command directly:

```bash
bqc sniff adapters sample.cbq --require-confident --emit-config sample.toml \
    --format json -o sample.adapters.json
bqc workflow sample.cbq --config sample.toml -o sample.clean.cbq
```

The emitted fragment looks like this:

```toml
# Written by `bqc sniff adapters`.
[adapter]
r1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
r2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
```

`--auto-detect` on the `adapter` command uses the same detector for its
recommendation instead of a report, so the two can never disagree.

### 5.10 `bqc sniff strand`

Infers RNA-seq library strandedness by mapping reads against a Salmon
transcriptome index. Strandedness is not a property of the reads on their
own — it is how they relate to oriented transcripts — so **a reference is
required** and there is no composition-based guess. Requires the `sniff-strand`
feature and a Salmon 2.x index (`salmon 1.x` pufferfish indexes are refused,
with the rebuild command in the message).

```bash
bqc sniff strand reads.cbq --index salmon-index --format json -o strand.json
bqc sniff strand reads.cbq --transcriptome transcripts.fa   # index built on the fly
```

```text
--index <PATH>                   Salmon 2.x transcriptome index directory
--transcriptome <PATH>           transcriptome FASTA; index built on the fly
--sample-size <INT>              records examined at most [1000000]
--target-informative <INT>       observations after which sampling may stop [50000]
--min-informative <INT>          observations below which no answer is given [5000]
--min-informative-fraction <FLOAT>  required informative share [0.05]
--stranded-threshold <FLOAT>     forward/reverse share at which a library is
                                 called stranded [0.80]
--unstranded-threshold <FLOAT>   forward/reverse difference below which a
                                 library is called unstranded [0.10]
```

`--index` and `--transcriptome` are mutually exclusive; exactly one is
required. **`--transcriptome`** takes a transcriptome FASTA instead of a ready
index and builds a Salmon index from it on the fly in a temporary directory,
discarded when the run ends. There is no persistent index cache — the index is
rebuilt on every run — so reach for `--index` when the same reference is used
repeatedly. The file must exist and be a FASTA; a missing path is a
configuration error.

Two related results are reported:

**Salmon's canonical library type** preserves the full mapping orientation:

```text
ISF ISR IU   OSF OSR OU   MSF MSR MU        paired
SF  SR  U                                   single-end
```

**The pipeline strandedness** — what most workflows actually consume:

```text
forward       forward fraction >= --stranded-threshold
reverse       reverse fraction >= --stranded-threshold
unstranded    |forward - reverse| < --unstranded-threshold
undetermined  none of the above, or not enough evidence
```

with `featurecounts_strand` and `htseq_stranded` alongside:

| Result         | featureCounts `-s` | HTSeq `--stranded` |
| -------------- | -----------------: | ------------------ |
| `unstranded`   | `0`                | `no`               |
| `forward`      | `1`                | `yes`              |
| `reverse`      | `2`                | `reverse`          |
| `undetermined` | `null`             | `null`             |

`null` for `undetermined` is deliberate: a downstream parameter must never be
manufactured from an answer nobody established.

**Evidence comes before inference.** Below the gates — `--min-informative`
(5000) observations and a 5% informative fraction — the answer is
`undetermined` with the reason `insufficient_mapping_evidence`, never a
library type. (Salmon's own inference returns *unstranded* from an empty
count array, which is right for a quantifier and wrong for a detector.)

For paired input, pair orientation is reported as measured — `inward`,
`outward`, `matching` or `undetermined` — with a warning when it is not
`inward`.

Example output (real run, first-strand paired library):

```text
Library type: ISR
Strandedness: reverse
Pair orientation: inward
Informative fragments: 20,000 (100.0%)
Forward fraction: 0.0%
Reverse fraction: 100.0%
Decision: confident
```

The report also records index provenance — reference count, k-mer length,
decoy state and content hashes — so a result can be matched back to the
index that produced it.

## 6. Sniff reports

Both sniff commands accept `--format text` (default), `json` or `tsv`, and
`-o` writes the report atomically (stdout otherwise).

* **JSON** is the stable pipeline interface. It carries a `schema_version`
  and a common envelope: `tool`, `input` (path, records, schema flags),
  `sample` (method, range, requested/selected counts), `parameters` (every
  threshold used, including the classification gates), `result` and
  `warnings`.
* **TSV** is for cohort aggregation: `sniff adapters` writes one row per
  candidate — `input, mate, decision, sequence, known_name, known_category,
  confidence, evidence_sources, supporting_reads, support_fraction,
  tail_connected_fraction, tail_enrichment, median_start,
  median_distance_to_end, exact_matches, substitution_matches,
  indel_matches, mean_error_rate, paired_overlap_support` — and `sniff
  strand` one summary row — `input, decision, salmon_library_type,
  strandedness, pair_orientation, forward_fraction, reverse_fraction,
  informative_records, informative_fraction, records_examined,
  featurecounts_strand, htseq_stranded`.

For adapter candidates supported only in the tail, the body rate is zero and
the enrichment ratio is mathematically unbounded. JSON represents that state as
`null`, TSV as `tail_only`, and text as `tail only`; no projection emits a
non-finite number. `known_category` is `adapter`, `primer` or `control` for a
catalogued sequence and `null`/`.` for a de novo candidate.

## 7. Outputs

Besides the accepted CBQ file, commands can write:

* **`--failed` / `--failed-mode`.** Rejected records, as CBQ. `--failed-mode
  original` (default) writes the untransformed record so rejected data stay
  recoverable; `processed` writes the transformed record instead.
* **`--failed-reasons`.** One row per mate, keyed by **original CBQ record
  index** (which is why the sidecar is needed after filtering — output indices
  are contiguous again and no longer correspond to input positions):

  ```text
  record_index  mate  status  reasons               original_length  adapter_trimmed_length  final_length  adapter_name  adapter_start
  10291         R1    FAIL    TOO_SHORT/TOO_MANY_N  151              40                      23            illumina      40
  10291         R2    PASS    PASS                  151              151                     147           .             .
  ```

* **`--pair-policy orphan`** (with `--orphan-prefix`) keeps the surviving mate
  of a broken pair instead of discarding it: `clean.cbq` holds pairs where
  both mates pass (paired schema), `PREFIX.R1.cbq` / `PREFIX.R2.cbq` the
  single surviving mates (single-end schema), and `rejected.cbq` pairs where
  neither mate passes. Every input record lands in exactly one destination.
  The orphan policy requires a paired input and a filter stage.

* **`--report`.** A structured JSON (default) or TSV report with the tool and
  `binseq` versions, input metadata and schema, the **fully resolved**
  configuration (including defaults and any auto-detected adapters), stage
  order, thread count, record and base counts, per-adapter hits, bases removed
  per operation, detection evidence, filter reason counts, accepted/rejected/
  orphan counts, output paths and throughput. A `correct` run adds correction
  counts, a linked run flank-match counts, a segment run source/fragment
  counts.
* **`--correction-log`** (from `correct` or a workflow) writes one row per
  corrected pair (`record_index, r1_header, r2_header, overlap_offset,
  overlap_length, overlap_mismatches, corrected_r1_bases, corrected_r2_bases,
  unresolved_mismatches, final_disposition`), or one per corrected base with
  `--correction-log-detail bases`.
* **`--segments`** (from `segment`) writes one row per fragment
  (`source_record_index, segment_index, source_mate, start, end, length,
  left_adapter, right_adapter, original_header, status, filter_reasons`),
  where `start`/`end` are coordinates in the source read after any trimming.

## 8. Exit codes

```text
0  the command completed
2  a command line, configuration or runtime error
3  the result was not confident, and --require-confident was given
```

Exit code 3 applies to `sniff adapters` and `sniff strand` only. Without
`--require-confident`, an inconclusive result is a *successful analysis* and
exits 0 — the command answered the question, and the answer was "not
determinable from this data".

## 9. Configuration files

`bqc workflow --config <PATH>` reads TOML. Command line arguments override
file values field by field, and unknown keys are rejected. Single-stage
commands take their options on the command line only.

```toml
threads = 8
steps = ["umi", "adapter", "trim", "filter"]   # optional; --no-* refines the default

[adapter]
r1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
r2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
min_overlap = 8
max_error_rate = 0.10
allow_indels = true                     # optional features
paired_overlap = true
paired_overlap_min_overlap = 30
auto_detect = true
detect_sample_size = 262144
detect_min_support = 0.01

[trim.quality_tail]                     # presence enables the operation
minimum_phred = 20
window = 4

[trim]
terminal_n = true

[filter]
min_length = 30
max_n = 5
qualified_quality = 15
max_unqualified_fraction = 0.40

[umi]
location = "read1"
length = 8
skip = 0
prefix = "UMI"
delimiter = ":"

[correction]
enabled = true
donor_quality = 30
recipient_quality = 14
log = "corrections.tsv"
log_detail = "reads"

[output]
failed = "rejected.cbq"
report = "clean.bqc.json"
pair_policy = "orphan"
orphan_prefix = "surviving"
```

Presence of a subtable enables the operation: `[trim.quality_tail]` turns on
quality-tail trimming, `[trim.poly_g]` turns on poly-G trimming.

A `[segment]` table in a workflow configuration is rejected rather than
ignored: `workflow` cannot segment, because segmentation changes output
cardinality.

## 10. A complete pipeline

```bash
# 1. What adapters are in this file? Write a report and a usable config.
bqc sniff adapters sample.cbq --require-confident \
    --emit-config sample.adapters.toml --format json -o sample.adapters.json
# 2. Is the library stranded? (needs a Salmon index)
bqc sniff strand sample.cbq --index reference.salmon \
    --require-confident --format json -o sample.strand.json
# 3. Clean the data in one pass using the detected adapters.
bqc workflow sample.cbq --config sample.adapters.toml \
    -o sample.clean.cbq -T 8
```

Exit codes make each step safe to gate: the pipeline stops with code 3 when a
sniff result is not confident, instead of proceeding on a guess.

## 11. Gotchas and safe defaults

* **Nothing biological happens unless requested.** `adapter` requires an
  adapter source (sequences, FASTA, `--paired-overlap` or `--auto-detect`);
  `trim` requires at least one trimming operation; `segment` requires a
  delimiter source; `filter` requires at least one predicate; a `workflow`
  with nothing configured is an error.
* **Quality operations on a quality-free CBQ are a hard error**, never
  silently treated as high quality.
* Options that only qualify another option are rejected on their own:
  `--detect-sample-size` without `--auto-detect`, `--poly-g-min-length`
  without `--poly-g`, `--paired-overlap-min-overlap` without
  `--paired-overlap`, `--failed-mode` without `--failed`, and mate-2 options
  (`--front-r2`, `--tail-r2`, `--max-length-r2`) on a single-end file.
* **Do not sum the sniff evidence sources.** Support is measured once by the
  verification pass; a read may appear in several sources, and the sources
  are provenance, not tallies.
* **Abundant does not mean adapter.** A genomic repeat can be as frequent in
  read tails as real read-through; what separates them is the consensus
  extension, which is why `sniff adapters` can report `inconclusive` on a
  repeat-rich file. Poly-A/poly-G tails are an artifact signal, never
  recommended.
* **Zero-record output caveat.** A CBQ file with zero records (e.g. everything
  filtered out) cannot be reopened with `binseq`'s `MmapReader`; `bqc`
  itself reads and writes such files correctly.
* **Command-specific constraints.** `segment` is single-end only, `correct`
  is paired-with-qualities only, and `--orphan-prefix` needs a paired input,
  a filter stage and `--pair-policy orphan`. Each violation is refused with a
  clear message before any processing starts.
* **An inconclusive sniff is information, not failure.** Exit code 0 with a
  report beats guessing; use `--require-confident` when the pipeline must
  branch.
