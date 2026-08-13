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
        BENCHMARK
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

# bqc vs. the ecosystem: end-to-end benchmark

Compares `bqc` against the four tools most often used for the same job —
`fastp`, `cutadapt`, `atropos`, `trimmomatic` — on identical paired-end data,
with one matched ("parity") configuration per tool. Everything here is
reproducible with `bench/run.sh`.

## Questions

1. How do the tools compare on adapter-removal fidelity (residue left behind)?
2. How do they react to edge cases: adapter dimers, indel-damaged adapters,
   short adapter stubs, poly-A/G tails, low-quality tails, all-N reads?
3. What is the throughput and memory trade-off, bqc's CBQ I/O included?

## Environment

| | |
|---|---|
| CPU | AMD Ryzen 7 5700X, 8 cores / 16 threads |
| RAM | 125 GB |
| OS | Linux |
| bqc | this tree, `cargo build --release` (edition 2024, rust 1.91+) |
| bqtools | 0.5.12 (`cargo install bqtools`; CBQ support needs ≥ 0.5) |
| fastp | 1.0.1 (env `fastp`) |
| cutadapt | 5.2 (env `cutadapt`) |
| atropos | 1.1.32 (env `cutadapt`) |
| trimmomatic | 0.41 (env `trimmomatic`) |

`trim_galore` 2.3.0 was installed and smoke-tested, then **discarded**: it is
single-threaded, accepts only its bundled adapter list (no explicit
`-a/-A`), and at ~1M pairs/min it alone would have doubled the benchmark wall
time. Its behavior is a thin wrapper over cutadapt, which is benchmarked
directly.

All tools are invoked through their conda environments except bqc/bqtools.

## Feature matrix

| feature | bqc | fastp | cutadapt | atropos | trimmomatic |
|---|---|---|---|---|---|
| explicit 3' adapter (paired) | `--adapter-r1/r2` | `--adapter_sequence(_r2)` | `-a/-A` | `-a/-A` | ILLUMINACLIP fasta |
| adapter auto-detection | `--auto-detect` (library + de novo) | built-in + poly-A/G heuristics | no | no | no |
| indel-tolerant adapter match | opt-in `--allow-indels` | limited | default (error-rate incl. indels) | default | no |
| paired-end mode | yes | yes | yes | yes | yes |
| 3' quality trim | `--quality-tail 20` | `--cut_right 20` | `-q 20` | `-q 20` | SLIDINGWINDOW:4:20 |
| min-length filter | `--min-length` | `--length_required` | `-m` | `--minimum-length` | MINLEN |
| poly-A/G trim | `--poly-a/--poly-g` | `--trim_poly_g` | adapter seqs | adapter seqs | no |
| multi-threaded | `-T 8` | `-w 8` | `-j 8` | `--threads 8` | `-threads 8` |
| machine-readable report | JSON | JSON + HTML | stderr | JSON | log |
| native input | CBQ (mmap) | gz FASTQ | gz FASTQ | gz FASTQ | gz FASTQ |

## Data

### Synthetic

`bench/gen.py` (stdlib-only, xorshift64*, fixed seed): identical reads are
served to every tool — FASTQ for the ecosystem, `bqtools encode -m cbq -S 4`
→ CBQ for bqc. Read model is Illumina 3' read-through: short inserts
continue into the adapter and flow-cell poly-A, the shape real data has.

| scenario | size | contents |
|---|---|---|
| `typical` | 2M pairs | 150 bp; 30% read-through (insert 10–149); 10% Q2 tails |
| `dimer` | 200k pairs | insert 0–19 bp: reads are nearly all adapter |
| `indel` | 200k pairs | read-through adapter carries 1 insertion + 1 substitution |
| `partial` | 200k pairs | reads end inside the adapter (6–15 bp stubs) |
| `polya` | 200k pairs | 50% of R1 end in poly-A, R2 in poly-G |
| `lowqual` | 200k pairs | last 40–80 bases at Q2 |
| `shortreads` | 200k pairs | 60–80 bp reads, 30% read-through |
| `alln` | 200k pairs | all-N reads, Q38 |

Sizes are set so a per-tool run lands in the ~20 s–2 min band: large enough
for real runtime differences (2M pairs ≈ 0.6 Gbp; CBQ mmap vs gzip I/O is
the point), small enough that the whole benchmark finishes in about an hour.
Correctness probes (the edge scenarios) are 10× smaller because their answer
is binary, not a runtime.

### Real

SRR8997011 (paired, Illumina; 31M pairs) subsampled to the first 2M pairs —
same per-run band as `typical`. FASTQ tools read the plain-FASTQ subsample;
the same reads are encoded to CBQ once (`bench/data/real.cbq`), encode cost
excluded from timing.

## Protocol

- Timing: `hyperfine --runs 3 --style none` per (scenario, tool) call; the
  table reports the median of the 3 runs. No warmup run — every tool sees
  the same (unwarmed) cache state. Peak RSS is measured separately with a
  single `/usr/bin/time -v` pass for `typical` and `real` only (10 runs
  total; hyperfine does not measure memory).
- Parity configuration, one per tool (adapters are the TruSeq 5' 34-mers
  `AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC` / `AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT`;
  quality trim Q20, min length 20, 8 threads):

  ```bash
  # bqc (CBQ in/out; decode for metric collection)
  bqc workflow in.cbq -o out.cbq -T 8 \
      --adapter-r1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
      --adapter-r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
      --quality-tail 20 --min-length 20

  fastp -i R1 -I R2 -o O1 -O O2 \
      --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
      --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
      --cut_right --cut_mean_quality 20 --length_required 20 -w 8

  cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
      -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -q 20 -m 20 -j 8 \
      -o O1 -p O2 R1 R2

  atropos trim -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
      -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -q 20 --minimum-length 20 \
      --threads 8 -o O1 -p O2 -pe1 R1 -pe2 R2

  trimmomatic PE -threads 8 -phred33 R1 R2 O1P O1U O2P O2U \
      ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:20
  ```

- Outputs are uncompressed FASTQ; bqc output decoded with `bqtools
  decode` (decode and metric passes are not timed). `conda run` startup
  (~0.5 s) is included in FASTQ-tool times and matters only at the 10k-pair
  smoke scale (`--quick`).

## Metrics

Collected by `bench/metrics.py` in one streaming pass over each output:
reads out per mate, bases out, mean read length, pairs kept (both / R1-only /
R2-only), fraction of Q≥30 bases, and adapter residue in ppm of reads
(records containing the adapter's 5' 20-mer; the indel scenario scans the
damaged adapter's 20-mer). 0-length and missing outputs are reported as zeros.

## Results

### Throughput and memory (median of 3, hyperfine)

`bench/results/summary.tsv` (reproduce: `bench/run.sh`).

| dataset | tool | wall s | RSS MB | reads out | bases out (M) | mean len | pairs kept | adapter ppm |
|---|---|---|---|---|---|---|---|---|
| typical | bqc | **1.59** | 176 | 1,957,469 | 251.7 | 128.6 | 1,957,469 | 0.0 |
| typical | fastp | 6.46 | 1235 | 1,957,469 | 251.3 | 128.4 | 1,957,469 | 0.0 |
| typical | cutadapt | 3.86 | 50 | 1,957,469 | 251.2 | 128.3 | 1,957,469 | 0.0 |
| typical | atropos | 16.65 | 60 | 1,957,469 | 251.2 | 128.3 | 1,957,469 | 0.0 |
| typical | trimmomatic | 12.95 | 1154 | 1,957,469 | 251.8 | 128.7 | 1,957,469 | 0.0 |
| real | bqc | **1.23** | 103 | 1,965,023 | 196.6 | 100.0 | 1,965,023 | 0.0 |
| real | fastp | 5.00 | 1240 | 1,828,671 | 180.9 | 98.3 | 1,828,671 | 0.0 |
| real | cutadapt | 3.08 | 51 | 1,964,259 | 196.2 | 99.7 | 1,964,259 | 0.0 |
| real | atropos | 15.43 | 62 | 1,964,259 | 196.2 | 99.7 | 1,964,259 | 0.0 |
| real | trimmomatic | 7.72 | 1149 | 1,921,288 / 1,863,229 | 188.1 / 181.4 | 97.6 | 1,827,052 | 0.0 |

On `typical`, all five tools keep the same 1,957,469 pairs at the same mean
length with zero adapter residue — the parity configs are genuinely
equivalent, so the runtime column is the comparison that matters. bqc is
the fastest (1.6 s; 2.4× cutadapt, 4× fastp, 8× trimmomatic, 10× atropos)
and sits mid-pack on memory: 176 MB vs cutadapt's 50 MB (lowest) and
fastp's 1.2 GB (highest). The real dataset shows the same ordering; fastp
keeps 7% fewer pairs (its per-base quality adaptation differs on real
data) and trimmomatic drops mates independently (94k R1-only, 36k R2-only).

### Edge-case reactions (200k pairs, 3 runs, median)

Adapter residue (ppm) and pair retention (%) per tool per scenario. 0 reads
out ⇒ `—`.

| scenario | metric | bqc | fastp | cutadapt | atropos | trimmomatic |
|---|---|---|---|---|---|---|
| dimer | reads kept | 0 | 0 | 0 | 0 | 0 |
| indel | reads kept / residue ppm | 200k / 431,300 | 200k / 431,288 | 185,662 / 0 | 185,662 / 0 | 187,094 / 41,757 |
| indel + `--allow-indels` | reads kept / residue ppm | **185,662 / 0** | | | | |
| partial | mean len (adapter stub cut?) | 140.8 | 139.5 | 139.5 | 139.5 | 150.0 (no cut) |
| polya | bases out (M) | 30.0 | 30.0 | 30.0 | 30.0 | 30.0 |
| lowqual | mean len (Q-trim cut?) | 92.5 | 90.5 | 90.4 | 90.4 | 90.5 |
| shortreads | reads kept | 200k | 200k | 200k | 200k | 200k |
| alln | reads kept | 200k | 0 | 200k | 200k | 0 |

- `indel`: the indel-tolerant matchers (cutadapt/atropos default error
  models, trimmomatic partially) find and remove the damaged adapter;
  bqc's default exact match and fastp's seed-based match pass it
  through (43% residue). With opt-in `--allow-indels`, bqc matches
  cutadapt/atropos exactly (185,662 pairs kept @85.0, 0 ppm vs their
  185,662 @84.8), confirming the option works as documented.
- `partial`: trimmomatic's ILLUMINACLIP seed-and-threshold heuristic does
  not clip ≤15 bp stubs; everyone else does (140.8 vs 139.5).
- `lowqual`: all tools trim the Q2 tails; bqc leaves ~2 bp more than the
  cutadapt family (92.5 vs 90.4–90.5).
- `alln`: fastp and trimmomatic drop all-N reads outright; bqc,
  cutadapt, and atropos keep them (documented per-tool policy, not a bug).


## UMI and deduplication

A focused benchmark for the two features without an ecosystem analog in the
parity config above: UMI extraction/relocation and exact deduplication.
Single-end synthetic data (2M reads, 150 bp, fixed seed, all-Q40 qualities)
so the comparison isolates the operation itself. Same environment as above
(Ryzen 7 5700X, 8 threads; fastp 1.0.1 via env `fastp`; bqc this tree).

### Data

`bench/scripts/gen_umidedup.py` (stdlib-only, fixed seed) writes one FASTQ per
scenario; `bqtools encode -m cbq` produces the CBQ bqc reads.

| dataset | reads | contents |
|---|---|---|
| `umi` | 2M | 8 bp random read1 UMI prefix + 142 bp insert |
| `dedup_0` … `dedup_90` | 2M | 150 bp reads tiled from a unique pool to the named duplication rate (0/1/10/50/90%) |

### Protocol

```bash
# UMI (read1, 8 bp)
bqc umi in.cbq -o out.cbq -T 8 --umi-location read1 --umi-length 8
fastp -i in.fq -o out.fq --umi --umi_loc read1 --umi_len 8 -w 8 \
    --disable_adapter_trimming --disable_trim_poly_g \
    --disable_quality_filtering --disable_length_filtering \
    --dont_eval_duplication

# Dedup
bqc dedup in.cbq -o out.cbq -T 8
fastp -i in.fq -o out.fq --dedup -w 8 \
    --disable_adapter_trimming --disable_trim_poly_g \
    --disable_quality_filtering --disable_length_filtering
```

Timing is `hyperfine --runs 3` (median) with peak RSS from a single
`/usr/bin/time -v` pass. fastp dedup runs at its default accuracy level 3 (a
fixed 4 GB Bloom); bqc at its default `--memory-mb 1024`.

### Results

#### UMI (8 bp read1, 2M reads)

| tool | wall s | RSS MB |
|---|---|---|
| bqc | **0.48** | 100 |
| fastp | 1.32 | 80 |

#### Dedup (2M reads)

| dup rate | bqc wall s | bqc RSS MB | fastp wall s | fastp RSS MB |
|---|---|---|---|---|
| 0% | 1.05 | 612 | 3.21 | 4204 |
| 1% | 1.07 | 691 | 3.32 | 4200 |
| 10% | 1.36 | 1078 | 3.43 | 4194 |
| 50% | 2.09 | 1473 | 3.61 | 4192 |
| 90% | 1.13 | 964 | 3.11 | 4181 |

Two differences follow directly from the designs:

* **Memory scales with candidate families.** bqc's two-Bloom + exact-arena
  design holds only the repeated families (612 MB at 0% duplication — the
  Bloom alone; 1.47 GB at 50% — 1M families in the arena). fastp pins 4.2 GB
  regardless of duplication. bqc's peak is a third of fastp's floor.
* **Exact vs approximate.** fastp's dedup is a Bloom filter and can delete a
  unique read on a hash collision; bqc verifies every candidate against exact
  bytes.

| rate | bqc removed | fastp removed | fastp false deletions |
|---|---|---|---|
| 10% | 200,000 | 200,168 | 168 |
| 50% | 1,000,000 | 1,000,065 | 65 |

bqc hits the exact count; fastp silently drops a few extra unique reads.

## Limitations

- bqc reads CBQ and writes CBQ; FASTQ tools read plain FASTQ (the real
  subsample is plain text). Times include each platform's I/O (CBQ mmap +
  zstd vs streaming text) but not conversion (encode/decode excluded; only
  bqc's own decode is not timed, the harness decodes post-run).
- Parity config is one point in a wide space: bqc's `--allow-indels` is
  off, so the indel scenario measures default behavior; cutadapt-family
  default error models include indels.
- Known tool quirks surfaced while building this: bqc writes a valid-but-
  empty CBQ when every read is rejected and `bqtools decode` crashes on it
  (harness falls back to zero metrics); trimmomatic refuses to guess quality
  encoding for uniformly-high synthetic quality (`-phred33` forced);
  trimmomatic's ILLUMINACLIP `2:30:10` does not clip adapter stubs shorter
  than its ~10 bp simple-clip threshold (the `partial` scenario), a
  deliberate behavior of its seed-and-threshold heuristic.
