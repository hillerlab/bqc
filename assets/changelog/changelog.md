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
        CHANGELOG
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

All notable changes to this project are documented in this file.

## 0.0.3 — 2026-08-13

### Added

- **`umi`** — UMI extraction and relocation, matching fastp's `--umi`: the
  six locations `read1`/`read2`/`index1`/`index2`/`per_index`/`per_read`, a
  fastp-compatible read-name tag (`--umi-prefix`/`--umi-delimiter`, inserted
  before the first space of both mate names), and `--umi-skip`. Read-derived
  UMIs physically remove the prefix; a read shorter than `length + skip` is
  an explicit error rather than a silently truncated UMI. Also a `[umi]` stage
  of `workflow`, running **before** correction so a read-head UMI never takes
  part in overlap inference, adapter matching, trimming or filtering.
- **`dedup`** — exact whole-dataset deduplication. Two passes: a Bloom
  discovery pass marks fingerprints that repeat, then an exact classifier
  re-reads in input order and byte-verifies every candidate, so hash or Bloom
  collisions can only add work, never drop a record. Parallel decode/hash and
  encode stages run around a serial ordered classifier, keeping the result
  exact, deterministic and first-occurrence preserving; memory scales with
  candidate families rather than unique reads (`--memory-mb`).
- Benchmarks for both commands (vs fastp `--umi` and `--dedup`) and the
  matching user-guide and README documentation.
  
## 0.0.2 — 2026-08-09

Auto-detection is no longer a gate on the run. Detection is advisory: a file
where nothing clears the confidence gates passes through untrimmed, and a
mixed library is trimmed with the strongest match instead of aborting.

### Changed

- `--auto-detect` no longer aborts on an inconclusive or mixed detection:
  a mate with no confident candidate contributes nothing (the file passes
  through untrimmed), and a mate with several unrelated candidates above the
  gates is trimmed with the strongest one.
- The adapter report names the confidence behind a non-confident picked
  sequence (`high`/`medium`/`low`), so a mixed choice is never presented as
  settled.
- An empty `AdapterStage` is the documented pass-through, confirming that a
  run which detects nothing still writes its output.
- `sniff adapters --emit-config` still requires a uniquely confident result;
  the "writes a config only when it can" contract is unchanged.
- Benchmark scripts moved under `bench/scripts` and the README highlights
  updated.

## 0.0.1 — 2026-08-06

Initial release of `bqc`, a CBQ-native all-in-one QC tool: it decodes a CBQ
file once, runs its stages in a single pass, and writes CBQ back out — no
FASTQ conversion, no intermediate files.

### Added

- **`adapter`** — 3' adapter removal with substitution or (`--allow-indels`)
  banded edit-distance matching, `--paired-overlap` insert-boundary inference,
  `--auto-detect` (shared with `sniff adapters`), and linked (amplicon/small-RNA)
  5'/3' flanks via `--linked-*`.
- **`trim`** — quality-window trimming (`--quality-front/tail/right`), fixed
  `--front`/`--tail` cuts, `--trim-terminal-n`, `--poly-g`/`--poly-x`,
  `--max-length`.
- **`filter`** — per-read predicates: length, N count/fraction, unqualified
  bases, mean quality, complexity.
- **`correct`** — paired-overlap base correction using the other mate's
  high-quality base, with a `--correction-log`.
- **`segment`** — splits reads at internal adapter occurrences; single-end only.
- **`workflow`** — any combination of correct → adapter → trim → filter in one
  pass, driven by TOML config or `--steps`.
- **`sniff adapters`** (reference-free) and **`sniff strand`** (RNA-seq
  orientation against a Salmon index, or a `--transcriptome` FASTA with the
  index built on the fly; optional `sniff-strand` feature).
- Deterministic, schema-preserving, atomic (temp-file + rename),
  bounded-memory, multi-threaded (`-T`) processing.
- Structured `--report` (JSON/TSV), `--failed` sidecars, `--span` record
  ranges, exit code 3 with `--require-confident`.

### Fixed

- Fix `map(...).unwrap_or(...)` clippy lints in the report input metadata
  (`src/report.rs`) and strand fold-back naming (`src/sniff/strand.rs`) by using
  `Result::map_or`.
- Collapse the canonical-base byte array into a byte-string literal
  (`src/trim.rs`) to satisfy the `byte-char-slices` clippy lint.
- Run the slow `determinism` and `miri` CI jobs only on pull requests and
  releases, not on every push; the rest of the suite still runs on push.
