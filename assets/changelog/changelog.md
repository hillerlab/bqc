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
