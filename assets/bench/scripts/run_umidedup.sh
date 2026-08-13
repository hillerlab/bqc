#!/usr/bin/env bash
# Copyright (c) 2026 The Hiller Lab at the Senckenberg Gesellschaft für Naturforschung
# Distributed under the terms of the GNU General Public License, Version 3.0.

# bench/run_umidedup.sh — UMI and dedup benchmark: bqc (CBQ) vs fastp (FASTQ).
# Hyperfine-timed (3 runs), RSS from a single /usr/bin/time -v pass. Reproduces
# the "UMI and deduplication" section of BENCHMARK.md.
#
#   bench/run_umidedup.sh            everything (2M reads, 6 scenarios)
#
# Data lives in bench/data, results in bench/results. Generation is skipped
# when the data already exists.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
BENCH="$(dirname "$HERE")"
BQC="$BENCH/../../bqc/target/release/bqc"
CONDA=/home/alejandro/.local/share/mamba/condabin/conda
DATA="$BENCH/data"
RESULTS="$BENCH/results"
THREADS=8
RUNS=3
N=2000000

mkdir -p "$DATA" "$RESULTS"

if [ ! -f "$DATA/umi.cbq" ]; then
    echo "[data] generating" >&2
    python3 "$HERE/gen_umidedup.py" "$DATA" -n "$N"
    echo "[data] encoding to CBQ" >&2
    for f in umi dedup_0 dedup_1 dedup_10 dedup_50 dedup_90; do
        bqtools encode -f q -m cbq "$DATA/$f.fq" -o "$DATA/$f.cbq" 2>/dev/null
    done
fi

# bench NAME CMD... -> prints "NAME  wall_s  rss_mb"
bench() {
    local name="$1"; shift
    local cmd="$*" json="$RESULTS/$name.json" wd="$RESULTS/$name.work"
    rm -rf "$wd"; mkdir -p "$wd"
    hyperfine --runs "$RUNS" --style none --export-json "$json" "$cmd" >/dev/null 2>&1
    local wall rss
    wall=$(jq -r '.results[0].median' "$json")
    /usr/bin/time -v -o "$RESULTS/$name.rss.log" "$@" >/dev/null 2>&1 || true
    rss=$(awk '/Maximum resident/{printf "%.0f", $6/1024}' "$RESULTS/$name.rss.log")
    printf '%s\t%s\t%s\n' "$name" "$wall" "$rss"
}

FASTP_OFF="--disable_adapter_trimming --disable_trim_poly_g --disable_quality_filtering --disable_length_filtering"

echo -e "scenario\twall_s\trss_mb"

bench umi.bqc    "$BQC" umi "$DATA/umi.cbq" -o "$RESULTS/umi.bqc.out.cbq" -T "$THREADS" --umi-location read1 --umi-length 8 --force
bench umi.fastp  "$CONDA" run -n fastp --no-capture-output fastp -i "$DATA/umi.fq" -o "$RESULTS/umi.fastp.out.fq" --umi --umi_loc read1 --umi_len 8 $FASTP_OFF --dont_eval_duplication -w "$THREADS" -j /dev/null -h /dev/null

for rate in 0 1 10 50 90; do
    bench "dedup_${rate}.bqc"   "$BQC" dedup "$DATA/dedup_${rate}.cbq" -o "$RESULTS/dedup_${rate}.bqc.out.cbq" -T "$THREADS" --force
    bench "dedup_${rate}.fastp" "$CONDA" run -n fastp --no-capture-output fastp -i "$DATA/dedup_${rate}.fq" -o "$RESULTS/dedup_${rate}.fastp.out.fq" --dedup $FASTP_OFF -w "$THREADS" -j /dev/null -h /dev/null
done
