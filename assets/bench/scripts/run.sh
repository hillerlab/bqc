#!/usr/bin/env bash
# Copyright (c) 2026 The Hiller Lab at the Senckenberg Gesellschaft für Naturforschung
# Distributed under the terms of the GNU General Public License, Version 3.0.

# bench/run.sh — one benchmark run: 8 synthetic scenarios + real SRR8997011,
# hyperfine-timed (3 runs), metrics from output FASTQs. Reproduces every
# number in BENCHMARK.md.
#
#   bench/run.sh            everything: synthetic (typical 2M / edge 200k
#                           pairs) + real (2M pairs) + summary.tsv
#   bench/run.sh --quick    smoke test: all scenarios at 10k pairs, 1 run,
#                           isolated in data/quick and results/quick
#
# Data lives in bench/data, results in bench/results. Generation is skipped
# when a scenario's .cbq exists; delete it to regenerate.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(dirname "$HERE")"
BQC="$ROOT/bqc/target/release/bqc"
CONDA=/home/alejandro/.local/share/mamba/condabin/conda
TRIMMO_ADAPTERS="$HOME/.local/share/mamba/envs/trimmomatic/share/trimmomatic/adapters/TruSeq3-PE-2.fa"
A1=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
A2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
D1=AGATCGTGAAGAGCACACGTCTGAACTCCAGTCAC  # damaged adapter (T @6), indel scenario
SCENARIOS=(typical dimer indel partial polya lowqual shortreads alln)
THREADS=8
RUNS=3
P_TYPICAL=2000000   # 2M pairs: per-run ~20 s-2 min, differences visible
P_EDGE=200000       # correctness probes; seconds suffice
P_REAL=2000000
MEMORY=1            # extra /usr/bin/time -v pass for RSS (typical + real only)

if [ "${1:-}" = "--quick" ]; then
    RUNS=1; P_TYPICAL=10000; P_EDGE=10000; P_REAL=0; MEMORY=0
    DATA="$HERE/data/quick"; RESULTS="$HERE/results/quick"
else
    DATA="$HERE/data"; RESULTS="$HERE/results"
fi

kmers_for() {  # adapter 20-mers scanned for residue; indel scans the damaged one
    case "$1" in
        typical|dimer|shortreads) echo "${A1:0:20},${A2:0:20}" ;;
        indel) echo "${D1:0:20},${A2:0:20}" ;;
        *) echo "" ;;
    esac
}

cmd_for() {  # cmd_for TOOL R1 R2 CBQ WD -> the tool's command line
    local tool="$1" r1="$2" r2="$3" cbq="$4" wd="$5"
    case "$tool" in
        bqc) echo "$BQC workflow '$cbq' -o '$wd/out.cbq' -T $THREADS --adapter-r1 '$A1' --adapter-r2 '$A2' --quality-tail 20 --min-length 20 --force" ;;
        fastp) echo "$CONDA run -n fastp --no-capture-output fastp -i '$r1' -I '$r2' -o '$wd/o1.fq' -O '$wd/o2.fq' --adapter_sequence '$A1' --adapter_sequence_r2 '$A2' --cut_right --cut_mean_quality 20 --length_required 20 -w $THREADS" ;;
        cutadapt) echo "$CONDA run -n cutadapt --no-capture-output cutadapt -a '$A1' -A '$A2' -q 20 -m 20 -j $THREADS -o '$wd/o1.fq' -p '$wd/o2.fq' '$r1' '$r2'" ;;
        atropos) echo "$CONDA run -n cutadapt --no-capture-output atropos trim -a '$A1' -A '$A2' -q 20 --minimum-length 20 --threads $THREADS -o '$wd/o1.fq' -p '$wd/o2.fq' -pe1 '$r1' -pe2 '$r2'" ;;
        trimmomatic) echo "$CONDA run -n trimmomatic --no-capture-output trimmomatic PE -threads $THREADS -phred33 '$r1' '$r2' '$wd/o1p.fq' '$wd/o1u.fq' '$wd/o2p.fq' '$wd/o2u.fq' ILLUMINACLIP:$TRIMMO_ADAPTERS:2:30:10 SLIDINGWINDOW:4:20 MINLEN:20" ;;
    esac
}

bench() {  # bench NAME TOOL R1 R2 CBQ RSS[0|1]
    local name="$1" tool="$2" r1="$3" r2="$4" cbq="$5" rss_extra="${6:-0}"
    local wd="$RESULTS/$name.$tool.work" json="$RESULTS/$name.$tool.json"
    # trimmomatic 0.41 has an intermittent FastqRecord threading
    # race (crashes mid-write, exit 1); retry — every other tool is
    # deterministic and only ever needs its one attempt.
    local attempt ok=0
    for attempt in 1 2 3; do
        rm -rf "$wd"; mkdir -p "$wd"
        echo "[$(date +%H:%M:%S)] $name $tool (attempt $attempt)"
        if hyperfine --runs "$RUNS" --style none --export-json "$json" \
                "$(cmd_for "$tool" "$r1" "$r2" "$cbq" "$wd")" 2>"$wd/stderr.log"; then
            ok=1; break
        fi
        echo "[$(date +%H:%M:%S)] $name $tool attempt $attempt failed"
    done
    if [ "$ok" = 1 ]; then
        if [ "$tool" = bqc ]; then
            # bqtools decode crashes on a valid-but-empty CBQ (all
            # reads rejected, e.g. dimer); emit zero metrics instead.
            if bqtools decode "$wd/out.cbq" --prefix "$wd/out" -f q 2>"$wd/decode.err"; then
                python3 "$HERE/metrics.py" one --tag "$name.bqc" \
                    --r1 "$wd/out_R1.fq" --r2 "$wd/out_R2.fq" --kmers "$(kmers_for "$name")" \
                    > "$RESULTS/$name.bqc.metrics"
            else
                printf '%s\t0\t0\t0\t0\t0.0\t0.0\t0\t0\t0\t0.0\t0.0\n' \
                    "$name.bqc" > "$RESULTS/$name.bqc.metrics"
            fi
        else
            local or1="$wd/o1.fq" or2="$wd/o2.fq" u1= u2=
            [ "$tool" = trimmomatic ] && or1="$wd/o1p.fq" && or2="$wd/o2p.fq" \
                && u1="--u1 $wd/o1u.fq" && u2="--u2 $wd/o2u.fq"
            python3 "$HERE/metrics.py" one --tag "$name.$tool" \
                --r1 "$or1" --r2 "$or2" $u1 $u2 --kmers "$(kmers_for "$name")" \
                > "$RESULTS/$name.$tool.metrics"
        fi
        if [ "$MEMORY" = 1 ] && [ "$rss_extra" = 1 ]; then
            /usr/bin/time -v -o "$wd/rss.log" bash -c "$(cmd_for "$tool" "$r1" "$r2" "$cbq" "$wd")" 2>/dev/null || true
            awk '/Maximum resident/{print $6/1024}' "$wd/rss.log" > "$RESULTS/$name.$tool.rss" || true
        fi
    else
        printf '%s\tFAILED\n' "$name.$tool" > "$RESULTS/$name.$tool.metrics"
    fi
}

run_scenario() {  # run_scenario NAME R1 R2 CBQ RSS[0|1]
    local name="$1" r1="$2" r2="$3" cbq="$4" rss_extra="${5:-0}"
    for tool in bqc fastp cutadapt atropos trimmomatic; do
        bench "$name" "$tool" "$r1" "$r2" "$cbq" "$rss_extra"
    done
}

make_data() {
    for s in "${SCENARIOS[@]}"; do
        if [ ! -f "$DATA/$s.cbq" ]; then
            local pairs=$P_EDGE
            [ "$s" = typical ] && pairs=$P_TYPICAL
            echo "[$(date +%H:%M:%S)] generating $s ($pairs pairs)"
            python3 "$HERE/gen.py" "$s" -n "$pairs" -o "$DATA/$s"
            bqtools encode -f q -S 4 -m cbq "$DATA/${s}_R1.fastq" "$DATA/${s}_R2.fastq" -o "$DATA/$s.cbq"
        fi
    done
}

run_real() {
    # bare `return` would propagate the failed test's status and
    # set -e would kill the script — return 0 explicitly.
    [ "$P_REAL" -gt 0 ] || return 0
    local real=/home/alejandro/Documents/projects/pipelines/big_samples/SRR8997011
    local r1="$DATA/real_R1.fq" r2="$DATA/real_R2.fq"
    # { zcat || true; } | head -N — with pipefail, head closing the
    # pipe SIGPIPEs zcat and set -e kills the run (data is complete though);
    # piping through gzip would also leave a footerless stream. Plain text
    # has no wrapper to corrupt; all tools read plain FASTQ.
    [ -f "$r1" ] || { zcat "$real"_1.fastq.gz 2>/dev/null || true; } | head -$((P_REAL * 4)) > "$r1"
    [ -f "$r2" ] || { zcat "$real"_2.fastq.gz 2>/dev/null || true; } | head -$((P_REAL * 4)) > "$r2"
    [ -f "$DATA/real.cbq" ] || bqtools encode -f q -S 4 -m cbq "$r1" "$r2" -o "$DATA/real.cbq"
    run_scenario real "$r1" "$r2" "$DATA/real.cbq" 1
}

mkdir -p "$DATA" "$RESULTS"
make_data
for s in "${SCENARIOS[@]}"; do
    rss_extra=0; [ "$s" = typical ] && rss_extra=1
    run_scenario "$s" "$DATA/${s}_R1.fastq" "$DATA/${s}_R2.fastq" "$DATA/$s.cbq" "$rss_extra"
done
run_real
python3 "$HERE/metrics.py" summarize "$RESULTS" | tee "$RESULTS/summary.tsv"
