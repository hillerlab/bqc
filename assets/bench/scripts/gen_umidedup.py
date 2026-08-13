#!/usr/bin/env python3
# Copyright (c) 2026 The Hiller Lab at the Senckenberg Gesellschaft für Naturforschung
# Distributed under the terms of the GNU General Public License, Version 3.0.

"""Deterministic single-end benchmark data for the UMI and dedup features.

Pure-ACGT 150 bp reads, fixed seed, all-Q40 qualities. `umi` carries an 8 bp
read1 UMI prefix; `dedup_<rate>` tiles a unique pool to the requested
duplication rate (0/1/10/50/90%). Stdlib only, fixed seed: re-running with the
same arguments gives the same bytes.
"""

import argparse
import os
import random

BASE = b"ACGT"


def seq(rng, n):
    return "".join(chr(BASE[b & 3]) for b in rng.randbytes(n))


def write_fastq(path, reads, read_len):
    qual = "I" * read_len
    chunk = []
    with open(path, "w", buffering=1 << 20) as fh:
        for i, s in enumerate(reads):
            chunk.append(f"@read_{i}\n{s}\n+\n{qual}\n")
            if len(chunk) >= 100_000:
                fh.write("".join(chunk))
                chunk.clear()
        if chunk:
            fh.write("".join(chunk))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("outdir")
    ap.add_argument("-n", type=int, default=2_000_000)
    ap.add_argument("--read-len", type=int, default=150)
    ap.add_argument("--umi-len", type=int, default=8)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    rng = random.Random(20_260_813)
    n = args.n
    rl = args.read_len

    reads = [seq(rng, args.umi_len) + seq(rng, rl - args.umi_len) for _ in range(n)]
    write_fastq(f"{args.outdir}/umi.fq", reads, rl)
    del reads

    for rate in (0, 1, 10, 50, 90):
        unique = max(1, round(n * (100 - rate) / 100))
        pool = [seq(rng, rl) for _ in range(unique)]
        reads = [pool[i % unique] for i in range(n)]
        write_fastq(f"{args.outdir}/dedup_{rate}.fq", reads, rl)
        del pool, reads


if __name__ == "__main__":
    main()
