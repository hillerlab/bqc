#!/usr/bin/env python3
# Copyright (c) 2026 The Hiller Lab at the Senckenberg Gesellschaft für Naturforschung
# Distributed under the terms of the GNU General Public License, Version 3.0.

"""Streams output FASTQs once and emits one TSV line of run metrics.

Usage:
  metrics.py --tag NAME --r1 F1 --r2 F2 [--u1 F] [--u2 F] [--kmers A[:20],B[:20]]

Reads are 4-line FASTQ. Pair accounting uses the header id (suffixes /1, /2,
whitespace, and trimmomatic's *_val_1 style mangling are stripped). Adapter
residue counts records containing any of the given 20-mers (sense strand;
adapters are never reverse-complemented by trimming tools in the reads).

`--summarize DIR` mode folds the per-run files into one median-per-tool table.
"""

import argparse
import glob
import json
import os
import re
import statistics


def parse_fastq(path):
    """Yields (id, seq, qual) tuples."""
    if path is None or not os.path.exists(path):
        return
    with open(path) as f:
        while True:
            name = f.readline()
            if not name:
                return
            seq = f.readline().strip()
            f.readline()  # +
            qual = f.readline().strip()
            yield name[1:].split()[0], seq, qual


def strip_mate(id_):
    id_ = re.sub(r"/[12]$", "", id_)
    return re.sub(r"_(val|trimmed)_[12]$", "", id_)


def measure(r1, r2, u1, u2, kmers):
    ids1, ids2 = set(), set()
    bases1 = bases2 = q30_1 = q30_2 = n1 = n2 = 0
    hits1 = hits2 = 0
    len1 = len2 = 0

    def scan(path, ids, bases, q30, n, hits, lens, r2mate):
        nonlocal bases1, bases2, q30_1, q30_2, n1, n2, hits1, hits2, len1, len2
        for id_, seq, qual in parse_fastq(path):
            ids.add(strip_mate(id_))
            if r2mate:
                bases2 += len(seq); q30_2 += sum(1 for q in qual if q >= "?")
                n2 += 1; len2 += len(seq)
                if any(k in seq for k in kmers):
                    hits2 += 1
            else:
                bases1 += len(seq); q30_1 += sum(1 for q in qual if q >= "?")
                n1 += 1; len1 += len(seq)
                if any(k in seq for k in kmers):
                    hits1 += 1

    scan(r1, ids1, bases1, q30_1, n1, hits1, len1, False)
    scan(r2, ids2, bases2, q30_2, n2, hits2, len2, True)
    scan(u1, ids1, 0, 0, 0, 0, 0, False)  # unpaired: only the id set matters
    scan(u2, ids2, 0, 0, 0, 0, 0, True)

    both = len(ids1 & ids2)
    only1 = len(ids1 - ids2)
    only2 = len(ids2 - ids1)
    q30frac = (q30_1 + q30_2) / (bases1 + bases2) if bases1 + bases2 else 0.0
    ml1 = len1 / n1 if n1 else 0.0
    ml2 = len2 / n2 if n2 else 0.0
    ppm = (hits1 + hits2) / (n1 + n2) * 1e6 if n1 + n2 else 0.0
    return {
        "reads1": n1, "reads2": n2, "bases1": bases1, "bases2": bases2,
        "meanlen1": ml1, "meanlen2": ml2, "pairs_both": both,
        "r1_only": only1, "r2_only": only2, "q30frac": q30frac,
        "adapter_ppm": ppm,
    }


def run_one(args):
    kmers = [k.strip() for k in (args.kmers or "").split(",") if k.strip()]
    m = measure(args.r1, args.r2, args.u1, args.u2, kmers)
    fields = [
        args.tag, m["reads1"], m["reads2"], m["bases1"], m["bases2"],
        f"{m['meanlen1']:.1f}", f"{m['meanlen2']:.1f}", m["pairs_both"],
        m["r1_only"], m["r2_only"], f"{m['q30frac']:.4f}", f"{m['adapter_ppm']:.1f}",
    ]
    print("\t".join(map(str, fields)))


def summarize(dir_):
    """Folds per-run NAME.TOOL.json (hyperfine) + .metrics files into one
    median-per-tool table; RSS column from NAME.TOOL.rss when measured."""
    print("scenario\ttool\twall_median_s\trss_MB\treads1\treads2\tbases1\tbases2\tmeanlen1\tmeanlen2\tpairs_both\tr1_only\tr2_only\tq30frac\tadapter_ppm")
    for json_path in sorted(glob.glob(os.path.join(dir_, "*.json"))):
        base = os.path.basename(json_path)[:-5]  # NAME.TOOL
        scenario, tool = base.rsplit(".", 1)
        mp = json_path[:-5] + ".metrics"
        if not os.path.exists(mp):
            continue
        metrics_line = open(mp).readline().rstrip("\n")
        if metrics_line.endswith("\tFAILED"):
            print(f"{scenario}\t{tool}\tFAILED")
            continue
        times = json.load(open(json_path))["results"][0]["times"]
        wall = f"{statistics.median(times):.2f}"
        rss_path = json_path[:-5] + ".rss"
        rss = open(rss_path).read().strip() if os.path.exists(rss_path) else "-"
        print(f"{scenario}\t{tool}\t{wall}\t{rss}\t{metrics_line.split(chr(9),1)[1]}")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="mode", required=True)
    one = sub.add_parser("one")
    one.add_argument("--tag", required=True)
    one.add_argument("--r1", required=True)
    one.add_argument("--r2", required=True)
    one.add_argument("--u1")
    one.add_argument("--u2")
    one.add_argument("--kmers", default="")
    s = sub.add_parser("summarize")
    s.add_argument("dir")
    args = ap.parse_args()
    if args.mode == "one":
        run_one(args)
    else:
        summarize(args.dir)


if __name__ == "__main__":
    main()
