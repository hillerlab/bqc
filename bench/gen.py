#!/usr/bin/env python3
# Copyright (c) 2026 The Hiller Lab at the Senckenberg Gesellschaft für Naturforschung
# Distributed under the terms of the GNU General Public License, Version 3.0.

"""Deterministic synthetic paired FASTQ generator for the bqc benchmark.

One scenario per invocation; writes PREFIX_R1.fastq / PREFIX_R2.fastq.
Stdlib only, fixed seed: re-running with the same arguments gives the same
bytes, so the benchmark is reproducible.

Reads model Illumina 3' read-through: when the insert is shorter than the
read, the read continues into the adapter and then into flow-cell poly-A,
the shape real data has.
"""

import argparse

# TruSeq R1/R2 adapters followed by the flow-cell poly-A a read runs into.
# Long enough that insert + read-through always covers a 150 bp read.
R1_READTHROUGH = (
    "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG" + "A" * 100
)
R2_READTHROUGH = (
    "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT" + "A" * 100
)


class Rng:
    """xorshift64*: deterministic and dependency-free."""

    def __init__(self, seed):
        self.state = seed | 1

    def next(self):
        x = self.state
        x ^= x >> 12
        x ^= (x << 25) & 0xFFFFFFFFFFFFFFFF
        x ^= x >> 27
        self.state = x
        return (x * 0x2545F4914F6CDD1D) & 0xFFFFFFFFFFFFFFFF

    def below(self, bound):
        return (self.next() >> 33) % bound

    def seq(self, n):
        return "".join("ACGT"[self.next() >> 62] for _ in range(n))


def damage(adapter, rng):
    """One insertion and one substitution inside the adapter."""
    a = list(adapter)
    a.insert(6, "T")
    pos = 20
    a[pos] = "ACGT".replace(a[pos], "")[rng.below(3)]
    return "".join(a)


def record(rng, name, read_len, insert_len, readthrough, q_low_tail=0):
    insert = rng.seq(min(insert_len, read_len))
    seq = (insert + readthrough + "A" * read_len)[:read_len]
    qual = [38] * read_len
    for i in range(read_len - q_low_tail, read_len):
        qual[i] = 2
    return name, seq, "".join(chr(q + 33) for q in qual)


def pair(rng, i, scenario, read_len=150):
    """Yields (name1, seq1, qual1, name2, seq2, qual2) for one pair."""
    n = f"@SIM:{scenario}:{i}"
    if scenario == "alln":
        seq = "N" * read_len
        q = "G" * read_len
        return n + "/1", seq, q, n + "/2", seq, q

    through1 = through2 = None
    tail = 0
    if scenario == "typical":
        if rng.below(10) < 3:
            through1, through2 = R1_READTHROUGH, R2_READTHROUGH
        if rng.below(10) == 0:
            tail = rng.below(40) + 20
    elif scenario == "dimer":
        through1, through2 = R1_READTHROUGH, R2_READTHROUGH
    elif scenario == "indel":
        through1, through2 = DAMAGED1, DAMAGED2
    elif scenario == "partial":
        k = rng.below(10) + 6  # 6-15 bp adapter stub
        through1, through2 = R1_READTHROUGH[:k], R2_READTHROUGH[:k]
    elif scenario == "polya":
        pass  # handled below
    elif scenario == "lowqual":
        tail = rng.below(40) + 40

    if scenario == "dimer":
        ins1 = ins2 = rng.below(20)
    elif scenario == "partial":
        ins1 = ins2 = read_len - len(through1)  # the read ends mid-adapter
    elif through1 is not None:
        ins1 = ins2 = rng.below(140) + 10
    else:
        ins1 = ins2 = read_len  # no read-through

    r1 = record(rng, n + "/1", read_len, ins1, through1 or "", tail)
    r2 = record(rng, n + "/2", read_len, ins2, through2 or "", tail)

    if scenario == "polya" and rng.below(2) == 0:
        a, g = rng.below(30) + 30, rng.below(30) + 30
        r1 = (r1[0], r1[1][: read_len - a] + "A" * a, r1[2])
        r2 = (r2[0], r2[1][: read_len - g] + "G" * g, r2[2])
    return r1 + r2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "scenario",
        choices=[
            "typical",
            "dimer",
            "indel",
            "partial",
            "polya",
            "lowqual",
            "shortreads",
            "alln",
        ],
    )
    ap.add_argument("-n", "--pairs", type=int, default=1_000_000)
    ap.add_argument("-o", "--out", required=True, help="output prefix")
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    global DAMAGED1, DAMAGED2
    rng = Rng(args.seed)
    DAMAGED1 = damage(R1_READTHROUGH, rng)
    DAMAGED2 = damage(R2_READTHROUGH, rng)

    read_lens = {"shortreads": None}
    with open(args.out + "_R1.fastq", "w") as f1, open(args.out + "_R2.fastq", "w") as f2:
        for i in range(args.pairs):
            rl = rng.below(20) + 60 if args.scenario == "shortreads" else 150
            name1, seq1, qual1, name2, seq2, qual2 = pair(rng, i, args.scenario, rl)
            f1.write(f"{name1}\n{seq1}\n+\n{qual1}\n")
            f2.write(f"{name2}\n{seq2}\n+\n{qual2}\n")


if __name__ == "__main__":
    main()
