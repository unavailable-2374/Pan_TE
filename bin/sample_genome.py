#!/usr/bin/env python3
"""
sample_genome.py — deterministic contiguous-window genome sampler for the large-genome
path of the mdl-repeat / Refiner_mdl pipeline.

For genomes too large to align or run mdl-repeat on whole (wheat 16 Gb, conifer 24 Gb),
draw a representative SAMPLE of random contiguous windows. mdl-repeat targets HIGH-copy
families; a family at genomic copy number C appears ~C·(S/G) times in a sample of size S
from a genome of size G, so at the default 5% a C≥100 family still leaves ≥5 sampled copies
— exactly the recurrence mdl needs. The complete genome is used only later (phase2b) for
the final copy count.

Strategy (matches design_v2/SCALABILITY_STRATEGY.md and the agent landing plan):
  * windows are CONTIGUOUS L-bp blocks (default 1 Mb) — preserves element-level context
    so a full TE sits inside a window, unlike scattered k-mers;
  * drawn WITHOUT replacement (>50%-overlap windows are rejected) and PROPORTIONAL to
    contig length (a 2× longer chromosome gets ~2× the windows);
  * total sampled size  S = clamp(fraction·G, min_bp, max_bp), capped at the samplable G;
  * contigs shorter than --min-contig-bp are skipped (assembly debris / unplaced scaffolds
    add noise, not families);
  * fully DETERMINISTIC: a fixed --seed RNG, so the same genome always yields the same
    sample (G10 reproducibility posture of the rest of the pipeline).

Memory-light: never loads the genome; reads contig lengths from the .fai (built with
samtools faidx if absent) and extracts each window with `samtools faidx genome region`.

Outputs the sample FASTA (windows renamed `<contig>_w<start>_<end>`, 1-based, BLAST-safe)
and, with --bed, a 0-based BED manifest of the sampled intervals for provenance / pre-mask.
"""

import argparse
import bisect
import os
import random
import subprocess
import sys
from collections import defaultdict


def _ensure_fai(genome: str, samtools: str) -> str:
    fai = genome + '.fai'
    if os.path.exists(fai) and os.path.getsize(fai) > 0:
        return fai
    try:
        subprocess.run([samtools, 'faidx', genome], check=True,
                       capture_output=True, text=True)
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        sys.exit(f"ERROR: samtools faidx failed on {genome}: {e}")
    if not os.path.exists(fai):
        sys.exit(f"ERROR: samtools faidx produced no index for {genome}")
    return fai


def _read_lengths(fai: str):
    """Return [(name, length)] from a .fai, preserving file order."""
    out = []
    with open(fai) as fh:
        for line in fh:
            p = line.rstrip('\n').split('\t')
            if len(p) >= 2:
                out.append((p[0], int(p[1])))
    return out


def _select_windows(lengths, L, S, min_contig, seed):
    """Pick contiguous windows proportional to contig length until ~S bp, deterministic.

    Returns (windows, samplable_bp) where windows = [(contig, start1, end1)] (1-based,
    inclusive). Windows overlapping an existing one by >50% are rejected (sampling without
    replacement); if S exceeds the samplable total, every eligible contig is taken whole.
    """
    eligible = [(n, l) for n, l in lengths if l >= min_contig]
    samplable = sum(l for _, l in eligible)
    if samplable == 0:
        return [], 0

    # S can't exceed what's available; if it (nearly) does, just take everything.
    if S >= samplable:
        return [(n, 1, l) for n, l in eligible], samplable

    rng = random.Random(seed)
    # Cumulative-length table for O(log n) length-weighted contig choice.
    names = [n for n, _ in eligible]
    lens = {n: l for n, l in eligible}
    cum = []
    acc = 0
    for _, l in eligible:
        acc += l
        cum.append(acc)
    total_w = cum[-1]

    used = defaultdict(list)   # contig -> [(s1,e1)]
    windows = []
    total = 0
    # Generous attempt ceiling so the without-replacement rejection can't spin forever
    # on a nearly-saturated contig set; relax the overlap rule in the final 10%.
    max_windows = int(S / L) + len(eligible) + 8
    max_attempts = max_windows * 20

    def _overlap_frac(c, s1, e1):
        worst = 0.0
        for us, ue in used[c]:
            lo, hi = max(s1, us), min(e1, ue)
            if hi >= lo:
                worst = max(worst, (hi - lo + 1) / (e1 - s1 + 1))
        return worst

    attempts = 0
    while total < S and len(windows) < max_windows and attempts < max_attempts:
        attempts += 1
        r = rng.random() * total_w
        c = names[bisect.bisect_right(cum, r)]
        clen = lens[c]
        if clen <= L:
            s1, e1 = 1, clen
        else:
            s0 = rng.randint(0, clen - L)   # 0-based start
            s1, e1 = s0 + 1, s0 + L
        strict = attempts < max_attempts * 0.9
        if strict and _overlap_frac(c, s1, e1) > 0.5:
            continue
        used[c].append((s1, e1))
        windows.append((c, s1, e1))
        total += (e1 - s1 + 1)
    return windows, samplable


def _extract(genome, windows, out_fa, samtools, batch=500):
    """Extract windows via samtools faidx (batched), renaming headers to BLAST-safe ids."""
    n_written = 0
    with open(out_fa, 'w') as out:
        for i in range(0, len(windows), batch):
            chunk = windows[i:i + batch]
            regions = [f"{c}:{s}-{e}" for c, s, e in chunk]
            try:
                r = subprocess.run([samtools, 'faidx', genome] + regions,
                                   check=True, capture_output=True, text=True)
            except (subprocess.CalledProcessError, FileNotFoundError) as e:
                sys.exit(f"ERROR: samtools faidx extraction failed: {e}")
            # samtools emits records in request order; rename by position in `chunk`.
            idx = -1
            for line in r.stdout.splitlines():
                if line.startswith('>'):
                    idx += 1
                    c, s, e = chunk[idx]
                    out.write(f">{c}_w{s}_{e}\n")
                    n_written += 1
                else:
                    out.write(line + '\n')
    return n_written


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('-g', '--genome', required=True, help='Full genome FASTA (indexed or indexable)')
    ap.add_argument('-o', '--out', required=True, help='Output sample FASTA')
    ap.add_argument('--bed', default='', help='Optional 0-based BED manifest of sampled windows')
    ap.add_argument('--chunk-len-bp', type=int, default=1_000_000, help='Window length L (default 1 Mb)')
    ap.add_argument('--fraction', type=float, default=0.05, help='Target sample fraction of samplable G (default 0.05)')
    ap.add_argument('--min-bp', type=int, default=200_000_000, help='Lower clamp on sample size S (default 200 Mb)')
    ap.add_argument('--max-bp', type=int, default=1_000_000_000, help='Upper clamp on sample size S (default 1 Gb)')
    ap.add_argument('--min-contig-bp', type=int, default=100_000, help='Skip contigs shorter than this (default 100 kb)')
    ap.add_argument('--seed', type=int, default=42, help='Deterministic RNG seed (default 42)')
    ap.add_argument('--samtools', default='samtools', help='samtools executable')
    args = ap.parse_args()

    if not os.path.exists(args.genome):
        sys.exit(f"ERROR: genome not found: {args.genome}")

    fai = _ensure_fai(args.genome, args.samtools)
    lengths = _read_lengths(fai)
    G = sum(l for _, l in lengths)
    eligible_bp = sum(l for _, l in lengths if l >= args.min_contig_bp)

    S = max(args.min_bp, min(args.max_bp, int(args.fraction * eligible_bp)))
    S = min(S, eligible_bp)

    windows, samplable = _select_windows(
        lengths, args.chunk_len_bp, S, args.min_contig_bp, args.seed)
    if not windows:
        sys.exit("ERROR: no samplable contigs (all below --min-contig-bp?); cannot sample")

    n = _extract(args.genome, windows, args.out, args.samtools)
    sampled_bp = sum(e - s + 1 for _, s, e in windows)

    if args.bed:
        with open(args.bed, 'w') as bf:
            for c, s, e in windows:
                bf.write(f"{c}\t{s-1}\t{e}\n")

    print(f"[sample_genome] genome G={G/1e9:.2f} Gb, samplable={samplable/1e9:.2f} Gb "
          f"(contigs >= {args.min_contig_bp/1e3:.0f} kb)", file=sys.stderr)
    print(f"[sample_genome] target S={S/1e6:.0f} Mb (frac {args.fraction}, "
          f"clamp [{args.min_bp/1e6:.0f},{args.max_bp/1e6:.0f}] Mb) -> "
          f"{n} windows of {args.chunk_len_bp/1e6:.1f} Mb, {sampled_bp/1e6:.0f} Mb sampled "
          f"({100*sampled_bp/G:.1f}% of genome) -> {args.out}", file=sys.stderr)


if __name__ == '__main__':
    main()
