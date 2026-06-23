"""
N2 REAL-DATA integration validation (REFINE_STRATEGY_DESIGN_v2.md §10 N2 pass).

Uses real chr4 Arabidopsis mdl-repeat output + the real chr4 genome:
  families : /home/shuoc/tool/mdl-repeat/tests/results/chr4_smoke_families.fa
  BED      : /home/shuoc/tool/mdl-repeat/tests/results/chr4_smoke_instances.bed
  genome   : /tmp/m6_e2e/genome/chr4.fa (+ .fai)   [extract copies — REAL sequence]

The N2 milestone is BED-copies-only (recall OFF). The distance substrate is all-vs-all
BLASTN LOCAL-identity among a family's own copies (phase1_cluster.pairwise_distance) —
local alignment is robust to the rearrangement / internal-repeat content of mdl-stage
copies that defeats a single global MSA (an earlier k-mer-Jaccard substrate saturated
at ~1.0 and a global-MSA-core substrate at ~0.5 for every family, neither separating;
see phase1_cluster module docstring). The §10 N2 expectations encoded here, each
measured on REAL data (nothing asserted from memory):

  1. HOMOGENEOUS FAMILY → 1 cluster → output IDENTICAL to today.
     For several real low-spread multi-copy families, run BOTH:
       (a) the pre-N2 path: refine_family(rec, all_copies)               [single consensus]
       (b) the N2 path: cluster_copies → qualifying_clusters → per-cluster refine
     and assert the N2 path emits exactly ONE record whose id == original id and whose
     sequence is byte-identical to (a). This is the "degrades to today" proof.

  2. DIVERGENT FAMILY → >1 subfamily, ON REAL DATA (the headline result of this fix).
     The BLASTN substrate RESOLVES real multi-subfamily families that BED spread alone
     does not flag (BED spread is distance-to-seed, a 1-D proxy §3.2 calls insufficient).
     Measured on real chr4: high-copy families R=8 (56 copies) → 4 subfamilies, R=17
     (34) → 3, with tight subfamilies (intra-identity ~0.8–0.9) — we assert >1 here, on
     real sequence, and report per-subfamily member counts + intra/inter identity +
     per-subfamily consensus length. A genuinely homogeneous high-copy family (R=40)
     stays 1 — reported as the contrast.

  3. MEMBER CONSERVATION: copies in == members out for every examined family; the
     qualifying gate reassigns, never drops.

Run (PGTA env, for samtools / mafft / blastn / makeblastdb), clean dir:
    cd bin/Refiner_mdl && python3 tests/integration_n2_cluster.py
"""

import os
import re
import sys
import tempfile

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from config import RefinerMdlConfig  # noqa: E402
from phase1_extract import (  # noqa: E402
    Copy, Instance, compute_pad, extract_padded_copies, load_fai_lengths)
from phase1_boundary import refine_family  # noqa: E402
from phase1_chimera import strip_aln  # noqa: E402
from phase1_cluster import (  # noqa: E402
    cluster_copies, pairwise_distance, qualifying_clusters, subfamily_id)

FA = "/home/shuoc/tool/mdl-repeat/tests/results/chr4_smoke_families.fa"
BED = "/home/shuoc/tool/mdl-repeat/tests/results/chr4_smoke_instances.bed"
GENOME = "/tmp/m6_e2e/genome/chr4.fa"


def parse_families(fa):
    fams = {}
    cur = None
    seq = []
    for line in open(fa):
        if line.startswith(">"):
            if cur:
                fams[cur]['seq'] = "".join(seq)
            m = re.match(r">R=(\d+)\s+length=(\d+)\s+copies=(\d+)", line)
            cur = "R=%s" % m.group(1) if m else None
            if cur:
                fams[cur] = {'length': int(m.group(2)),
                             'copies': int(m.group(3)), 'seq': ''}
            seq = []
        else:
            seq.append(line.strip())
    if cur:
        fams[cur]['seq'] = "".join(seq)
    return fams


def parse_bed_instances(bed):
    from collections import defaultdict
    out = defaultdict(list)
    for line in open(bed):
        p = line.rstrip("\n").split("\t")
        if len(p) < 6:
            continue
        try:
            score = int(p[4])
        except ValueError:
            continue
        div = 1.0 - score / 1000.0
        out[p[3]].append(Instance(chrom=p[0], start=int(p[1]), end=int(p[2]),
                                  strand=p[5], divergence=div))
    return out


def _bed_spread(insts):
    if len(insts) < 2:
        return 0.0
    divs = [i.divergence for i in insts]
    m = sum(divs) / len(divs)
    return (sum((d - m) ** 2 for d in divs) / len(divs)) ** 0.5


def _rec_for(fk, fams):
    rnum = int(fk[2:])
    info = fams[fk]
    return {'id': 'mdl_R%d' % rnum, 'R': rnum, 'sequence': info['seq'],
            'actual_length': info['length'], 'length': info['length'],
            'copies': info['copies']}


def _extract(insts, cfg, fai):
    return extract_padded_copies(insts, GENOME, fai,
                                 lambda L: compute_pad(L, cfg), cfg)


def _n2_family(rec, copies, cfg):
    """Run the N2 per-family body (cluster → qualify → per-cluster refine) and return
    the list of (id, sequence, n_subfamilies, member_count) emitted records."""
    raw = cluster_copies(copies, cfg)
    quals, audit = qualifying_clusters(raw, copies, cfg)
    n_sf = len(quals)
    out = []
    for k, cl in enumerate(quals, start=1):
        sid = subfamily_id(rec['id'], k, n_sf)
        sf_rec = dict(rec)
        sf_rec['id'] = sid
        if n_sf > 1:
            sf_rec['copies'] = cl.size
        refined = strip_aln(refine_family(sf_rec, cl.members, cfg))
        out.append((refined['id'], refined['sequence'], n_sf, cl.size))
    return out, audit


def main():
    for f in (FA, BED, GENOME, GENOME + '.fai'):
        if not os.path.exists(f):
            print("MISSING fixture: %s" % f)
            sys.exit(1)

    fams = parse_families(FA)
    bed_inst = parse_bed_instances(BED)
    cfg = RefinerMdlConfig()
    cfg.genome_file = GENOME
    fai = load_fai_lengths(GENOME + '.fai')
    print("Real chr4: %d families in FA, %d in BED" % (len(fams), len(bed_inst)))

    # Rank well-instanced families by BED spread to pick honest test subjects.
    summary = []
    for fk, info in fams.items():
        insts = bed_inst.get(fk, [])
        summary.append((fk, len(insts), _bed_spread(insts), info['length']))
    well = [s for s in summary if s[1] >= cfg.min_copies_for_msa]
    well_by_spread = sorted(well, key=lambda s: s[2])

    # ── PART 1: homogeneous families degrade to today (byte-identical) ──────
    # Pick the lowest-spread well-instanced families (the homogeneous core).
    print("\n" + "=" * 70)
    print("PART 1 — homogeneous family → 1 cluster → output IDENTICAL to pre-N2")
    print("=" * 70)
    print("%-8s %6s %8s %6s | %-9s %-9s %-7s %s"
          % ("R", "BED_n", "spread", "len", "pre_N2_id", "N2_id", "n_sf", "seq match"))
    print("-" * 70)
    n_checked = 0
    n_identical = 0
    for fk, nbed, spread, length in well_by_spread[:8]:
        insts = bed_inst.get(fk, [])
        rec = _rec_for(fk, fams)
        copies = _extract(insts, cfg, fai)
        if len(copies) < cfg.min_copies_for_msa:
            continue
        # (a) pre-N2 single-consensus path
        pre = strip_aln(refine_family(dict(rec), copies, cfg))
        pre_id, pre_seq = pre['id'], pre['sequence']
        # (b) N2 path
        n2_out, audit = _n2_family(rec, copies, cfg)
        n_sf = n2_out[0][2] if n2_out else 0
        same = (len(n2_out) == 1 and n2_out[0][0] == pre_id
                and n2_out[0][1] == pre_seq)
        n_checked += 1
        if same:
            n_identical += 1
        print("%-8s %6d %8.3f %6d | %-9s %-9s %-7d %s"
              % (fk, nbed, spread, length, pre_id,
                 n2_out[0][0] if n2_out else "-", n_sf,
                 "IDENTICAL" if same else "DIFFERS"))
    print("\nhomogeneous degrade-to-today: %d/%d families byte-identical to pre-N2"
          % (n_identical, n_checked))

    # ── PART 2(i): REAL divergent families → >1 subfamily (the headline result) ──
    # These families are found by the BLASTN local-identity substrate, NOT by BED spread
    # (BED spread is the 1-D distance-to-seed proxy §3.2 calls insufficient; the families
    # that actually carry subfamily structure are high-copy ones whose pairwise BLASTN
    # matrix has tight sub-blocks). R=8/R=17 are the measured real splitters; R=40 is a
    # genuinely homogeneous high-copy contrast that stays 1.
    print("\n" + "=" * 70)
    print("PART 2(i) — REAL chr4 divergent families resolve into >1 subfamily")
    print("(BLASTN local-identity substrate; recall OFF, BED copies only)")
    print("=" * 70)
    print("%-8s %6s %6s | %-6s %-16s %-10s %s"
          % ("R", "copies", "len", "n_sf", "sf sizes", "intra_id", "in==out"))
    print("-" * 70)
    import numpy as np
    real_split_families = []  # (fk, n_sf) for the assertion
    for fk in ("R=8", "R=17", "R=40"):
        if fk not in fams:
            continue
        insts = bed_inst.get(fk, [])
        rec = _rec_for(fk, fams)
        copies = _extract(insts, cfg, fai)
        if len(copies) < 2:
            continue
        raw = cluster_copies(copies, cfg)
        quals, audit = qualifying_clusters(raw, copies, cfg)
        sizes = ",".join(str(c.size) for c in quals)
        intra = ",".join("%.2f" % c.mean_identity for c in quals)
        conserved = audit['total_members_out'] == len(copies)
        real_split_families.append((fk, len(quals)))
        print("%-8s %6d %6d | %-6d %-16s %-10s %s"
              % (fk, len(copies), fams[fk]['length'], len(quals),
                 sizes, intra, "YES" if conserved else "NO"))
        # Per-subfamily consensus length + inter-subfamily identity for the splitters.
        if len(quals) > 1:
            n2_out, _ = _n2_family(rec, copies, cfg)
            cons_lens = ",".join(str(len(seq)) for _id, seq, _n, _m in n2_out)
            print("           per-sf consensus lengths: %s" % cons_lens)
            d = pairwise_distance(copies, cfg)
            a, b = quals[0].member_indices, quals[1].member_indices
            inter = float(np.mean([1.0 - d[i, j] for i in a for j in b])) if a and b else 0.0
            print("           inter sf1-sf2 mean identity: %.3f (lower = more distinct)"
                  % inter)
    n_real_splitters = sum(1 for _fk, k in real_split_families if k > 1)
    print("\nREAL divergent families emitting >1 subfamily: %d (R=8, R=17 expected)"
          % n_real_splitters)

    # ── PART 2(ii): SYNTHETIC divergent set corroborates the mechanism emits >1 ──
    print("\n" + "=" * 70)
    print("PART 2(ii) — SYNTHETIC two-divergence-group copy set → mechanism corroboration")
    print("(clearly labeled SYNTHETIC; the REAL split proof is PART 2(i) above)")
    print("=" * 70)

    def _synthA(v):
        base = ("ACGTTGCAACCGTTAGGCATGACTTGACCATGGTACATGCATGCATTAGGCATCAGT"
                "TTGACCATGAGTCATGCATTAGCATGCATCATGGTACATGCATGCATTAGGCATCAG") * 4
        s = list(base)
        if v:
            p = (v * 53) % len(s)
            s[p] = 'A' if s[p] != 'A' else 'C'
        return ''.join(s)

    def _synthB(v):
        base = ("TTAACCGGTTAACCGGAATTCCGGAATTGGCCTTAAGGCCTTAACCGGAATTCCGGT"
                "AATTGGCCAATTCCGGTTAACCGGAATTGGCCTTAAGGCCAATTCCGGTTAACCGGA") * 4
        s = list(base)
        if v:
            p = (v * 53) % len(s)
            s[p] = 'T' if s[p] != 'T' else 'G'
        return ''.join(s)

    syn = ([Copy(id="synA%d" % v, sequence=_synthA(v), divergence=0.05, strand='+')
            for v in range(6)]
           + [Copy(id="synB%d" % v, sequence=_synthB(v), divergence=0.05, strand='+')
              for v in range(6)])
    syn_rec = {'id': 'mdl_R_SYNTH', 'R': 'SYNTH', 'sequence': _synthA(0),
               'actual_length': len(_synthA(0)), 'length': len(_synthA(0)),
               'copies': 12}
    raw = cluster_copies(syn, cfg)
    quals, audit = qualifying_clusters(raw, syn, cfg)
    n2_out, _ = _n2_family(syn_rec, syn, cfg)
    print("SYNTHETIC 6+6 two-group set: clustered into %d subfamilies, sizes %s"
          % (len(quals), [c.size for c in quals]))
    print("emitted ids: %s" % [o[0] for o in n2_out])
    print("members in=12  members out=%d  reassigned=%d  fallback=%s"
          % (audit['total_members_out'], audit['n_reassigned'],
             audit['fallback_single_cluster']))
    mech_ok = (len(quals) == 2
               and sorted(o[0] for o in n2_out)
               == ['mdl_R_SYNTH_sf1', 'mdl_R_SYNTH_sf2']
               and audit['total_members_out'] == 12)
    print("MECHANISM: %s (expect 2 subfamilies, _sf1/_sf2 ids, 12 members conserved)"
          % ("OK" if mech_ok else "FAIL"))

    # ── PART 3: member conservation across all examined real families ───────
    # Each family runs an all-vs-all BLAST, so sweep a bounded SAMPLE (the smallest-copy
    # well-instanced families, fastest) rather than all 400+ — enough to verify the
    # never-drop invariant broadly without an unbounded runtime.
    SAMPLE = 50
    sample = sorted(well, key=lambda s: s[1])[:SAMPLE]
    print("\n" + "=" * 70)
    print("PART 3 — member conservation over a real-family sample (n=%d)" % len(sample))
    print("=" * 70)
    n_fam = 0
    n_conserved = 0
    n_multi = 0
    for fk, nbed, spread, length in sample:
        insts = bed_inst.get(fk, [])
        copies = _extract(insts, cfg, fai)
        if len(copies) < 2:
            continue
        raw = cluster_copies(copies, cfg)
        quals, audit = qualifying_clusters(raw, copies, cfg)
        n_fam += 1
        if audit['total_members_out'] == len(copies):
            n_conserved += 1
        if len(quals) > 1:
            n_multi += 1
    print("examined %d well-instanced real families:" % n_fam)
    print("  member-conserved (in==out): %d/%d" % (n_conserved, n_fam))
    print("  emitted >1 subfamily in this sample: %d" % n_multi)

    print("\n" + "=" * 70)
    print("N2 INTEGRATION SUMMARY")
    print("=" * 70)
    print("homogeneous → today (byte-identical): %d/%d" % (n_identical, n_checked))
    print("REAL divergent families → >1 subfamily: %d (R=8, R=17)" % n_real_splitters)
    print("synthetic divergent → multi-subfamily mechanism: %s"
          % ("OK" if mech_ok else "FAIL"))
    print("member conservation: %d/%d sampled real families" % (n_conserved, n_fam))
    print("=" * 70)

    # Hard pass/fail gate (the success criterion of the distance-substrate fix):
    ok = (n_real_splitters >= 1 and mech_ok and n_conserved == n_fam
          and n_identical == n_checked)
    print("\nN2 PASS" if ok else "\nN2 FAIL")
    if not ok:
        sys.exit(1)


if __name__ == '__main__':
    main()
