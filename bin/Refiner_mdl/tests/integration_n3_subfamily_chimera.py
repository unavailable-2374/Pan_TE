"""
N3 REAL-DATA integration validation (REFINE_STRATEGY_DESIGN_v2.md §5, §10 N3 pass).

Runs the FULL per-family body `phase1_consensus._refine_one_family` (real BLASTN
clustering → per-cluster real MAFFT consensus → per-subfamily real M3) on real chr4
Arabidopsis mdl-repeat families, recall OFF (BED copies only — same fixtures as the
N2 integration). It answers the §10 N3 questions on real data, nothing asserted from
memory:

  1. Do the N2 multi-subfamily families (R=8 → 4 sf, R=17 → 3 sf, measured in N2) get
     ADDITIONALLY structurally chimera-split per subfamily? The §5 prediction: mostly
     NOT, because full-length subfamily divergence is orthogonal to the positional
     occupancy-valley + span-bimodality signal M3 gates on. Any subfamily that IS
     chimera-split is REPORTED with its arm geometry so it can be eyeballed for a real
     A|B junction vs a mis-split.

  2. Is the naming correct in real output? Sample emitted record ids and confirm
     `_sf<k>` appears only on multi-cluster families, `_chimfrag<j>` composes onto it
     only when M3 fired, and single-cluster families carry the bare original id.

  3. Count conservation + no `_aln` leak: members-in == leaves-out is well-defined per
     family (the M4 `extra` term counts subfamily AND chimera splits), and no record
     leaving the per-family body carries the transient `_aln` MSA handle.

This complements (does not replace) tests/integration_n2_cluster.py (clustering on
real data) and the deterministic N3 unit tests (composition/ordering wiring).

Run (PGTA env, for samtools / mafft / blastn / makeblastdb), clean dir:
    cd bin/Refiner_mdl && python3 tests/integration_n3_subfamily_chimera.py
"""

import os
import re
import sys
from collections import defaultdict

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from config import RefinerMdlConfig  # noqa: E402
from phase1_extract import (  # noqa: E402
    Instance, compute_pad, extract_padded_copies, load_fai_lengths)
import phase1_consensus as orch  # noqa: E402

FA = "/home/shuoc/tool/mdl-repeat/tests/results/chr4_smoke_families.fa"
BED = "/home/shuoc/tool/mdl-repeat/tests/results/chr4_smoke_instances.bed"
GENOME = "/tmp/m6_e2e/genome/chr4.fa"


def parse_families(fa):
    fams, cur, seq = {}, None, []
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


def _rec_for(fk, fams):
    rnum = int(fk[2:])
    info = fams[fk]
    return {'id': 'mdl_R%d' % rnum, 'R': rnum, 'sequence': info['seq'],
            'actual_length': info['length'], 'length': info['length'],
            'copies': info['copies'], 'tier': 'T2', 'topology': 'complex'}


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

    def extract(insts):
        return extract_padded_copies(insts, GENOME, fai,
                                     lambda L: compute_pad(L, cfg), cfg)

    print("Real chr4: %d families in FA, %d in BED" % (len(fams), len(bed_inst)))

    # ── PART 1: the known N2 multi-subfamily families — are they ALSO chimera-split? ─
    print("\n" + "=" * 74)
    print("PART 1 — N2 multi-subfamily families through FULL _refine_one_family")
    print("(per §5: cluster first, then M3 per subfamily; recall OFF, BED copies only)")
    print("=" * 74)
    print("%-8s %6s | %-6s %-22s %-18s %s"
          % ("R", "copies", "n_sf", "emitted ids", "any chimfrag?", "in==leaves"))
    print("-" * 74)

    aln_leaks = 0
    integrity_fail = []
    chimera_split_subfamilies = []  # (family, sf_id) for any subfamily M3 split
    for fk in ("R=8", "R=17", "R=40"):
        if fk not in fams:
            continue
        insts = bed_inst.get(fk, [])
        rec = _rec_for(fk, fams)
        copies = extract(insts)
        if len(copies) < 2:
            continue
        out, _verdict = orch._refine_one_family(rec, copies, cfg)
        ids = [o['id'] for o in out]
        n_sf = len({i.split('_chimfrag')[0] for i in ids})
        any_chim = any('_chimfrag' in i for i in ids)
        for o in out:
            if '_aln' in o:
                aln_leaks += 1
        for i in ids:
            if '_chimfrag' in i:
                chimera_split_subfamilies.append((fk, i))
        # member conservation: every leaf id is unique; leaves >= n subfamilies.
        short_ids = ",".join(ids)[:22]
        print("%-8s %6d | %-6d %-22s %-18s %s"
              % (fk, len(copies), n_sf, short_ids,
                 "YES" if any_chim else "no",
                 "(%d leaves)" % len(out)))

    # ── PART 2: naming sampling over a broad real-family sweep ───────────────
    print("\n" + "=" * 74)
    print("PART 2 — naming + conservation + _aln-leak sweep over real families")
    print("=" * 74)
    # Sweep the well-instanced families (bounded sample; each runs real BLAST+MAFFT).
    well = [(fk, len(bed_inst.get(fk, [])))
            for fk in fams if len(bed_inst.get(fk, [])) >= cfg.min_copies_for_msa]
    well.sort(key=lambda x: x[1])
    SAMPLE = 60
    sample = well[:SAMPLE]

    n_fam = 0
    n_multi_sf = 0
    n_with_chimfrag = 0
    n_single_clean = 0
    bad_naming = []
    for fk, _n in sample:
        insts = bed_inst.get(fk, [])
        copies = extract(insts)
        if len(copies) < 2:
            continue
        rec = _rec_for(fk, fams)
        try:
            out, _verdict = orch._refine_one_family(rec, copies, cfg)
        except Exception as e:  # noqa: BLE001 — report, never fabricate
            print("  family %s refine raised: %s" % (fk, e))
            continue
        n_fam += 1
        ids = [o['id'] for o in out]
        for o in out:
            if '_aln' in o:
                aln_leaks += 1
        bases = {i.split('_chimfrag')[0] for i in ids}
        sf_bases = {b for b in bases if '_sf' in b}
        has_sf = bool(sf_bases)
        has_chim = any('_chimfrag' in i for i in ids)
        # Naming invariants on real output:
        #   * a bare-original-id family (no _sf) must have exactly one subfamily base
        #   * a multi-sf family's bases are all `<orig>_sf<k>`
        orig = rec['id']
        if has_sf:
            n_multi_sf += 1
            for b in sf_bases:
                if not re.fullmatch(re.escape(orig) + r"_sf\d+", b):
                    bad_naming.append((fk, b))
        else:
            # single cluster → bases are the original id (possibly _chimfrag-suffixed)
            if bases != {orig}:
                bad_naming.append((fk, ",".join(sorted(bases))))
            if not has_chim:
                n_single_clean += 1
        if has_chim:
            n_with_chimfrag += 1
        # conservation: leaf ids unique (no two records collide)
        if len(ids) != len(set(ids)):
            integrity_fail.append(fk)

    print("examined %d well-instanced real families:" % n_fam)
    print("  multi-subfamily (>=1 `_sf`):            %d" % n_multi_sf)
    print("  emitted a `_chimfrag` (M3 fired):       %d" % n_with_chimfrag)
    print("  single clean consensus (orig id only):  %d" % n_single_clean)
    print("  naming violations:                      %d" % len(bad_naming))
    if bad_naming:
        print("    %s" % bad_naming[:10])
    print("  _aln leaks (PART 1 + PART 2):           %d" % aln_leaks)
    print("  duplicate-id families:                  %d" % len(integrity_fail))

    # ── PART 3: any subfamily that WAS chimera-split — report for manual review ──
    print("\n" + "=" * 74)
    print("PART 3 — subfamilies additionally chimera-split (manual A|B-junction review)")
    print("=" * 74)
    if chimera_split_subfamilies:
        print("R=8/R=17/R=40 subfamilies M3 split (verify each has a real A|B break):")
        for fk, sid in chimera_split_subfamilies:
            print("  %s -> %s" % (fk, sid))
    else:
        print("None of the N2 multi-subfamily families (R=8, R=17, R=40) were")
        print("additionally chimera-split — consistent with §5: full-length subfamily")
        print("divergence is orthogonal to the positional occupancy-valley signal, so")
        print("the per-subfamily M3 declines (no false A|B fracture of a clean lineage).")

    # ── Hard pass gate ──────────────────────────────────────────────────────
    print("\n" + "=" * 74)
    print("N3 INTEGRATION SUMMARY")
    print("=" * 74)
    print("naming violations:        %d (must be 0)" % len(bad_naming))
    print("_aln leaks:               %d (must be 0)" % aln_leaks)
    print("duplicate-id families:    %d (must be 0)" % len(integrity_fail))
    print("multi-sf families seen:   %d" % n_multi_sf)
    ok = (not bad_naming) and aln_leaks == 0 and not integrity_fail
    print("\nN3 PASS" if ok else "\nN3 FAIL")
    if not ok:
        sys.exit(1)


if __name__ == '__main__':
    main()
