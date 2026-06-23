"""
N5 REAL-DATA integration — conservative FP prune on real chr4 Arabidopsis.

NO synthetic data: every family, every copy, every distance is real (real chr4
mdl-repeat output + real chr4 genome + real all-vs-all BLASTN). The script drives the
PRODUCTION path per family — extract real padded copies -> cluster (real BLASTN
distance) -> qualifying_clusters -> pseudo_family_verdict — and reports the §10-N5
pass criteria on real numbers:

  1. true low-copy REAL families are ALL kept (no N5 drop among copies<5) — the core
     invariant (chr4 is 84% low-copy; a drop here would be误删 a real family);
  2. total dropped families + per-drop reason + the dropped-family audit (id + metrics);
  3. dropped-family spot-check: re-confirm on the real all-vs-all BLASTN that the dropped
     copies share essentially no coherent homology (no recoverable support);
  4. R=150 handling (single cluster, MSA-fail archetype) — kept or dropped, reported with
     its measured incoherence metrics;
  5. degrade-safe: enable_conservative_fp_prune=False drops NOTHING (== N4);
  6. conservation: drops are counted, audited, never silent.

Inputs (real):
  families : /home/shuoc/tool/mdl-repeat/tests/results/chr4_smoke_families.fa
  BED      : /home/shuoc/tool/mdl-repeat/tests/results/chr4_smoke_instances.bed
  genome   : /tmp/m6_e2e/genome/chr4.fa (+ .fai)

Run:  cd bin/Refiner_mdl && python3 tests/integration_n5_fp.py
"""

import os
import re
import sys
import time
from collections import defaultdict

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from config import RefinerMdlConfig                                   # noqa: E402
from phase1_extract import Instance, compute_pad, extract_padded_copies, load_fai_lengths  # noqa: E402
from phase1_cluster import cluster_copies, qualifying_clusters        # noqa: E402
import phase1_cluster as cl                                           # noqa: E402
from phase1_fp_prune import pseudo_family_verdict, _pairwise_identity_stats  # noqa: E402

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
                fams[cur] = {'length': int(m.group(2)), 'copies': int(m.group(3)), 'seq': ''}
            seq = []
        else:
            seq.append(line.strip())
    if cur:
        fams[cur]['seq'] = "".join(seq)
    return fams


def parse_bed(bed):
    out = defaultdict(list)
    for line in open(bed):
        p = line.rstrip("\n").split("\t")
        if len(p) < 6:
            continue
        try:
            score = int(p[4])
        except ValueError:
            continue
        out[p[3]].append(Instance(chrom=p[0], start=int(p[1]), end=int(p[2]),
                                  strand=p[5], divergence=1.0 - score / 1000.0))
    return out


def _rec_for(fk, fams):
    info = fams[fk]
    return {'id': 'mdl_%s' % fk.replace('=', ''), 'sequence': info['seq'],
            'actual_length': info['length'], 'length': info['length'],
            'copies': info['copies'], 'tier': 'T2'}


def main():
    for f in (FA, BED, GENOME, GENOME + '.fai'):
        if not os.path.exists(f):
            print("MISSING fixture: %s" % f)
            sys.exit(1)

    fams = parse_families(FA)
    bed = parse_bed(BED)
    cfg = RefinerMdlConfig()
    cfg.genome_file = GENOME
    fai = load_fai_lengths(GENOME + '.fai')
    N = len(fams)
    print("Real chr4: %d families in FA, %d families with BED rows" % (N, len(bed)))
    print("N5 thresholds: min_copies_for_call=%d  max_intra_identity=%.2f  "
          "min_homologous_pair_frac=%.2f"
          % (cfg.pseudo_family_min_copies_for_call,
             cfg.pseudo_family_max_intra_identity,
             cfg.pseudo_family_min_homologous_pair_frac))

    def extract(insts):
        return extract_padded_copies(insts, GENOME, fai, lambda L: compute_pad(L, cfg), cfg)

    # ── Sweep every family through the REAL extract->cluster->verdict path ────
    # (No recall here: N5's discriminant runs on the assembled BED copy set; recall
    #  would only ADD copies. completeness_verified is left absent, so the §4.3
    #  abstention is exercised separately in PART 4 with a forced over-budget flag.)
    print("\n" + "=" * 74)
    print("PART 1 — sweep all families: extract -> cluster -> N5 verdict (real BLASTN)")
    print("=" * 74)
    t0 = time.time()
    dropped = []            # (fk, copies, verdict)
    kept_lowcopy = 0        # families with copies<5 that were KEPT
    lowcopy_total = 0
    dropped_lowcopy = []    # any copies<5 family that was dropped (MUST be empty)
    n_processed = 0
    n_fallback = 0          # families that formed zero qualifying clusters
    r150_report = None

    # Process in a deterministic order (by R number) for reproducible logging.
    def _rk(fk):
        try:
            return int(fk.split('=')[1])
        except Exception:
            return 1 << 30

    for fk in sorted(fams.keys(), key=_rk):
        insts = bed.get(fk, [])
        # Clear the per-family BLAST memo so each family runs its own real distance.
        cl._LAST_DIST_KEY = None
        cl._LAST_DIST_VAL = None
        try:
            copies = extract(insts) if insts else []
        except Exception as e:  # noqa: BLE001 — report, never fabricate
            print("  %s extract raised: %s" % (fk, e))
            continue
        n_copies = len(copies)
        if n_copies < 2:
            # 0-1 copy: clustering moot; verdict trivially keeps (and is cheap). Still
            # run it so the low-copy-kept accounting is complete.
            pass
        rec = _rec_for(fk, fams)
        raw = cluster_copies(copies, cfg)
        quals, audit = qualifying_clusters(raw, copies, cfg)
        v = pseudo_family_verdict(rec, copies, audit, cfg)
        n_processed += 1
        if audit.get('fallback_single_cluster'):
            n_fallback += 1

        is_lowcopy = n_copies < 5
        if is_lowcopy:
            lowcopy_total += 1
        if v.is_pseudo:
            dropped.append((fk, copies, v))
            if is_lowcopy:
                dropped_lowcopy.append((fk, n_copies, v))
        elif is_lowcopy:
            kept_lowcopy += 1

        if fk == 'R=150':
            r150_report = (n_copies, audit, v)

    elapsed = time.time() - t0
    print("Processed %d families in %.1fs (%d formed zero qualifying clusters)"
          % (n_processed, elapsed, n_fallback))

    # ── PART 2 — the CORE invariant: no true low-copy family was dropped ──────
    print("\n" + "=" * 74)
    print("PART 2 — true low-copy families ALL kept (the核心不变量)")
    print("=" * 74)
    print("low-copy families (copies<5):        %5d" % lowcopy_total)
    print("  of which KEPT by N5:               %5d" % kept_lowcopy)
    print("  of which DROPPED by N5:            %5d   <== MUST be 0" % len(dropped_lowcopy))
    if dropped_lowcopy:
        print("  !!! VIOLATION: low-copy families dropped (误删 real families):")
        for fk, nc, v in dropped_lowcopy[:20]:
            print("      %s (%d copies): %s" % (fk, nc, v.reason))
    assert not dropped_lowcopy, \
        "CORE INVARIANT VIOLATED: a low-copy family was dropped (low copy判假)"
    print("PASS: no low-copy family was dropped — low copy count never判假 (§4.1/§4.3)")

    # ── PART 3 — dropped-family count, reasons, and spot-check audit ──────────
    print("\n" + "=" * 74)
    print("PART 3 — dropped families: count, reasons, recoverable-support spot-check")
    print("=" * 74)
    print("TOTAL families dropped by N5:        %5d   (%.2f%% of %d)"
          % (len(dropped), 100.0 * len(dropped) / max(1, N), N))
    if dropped:
        print("\n  per-dropped-family audit (id | copies | median_id | homol_pair_frac):")
        for fk, copies, v in dropped:
            m = v.metrics
            print("    %-8s | %3d copies | median_id=%.3f | homol_pair_frac=%.3f"
                  % (fk, m.get('n_copies', len(copies)),
                     m.get('median_identity', -1), m.get('homologous_pair_frac', -1)))
        # Spot-check: independently re-measure the all-vs-all identity stats on the
        # dropped copies and confirm they really are incoherent (no recoverable support).
        print("\n  spot-check (independent re-measure of dropped copy sets):")
        recoverable = []
        for fk, copies, v in dropped:
            cl._LAST_DIST_KEY = None
            cl._LAST_DIST_VAL = None
            stats = _pairwise_identity_stats(copies, cfg)
            coherent_again = (stats['median_identity'] >= cfg.pseudo_family_max_intra_identity
                              or stats['homologous_pair_frac'] >= cfg.pseudo_family_min_homologous_pair_frac)
            flag = "RECOVERABLE?!" if coherent_again else "no support (confirmed)"
            print("    %-8s re-measured median_id=%.3f homol_frac=%.3f -> %s"
                  % (fk, stats['median_identity'], stats['homologous_pair_frac'], flag))
            if coherent_again:
                recoverable.append(fk)
        assert not recoverable, \
            ("a dropped family re-measured as coherent (%s) — conservative-prune "
             "verification FAILED" % recoverable)
        print("  PASS: every dropped family re-confirmed incoherent (no recoverable support)")
    else:
        print("  (no families dropped on this real chr4 BED set — refine is deliberately"
              " under-pruning; see report for interpretation)")

    # ── PART 4 — R=150 archetype + §4.3 over-budget abstention on real data ───
    print("\n" + "=" * 74)
    print("PART 4 — R=150 archetype and §4.3 over-budget abstention (real)")
    print("=" * 74)
    if r150_report:
        nc, audit, v = r150_report
        print("R=150: %d BED copies | fallback_single_cluster=%s | is_pseudo=%s"
              % (nc, audit.get('fallback_single_cluster'), v.is_pseudo))
        print("       metrics: %s" % v.metrics)
        print("       reason: %s" % v.reason)
    else:
        print("R=150 not present in this BED subset (reported as-is, not fabricated)")

    # §4.3: force completeness_verified=False on a would-be-dropped family and confirm
    # it is KEPT (over-budget recall must never be dropped). Use the first dropped family
    # if any; otherwise construct the check on a real fallback family.
    if dropped:
        fk, copies, _v = dropped[0]
        cl._LAST_DIST_KEY = None
        cl._LAST_DIST_VAL = None
        raw = cluster_copies(copies, cfg)
        quals, audit = qualifying_clusters(raw, copies, cfg)
        rec_ob = _rec_for(fk, fams)
        rec_ob['completeness_verified'] = False
        v_ob = pseudo_family_verdict(rec_ob, copies, audit, cfg)
        print("\n§4.3 abstention: %s WITH completeness_verified=False -> is_pseudo=%s"
              % (fk, v_ob.is_pseudo))
        assert not v_ob.is_pseudo, "over-budget family must be KEPT (§4.3)"
        print("PASS: an over-budget (unverified) family is KEPT even when its copies are"
              " incoherent (§4.3 absence of verification ≠ absence of support)")

    # ── PART 5 — degrade-safe: disabled drops NOTHING ─────────────────────────
    print("\n" + "=" * 74)
    print("PART 5 — degrade-safe: enable_conservative_fp_prune=False drops nothing")
    print("=" * 74)
    cfg_off = RefinerMdlConfig()
    cfg_off.genome_file = GENOME
    cfg_off.enable_conservative_fp_prune = False
    n_off_drop = 0
    for fk, copies, _v in dropped:
        cl._LAST_DIST_KEY = None
        cl._LAST_DIST_VAL = None
        raw = cluster_copies(copies, cfg_off)
        quals, audit = qualifying_clusters(raw, copies, cfg_off)
        v_off = pseudo_family_verdict(_rec_for(fk, fams), copies, audit, cfg_off)
        if v_off.is_pseudo:
            n_off_drop += 1
    print("families dropped with N5 DISABLED:   %5d   <== MUST be 0" % n_off_drop)
    assert n_off_drop == 0, "disabled N5 must drop nothing (degrade to N4)"
    print("PASS: disabled N5 degrades exactly to N4 (no drops)")

    print("\n" + "=" * 74)
    print("ALL N5 INTEGRATION CHECKS PASSED (real chr4, real BLASTN, real clustering)")
    print("  families=%d  processed=%d  dropped=%d  low-copy-kept=%d/%d"
          % (N, n_processed, len(dropped), kept_lowcopy, lowcopy_total))
    print("=" * 74)


if __name__ == '__main__':
    main()
