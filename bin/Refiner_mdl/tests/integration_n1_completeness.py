"""
N1 REAL-DATA integration validation (REFINE_STRATEGY_DESIGN_v2.md §10 N1 pass).

Uses real chr4 Arabidopsis mdl-repeat output:
  families : /home/shuoc/tool/mdl-repeat/tests/results/chr4_smoke_families.fa
  BED      : /home/shuoc/tool/mdl-repeat/tests/results/chr4_smoke_instances.bed
  genome DB: /tmp/m6_e2e/refiner_output/temp_work/genome_blastdb/genome (blastn
             ground-truth ONLY — never consulted by the N1 flag decision)

Steps:
  1. Parse families + run the extended write_instance_index to emit the real
     completeness_presignals.tsv from a single BED pass.
  2. Run assess_completeness on ALL surviving families; report flagged count/fraction
     and per-signal breakdown.  This decision touches NO genome (verified).
  3. For a sample (R=150 known divergent-tail + homogeneous controls), independently
     measure the divergent tail with blastn and tabulate flag vs ground-truth.

Run:  cd bin/Refiner_mdl && python3 tests/integration_n1_completeness.py
(requires PGTA conda env active for blastn)
"""

import os
import re
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from config import RefinerMdlConfig  # noqa: E402
from phase1_extract import Instance  # noqa: E402
import phase0_triage  # noqa: E402
from phase1_completeness import assess_completeness, load_presignals  # noqa: E402

FA = "/home/shuoc/tool/mdl-repeat/tests/results/chr4_smoke_families.fa"
BED = "/home/shuoc/tool/mdl-repeat/tests/results/chr4_smoke_instances.bed"
DB = "/tmp/m6_e2e/refiner_output/temp_work/genome_blastdb/genome"


def parse_families(fa):
    """{R=N: {'seq':..., 'length':int, 'copies':int}} from mdl-repeat FASTA."""
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
    """{R=N: [Instance,...]} from the real BED."""
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


def blastn_tail(seq, db, fam):
    """Return (n_hits>=70, n_divergent_70_90) for one family's consensus."""
    with tempfile.NamedTemporaryFile('w', suffix='.fa', delete=False) as fh:
        fh.write(">%s\n%s\n" % (fam, seq))
        q = fh.name
    try:
        cmd = ['blastn', '-query', q, '-db', db,
               '-outfmt', '6 pident qcovhsp', '-evalue', '1e-5',
               '-max_target_seqs', '500', '-num_threads', '4']
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        h70 = 0
        div = 0
        for line in r.stdout.splitlines():
            parts = line.split('\t')
            if len(parts) < 2:
                continue
            pid = float(parts[0])
            if pid >= 70:
                h70 += 1
                if pid < 90:
                    div += 1
        return h70, div
    finally:
        os.unlink(q)


def main():
    fams = parse_families(FA)
    bed_inst = parse_bed_instances(BED)
    print("Real chr4: %d families in FA, %d in BED" % (len(fams), len(bed_inst)))

    cfg = RefinerMdlConfig()

    # ── Step 1: real presignals from the extended single BED pass ──────────
    # Build minimal Phase-0-style records (id, R, sequence, actual_length).
    records = []
    for fk, info in fams.items():
        rnum = int(fk[2:])
        records.append({
            'id': 'mdl_R%d' % rnum, 'R': rnum,
            'sequence': info['seq'], 'actual_length': info['length'],
            'length': info['length'], 'copies': info['copies'],
        })
    tmp = tempfile.mkdtemp(prefix='n1_int_')
    idx = os.path.join(tmp, 'instances_index.tsv')
    pre = os.path.join(tmp, 'completeness_presignals.tsv')
    phase0_triage.write_instance_index(records, BED, idx,
                                       presignal_path=pre, config=cfg)
    presig = load_presignals(pre)
    print("Presignal rows: %d (one per record)" % len(presig))

    # ── Step 2: assess_completeness on ALL families (zero genome alignment) ─
    verdicts = {}
    n_flag = 0
    sig_counts = {'bed_divergence_spread': 0, 'copy_count_vs_length_class': 0,
                  'under_instanced': 0}
    for rec in records:
        fk = "R=%d" % rec['R']
        insts = bed_inst.get(fk, [])
        v = assess_completeness(rec, insts, cfg)
        verdicts[fk] = v
        if v.incomplete_recall_eligible:
            n_flag += 1
        for k in sig_counts:
            if v.flags.get(k):
                sig_counts[k] += 1
    n = len(records)
    print("\n=== Step 2: completeness flags over ALL %d families ===" % n)
    print("flagged incomplete (recall-eligible): %d (%.1f%%)"
          % (n_flag, 100.0 * n_flag / n))
    print("per-signal fired counts:")
    for k, c in sig_counts.items():
        print("  %-28s %5d (%.1f%%)" % (k, c, 100.0 * c / n))

    # ── Step 3: sample flag vs blastn divergent-tail ground truth ──────────
    # R=150 is the known divergent-tail case; R=8/R=113/R=847 are homogeneous
    # low-spread controls picked from the real data.
    sample = ["R=150", "R=8", "R=113", "R=847", "R=2", "R=60"]
    print("\n=== Step 3: sample flag vs blastn ground-truth (real) ===")
    hdr = ("%-8s %6s %8s %8s | %7s %9s | %-8s %s"
           % ("R", "BED_n", "div_sprd", "len",
              "bl>=70", "bl70-90", "FLAG", "fired_signals"))
    print(hdr)
    print("-" * len(hdr))
    rows = []
    for fk in sample:
        if fk not in fams:
            continue
        v = verdicts[fk]
        insts = bed_inst.get(fk, [])
        sp = v.bed_divergence_spread
        h70, div = blastn_tail(fams[fk]['seq'], DB, fk)
        # ground truth: divergent tail iff blastn>=70 hits exceed BED count
        # AND there is a divergent (70-90%) component.
        bedn = len(insts)
        has_tail = (h70 > bedn) and (div > 0)
        fired = ",".join(k for k, val in v.flags.items() if val) or "-"
        print("%-8s %6d %8s %8d | %7d %9d | %-8s %s"
              % (fk, bedn,
                 ("%.3f" % sp) if sp is not None else "NA",
                 fams[fk]['length'], h70, div,
                 "YES" if v.incomplete_recall_eligible else "no", fired))
        rows.append((fk, v.incomplete_recall_eligible, has_tail))

    print("\n=== concordance (flag vs blastn-tail ground truth) ===")
    for fk, flag, tail in rows:
        verdict = "CONCORDANT" if flag == tail else "DISCORDANT"
        print("  %-8s flag=%-3s  blastn_tail=%-3s  -> %s"
              % (fk, "YES" if flag else "no", "YES" if tail else "no", verdict))


if __name__ == '__main__':
    main()
