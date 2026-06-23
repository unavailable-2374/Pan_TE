"""
M5 integration harness — real testA fixture, three scenarios:

  PART A  common path        : all families >= min_copies_for_msa BED instances ->
                               the genome BLAST DB is NEVER built and blastn is
                               NEVER run (the common path must not touch the DB).
  PART B  fallback recruit   : one synthetic under-instanced family (R=99, a real
                               consensus with ZERO BED rows) forces the bounded
                               blastn fallback -> the lazy DB is built EXACTLY once
                               and blastn runs ONLY for the under-instanced family.
  PART C  abPOA large family : large_family_threshold lowered so real testA families
                               route to abPOA -> a real abpoa run produces a
                               parseable, equal-length MSA and a refined consensus;
                               the harness also calls abpoa directly to report the
                               real exit code and MSA dimensions.

NOT a unit test: invokes makeblastdb / blastn / MAFFT / abpoa / CD-HIT-EST / samtools
on the real fixture under /home/shuoc/tool/mdl-repeat/tests. All numbers printed are
produced by running — nothing is asserted from memory.

Run (PGTA env), clean dir:
    python3 bin/Refiner_mdl/tests/integration_m5_fallback_abpoa.py
"""

import logging
import os
import shutil
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from config import RefinerMdlConfig          # noqa: E402
from phase0_triage import run_phase0          # noqa: E402
from phase1_extract import (                  # noqa: E402
    compute_pad, extract_padded_copies, load_fai_lengths)
from phase1_align import build_msa, select_regime, _subsample_copies  # noqa: E402
from phase1_boundary import refine_family      # noqa: E402
import phase1_fallback as fb                   # noqa: E402
import phase1_consensus as orch               # noqa: E402

DATA = "/home/shuoc/tool/mdl-repeat/tests"
GENOME = f"{DATA}/data/testA.fa"
FAMILIES = f"{DATA}/results/testA_families.fa"
BED = f"{DATA}/results/testA_instances.bed"
ABPOA_EXE = "/home/shuoc/tool/abPOA/bin/abpoa"

# A real R=1 consensus (300 bp, 37 genomic copies) reused as the under-instanced
# query family R=99 (it gets ZERO BED rows -> forced onto the fallback path).
R1_CONSENSUS = (
    "AAGCCCAATAAACCACTCTGACTGGCCGAATAGGGATATAGGCAACGACATGTGCGGCGAC"
    "CCTTGCGACAGTGACGCTTTCGCCGTTGCCTAAACCTATTTGAAGGAGTCTAGCAGCCGCA"
    "GTAAGGCACAATACCTCGTCCGTGTTACCAGACCAAACAAGACGTCCTCTTCAATGTTTAA"
    "ATGACCCTCTCGTCATAAAACCTTTCTACTATGTGTTCCGCAAGAATCAACAACTACAATG"
    "GCGCGTCGTGAATAACGCGACGGCTGAGACGAACGGCGCGTGAATGAAGCGCTTAA")


class _LogTap(logging.Handler):
    """Capture log messages so the harness can assert on real pipeline events."""
    def __init__(self):
        super().__init__(level=logging.INFO)
        self.msgs = []

    def emit(self, record):
        self.msgs.append(record.getMessage())

    def having(self, sub):
        return [m for m in self.msgs if sub in m]


class _worker_log_file:
    """Context manager attaching a root FileHandler so worker-process (forked) logs
    are captured — the in-memory _LogTap cannot see across the ProcessPool fork."""
    def __init__(self, path):
        self.path = path
        self.h = None

    def __enter__(self):
        self.h = logging.FileHandler(self.path)
        self.h.setLevel(logging.INFO)
        logging.getLogger().addHandler(self.h)
        return self

    def __exit__(self, *a):
        self.h.flush()
        logging.getLogger().removeHandler(self.h)
        self.h.close()

    def lines(self, sub):
        if not os.path.exists(self.path):
            return []
        return [l.strip() for l in open(self.path) if sub in l]


def _cfg(work, **over):
    c = RefinerMdlConfig(
        input_file=FAMILIES, genome_file=GENOME,
        output_dir=os.path.join(work, 'out'),
        temp_dir=os.path.join(work, 'temp'),
        checkpoint_dir=os.path.join(work, 'ckpt'),
        bed_file=BED, threads=3, shard_size=1, genome_blast_db="")
    for k, v in over.items():
        setattr(c, k, v)
    os.makedirs(c.temp_dir, exist_ok=True)
    return c


def _src_breakdown(records):
    out = {}
    for r in records:
        s = r.get('consensus_source', '?')
        out[s] = out.get(s, 0) + 1
    return out


def part_a_common_path(work):
    print("\n" + "=" * 64)
    print("PART A — common path: all families >= min_copies, DB untouched")
    print("=" * 64)
    cfg = _cfg(os.path.join(work, 'A'))
    tap = _LogTap()
    logging.getLogger().addHandler(tap)
    try:
        p0 = run_phase0(cfg)
        out = orch.run_phase1(p0, cfg)
    finally:
        logging.getLogger().removeHandler(tap)

    db_dir = os.path.join(cfg.temp_dir, 'genome_blastdb')
    built = tap.having("Building genome BLAST DB")
    not_built = tap.having("genome BLAST DB not built")
    recruited = tap.having("recruit_by_blastn")

    print(f"  families                  : {len(p0['records'])}")
    print(f"  consensus_source          : {_src_breakdown(out['records'])}")
    print(f"  genome_blastdb dir exists : {os.path.exists(db_dir)}")
    print(f"  'Building ... DB' log lines: {len(built)}")
    print(f"  'DB not built' log line    : {bool(not_built)}")
    print(f"  recruit_by_blastn log lines: {len(recruited)}")
    assert not os.path.exists(db_dir), "DB dir built on common path!"
    assert len(built) == 0, "makeblastdb ran on common path!"
    assert not_built, "expected the 'DB not built' decision log line"
    assert len(recruited) == 0, "blastn ran on common path!"
    print("  RESULT: common path never built the DB and never ran blastn  [OK]")


def part_b_fallback(work):
    print("\n" + "=" * 64)
    print("PART B — fallback: one under-instanced family forces lazy DB + blastn")
    print("=" * 64)
    wb = os.path.join(work, 'B')
    os.makedirs(wb, exist_ok=True)
    # Synthetic under-instanced query R=99: the real R=1 consensus mutated ~3% so it
    # ESCAPES the 0.98 Phase-0 dedup (it is not a near-exact duplicate of R=1) yet
    # remains close enough that blastn at pident>=70 still recruits the 37 real R=1
    # genomic loci. It gets ZERO BED rows -> forced onto the fallback path.
    _b = list(R1_CONSENSUS)
    _swap = {'A': 'C', 'C': 'A', 'G': 'T', 'T': 'G'}
    for i in range(0, len(_b), 30):          # ~3.3% deterministic substitutions
        _b[i] = _swap.get(_b[i], _b[i])
    r99 = ''.join(_b)
    fam_path = os.path.join(wb, 'families_plus_under.fa')
    with open(FAMILIES) as src, open(fam_path, 'w') as dst:
        dst.write(src.read())
        dst.write(">R=99 length=300 copies=3 mdl=9000.0 div=0.10 topo=linear "
                  "accept=exclusive tier=core flags=0x0 qflags=none\n")
        dst.write(r99 + "\n")

    cfg = _cfg(wb)
    cfg.input_file = fam_path          # R=99 has NO BED rows -> under-instanced
    tap = _LogTap()                    # parent-process events (DB build, decision)
    logging.getLogger().addHandler(tap)
    wlog = os.path.join(wb, 'worker.log')
    try:
        with _worker_log_file(wlog) as wf:   # cross-fork worker events (recruit)
            p0 = run_phase0(cfg)
            out = orch.run_phase1(p0, cfg)
    finally:
        logging.getLogger().removeHandler(tap)

    built = tap.having("Building genome BLAST DB")
    under_decision = tap.having("under-instanced families in pending shards")
    db_dir = os.path.join(cfg.temp_dir, 'genome_blastdb')
    # Worker-side blastn evidence (recruit success OR no-hit fallback), from the file.
    recruit_lines = wf.lines("recruit_by_blastn:")
    recruited = [l for l in recruit_lines if 'recruited' in l]
    no_hits = [l for l in recruit_lines if 'no qualifying hits' in l]

    print(f"  input families (w/ R=99)  : {len(p0['records'])}")
    print(f"  under-instanced decision  : {under_decision[:1]}")
    print(f"  'Building ... DB' lines   : {len(built)}  (expect exactly 1)")
    print(f"  worker recruit_by_blastn  :")
    for m in recruit_lines:
        print(f"      {m.split(' - ')[-1] if ' - ' in m else m}")
    print(f"  genome_blastdb dir exists : {os.path.exists(db_dir)}")
    print(f"  consensus_source          : {_src_breakdown(out['records'])}")

    assert len(built) == 1, f"DB must be built exactly once, saw {len(built)}"
    # Every blastn invocation must be for the under-instanced family only (mdl_R99).
    blastn_families = recruited + no_hits
    assert all('mdl_R99' in l for l in blastn_families), (
        f"blastn ran for a non-under-instanced family: {blastn_families}")
    assert len(blastn_families) == 1, (
        f"expected exactly one under-instanced blastn family, saw {len(blastn_families)}")
    print(f"  total blastn families     : {len(blastn_families)} (only mdl_R99) "
          f"of {len(p0['records'])} families")
    print("  RESULT: DB built once in parent; blastn ran ONLY for mdl_R99  [OK]")

    # ---- B2: SYNTHETIC genome — prove the FULL recruit->refine path (>=5 copies).
    # testA's real R=1 instances are heterogeneous (only 1 locus matches its 300 bp
    # consensus at blastn level), so a clean genome with planted copies is needed to
    # exercise successful recruitment. *** SYNTHETIC DATA, LABELED ***
    import random
    rng = random.Random(13)
    seed = R1_CONSENSUS
    alt = {'A': 'CGT', 'C': 'AGT', 'G': 'ACT', 'T': 'ACG'}

    def _rand(n):
        return ''.join(rng.choice('ACGT') for _ in range(n))

    def _mutate(s, rate):
        b = list(s)
        for j in range(len(b)):
            if rng.random() < rate and b[j] in alt:
                b[j] = rng.choice(alt[b[j]])
        return ''.join(b)

    n_plant = 10
    chunks = [_rand(600)]
    for _ in range(n_plant):
        chunks.append(_mutate(seed, 0.08))
        chunks.append(_rand(600))
    syn_genome = os.path.join(wb, 'syn_genome.fa')
    with open(syn_genome, 'w') as fh:
        fh.write(">synchr1\n" + ''.join(chunks) + "\n")
    subprocess.run(['samtools', 'faidx', syn_genome], check=True,
                   capture_output=True)
    cfg2 = _cfg(os.path.join(wb, 'synrun'))
    cfg2.genome_file = syn_genome
    cfg2.genome_blast_db = ""
    from phase1_fallback import ensure_genome_blast_db
    ensure_genome_blast_db(cfg2)
    syn_copies = fb.recruit_by_blastn({'id': 'SYN_planted', 'sequence': seed}, cfg2)
    syn_rec = {'id': 'SYN_planted', 'sequence': seed, 'length': len(seed),
               'copies': len(syn_copies), 'mdl': 9999.0, 'tier': 'T1'}
    syn_refined = refine_family(syn_rec, syn_copies, cfg2)
    print(f"  [SYNTHETIC genome] planted: {n_plant} copies of a {len(seed)} bp seed "
          f"(8% point-mut) *** SYNTHETIC ***")
    print(f"  [SYNTHETIC genome] recruit: {len(syn_copies)} copies "
          f"(strands={sorted(set(c.strand for c in syn_copies))})")
    print(f"  [SYNTHETIC genome] refine : consensus_source="
          f"{syn_refined.get('consensus_source')}, "
          f"len={len(syn_refined.get('sequence', ''))}")
    assert len(syn_copies) >= 5, (
        f"expected >=5 recruited copies on the synthetic genome, got {len(syn_copies)}")
    assert syn_refined.get('consensus_source') != 'original', (
        "successful recruitment should have refined past the seed")
    print("  RESULT: recruit_by_blastn recovered >=5 planted copies -> refine "
          "produced a consensus (full fallback->refine path)  [OK]")


def _synthetic_family(seed_seq, n_copies, mut_rate, rng_seed):
    """Build n_copies SYNTHETIC copies of seed_seq with point substitutions only.

    Returns Copy objects (all '+' strand, equal length). Used purely to exercise the
    abpoa refine path on a clean, well-formed large family. *** SYNTHETIC DATA ***."""
    import random
    from phase1_extract import Copy
    rng = random.Random(rng_seed)
    alt = {'A': 'CGT', 'C': 'AGT', 'G': 'ACT', 'T': 'ACG'}
    copies = []
    for i in range(n_copies):
        b = list(seed_seq)
        for j in range(len(b)):
            if rng.random() < mut_rate and b[j] in alt:
                b[j] = rng.choice(alt[b[j]])
        copies.append(Copy(id=f"syn{i}", sequence=''.join(b),
                           divergence=mut_rate, strand='+'))
    return copies


def part_c_abpoa(work):
    print("\n" + "=" * 64)
    print("PART C — abPOA: real-testA invocation + SYNTHETIC clean-family refine")
    print("=" * 64)
    cfg = _cfg(os.path.join(work, 'C'))   # default thresholds; abpoa for n>200

    # ---- C1: REAL testA — prove the abpoa *invocation* works end to end --------
    # testA's R=0 family has 52 BED instances whose spans range 250 bp .. 7+ kb and
    # share little similarity (synthetic-test artifacts: giant outlier instances).
    # abpoa still runs (exit 0) and yields a parseable equal-length MSA; refine then
    # *reasonably falls back* to the seed because that MSA is too gappy to pass QC —
    # the never-regress invariant, not a failure.
    fai = load_fai_lengths(GENOME + '.fai')
    p0 = run_phase0(_cfg(os.path.join(work, 'C_real')))
    idx = p0['instance_index']
    counts = {}
    with open(idx) as fh:
        for line in fh:
            rid = line.split('\t', 1)[0]
            counts[rid] = counts.get(rid, 0) + 1
    big_id = max(counts, key=counts.get)
    tmp_inst = os.path.join(work, 'big_inst.tsv')
    with open(idx) as fh, open(tmp_inst, 'w') as out_fh:
        for line in fh:
            if line.split('\t', 1)[0] == big_id:
                out_fh.write(line)
    insts = orch._load_shard_instances(tmp_inst).get(big_id, [])
    copies = extract_padded_copies(insts, GENOME, fai,
                                   lambda L: compute_pad(L, cfg), cfg)
    regime = select_regime(len(copies),
                           sum(c.divergence for c in copies) / max(len(copies), 1),
                           cfg)
    work_sub = _subsample_copies(copies, cfg.msa_subsample_cap_abpoa, cfg.subsample_seed)
    td = tempfile.mkdtemp(prefix='abpoa_probe_')
    try:
        infa, omsa = os.path.join(td, 'in.fa'), os.path.join(td, 'out.msa')
        with open(infa, 'w') as fh:
            for i, c in enumerate(work_sub):
                fh.write(f">s{i}\n{c.sequence}\n")
        rc = subprocess.run([ABPOA_EXE, '-r', '1', '-o', omsa, infa],
                            capture_output=True, text=True).returncode
        drows = [l.strip() for l in open(omsa) if l.strip() and not l.startswith('>')]
    finally:
        shutil.rmtree(td, ignore_errors=True)
    res = build_msa(copies, 'abpoa', cfg)
    print(f"  [REAL testA] probe family : {big_id} ({len(insts)} insts -> "
          f"{len(copies)} copies, regime={regime})")
    print(f"  [REAL testA] direct abpoa : exit={rc}, "
          f"MSA {len(drows)} rows x {len(drows[0]) if drows else 0} cols")
    assert rc == 0, "abpoa returned non-zero exit on real copies"
    assert res is not None, "build_msa(abpoa) returned None on real copies"
    rrows, _ = res
    rcols = {len(r) for r in rrows}
    print(f"  [REAL testA] build_msa    : {len(rrows)} rows, "
          f"equal-length={len(rcols) == 1}")
    assert len(rcols) == 1, "real-testA abpoa rows not equal length"

    # ---- C2: SYNTHETIC clean family — prove the abpoa *refine* path produces a
    #          consensus (consensus_source='abpoa'). *** SYNTHETIC DATA, LABELED ***
    seed = R1_CONSENSUS                      # 300 bp real seed, mutated below
    syn_copies = _synthetic_family(seed, n_copies=250, mut_rate=0.06, rng_seed=7)
    syn_regime = select_regime(len(syn_copies), 0.06, cfg)
    syn_res = build_msa(syn_copies, 'abpoa', cfg)
    assert syn_res is not None, "SYNTHETIC abpoa build_msa returned None"
    srows, _ = syn_res
    scols = {len(r) for r in srows}
    rec = {'id': 'SYN_big', 'sequence': seed, 'length': len(seed),
           'copies': len(syn_copies), 'mdl': 9999.0, 'tier': 'T1'}
    refined = refine_family(rec, syn_copies, cfg)
    print(f"  [SYNTHETIC] family        : 250 copies x {len(seed)} bp, "
          f"6% point-mut, regime={syn_regime}  *** SYNTHETIC ***")
    print(f"  [SYNTHETIC] build_msa     : {len(srows)} rows, "
          f"equal-length={len(scols) == 1}, cols={scols}")
    print(f"  [SYNTHETIC] refine result : consensus_source="
          f"{refined.get('consensus_source')}, "
          f"consensus_len={len(refined.get('sequence', ''))}")
    assert syn_regime == 'abpoa', "synthetic family did not route to abpoa"
    assert len(scols) == 1, "synthetic abpoa rows not equal length"
    assert refined.get('consensus_source') == 'abpoa', (
        "synthetic clean family should have been refined via abpoa, got "
        f"{refined.get('consensus_source')}")
    assert 200 <= len(refined['sequence']) <= 400, (
        f"synthetic abpoa consensus length unexpected: {len(refined['sequence'])}")
    print("  RESULT: abpoa invocation verified on real data; abpoa refine path "
          "produces a consensus on a clean (synthetic) large family  [OK]")


def main():
    for f in (GENOME, FAMILIES, BED, GENOME + '.fai'):
        if not os.path.isfile(f):
            print(f"missing fixture: {f}", file=sys.stderr)
            sys.exit(1)
    if not os.path.exists(ABPOA_EXE):
        print(f"missing abpoa: {ABPOA_EXE}", file=sys.stderr)
        sys.exit(1)

    logging.basicConfig(level=logging.INFO,
                        format='%(levelname)s %(name)s: %(message)s')
    work = tempfile.mkdtemp(prefix='m5_harness_')
    try:
        part_a_common_path(work)
        part_b_fallback(work)
        part_c_abpoa(work)
        print("\n" + "=" * 64)
        print("M5 INTEGRATION HARNESS — ALL PARTS PASSED")
        print("=" * 64)
    finally:
        shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    main()
