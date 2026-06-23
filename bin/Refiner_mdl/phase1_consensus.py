"""
Phase 1 (rewrite) — Module M4: orchestrator (shard plan, BED routing, parallel
per-shard refine, concat, post-refine merge).

This replaces the old monolithic whole-genome-BLAST `phase1_consensus.py`. The
per-family refinement logic now lives in M1-M3 (`phase1_extract` / `phase1_align` /
`phase1_boundary` / `phase1_chimera`); this module only *schedules* it:

  1. make_shards   — partition records into shards of `shard_size`, distributing
                     heavy (len*copies) families across shards for load balance.
  2. route_instances — stream the Phase 0 instance index ONCE (it is pre-sorted by
                     record_id) and append each family's BED rows to its owner
                     shard's instance file. Streaming group-by, no in-RAM index.
  3. _process_shard — one worker per shard (ProcessPoolExecutor). Per family:
                     extract padded copies (scaffolds: rebuild spanning loci) →
                     CLUSTER into subfamilies (N2) → per qualifying cluster
                     refine_family → (chimera-eligible?) detect_chimera. A
                     homogeneous family yields ONE cluster and reproduces today's
                     single-consensus output verbatim (original id); a divergent
                     family yields several `mdl_R<N>_sf<k>` subfamilies. Every
                     family is wrapped in try/except: on ANY error the original
                     mdl-repeat seed is kept — never a faked success. Shards whose
                     `.done` marker exists are skipped (resume).
  4. concat + merge_post_refine — collect refined shards, then a single CD-HIT-EST
                     at `post_refine_merge_identity` (the looser merge DEFERRED from
                     Phase 0) collapses boundary-convergent duplicates.

Data-integrity invariant (logged at every stage): families are never silently
dropped. Output count == input count + fragments_added − post-merge merged, where
fragments_added now counts BOTH structural chimera fragments AND subfamily splits (N2:
a divergent family emitting k>1 subfamily consensi adds k-1, §8.2). Each term is
tracked from an independent source so the reconciliation is a real cross-check, not a
tautology — a silently-dropped subfamily member still trips the assertion.

See REFINE_IMPLEMENTATION_PLAN.md §1.6, §1.7, §8, §11, §12.M4.
"""

import json
import logging
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from statistics import mean
from types import SimpleNamespace
from typing import Dict, List, Optional

import numpy as np

from phase0_triage import _run_cdhit
from phase1_align import select_regime  # noqa: F401 (kept for symmetry / re-export)
from phase1_boundary import refine_family
from phase1_chimera import detect_chimera, strip_aln
from phase1_cluster import (
    cluster_copies,
    qualifying_clusters,
    prefilter_copies_for_clustering,
    subfamily_id,
    subfamily_lineage,
)
from phase1_completeness import recall_eligible
from phase1_fp_prune import pseudo_family_verdict, PseudoFamilyVerdict
from phase1_extract import (
    Instance,
    compute_pad,
    extract_padded_copies,
    load_fai_lengths,
)
from phase1_fallback import ensure_genome_blast_db, recruit_by_blastn
from phase1_seed_chimera import (
    _rmblastn_hsps, hsps_to_instances_loci, detect_seed_breakpoints, split_seed,
)
from phase1_copy_trim import trim_copies_structural


def _rmblastn_units(rec: Dict, cfg, fai) -> "List[Tuple[Dict, list]]":
    """#1 + N7: recruit the family's copies by rmblastn (instead of mdl BED) and, if the
    seed is a nested chimera, split it. Returns a list of (record, copies) units — one
    unit normally, >=2 when the seed was de-nested. Returns [] (caller falls back to BED)
    when rmblastn yields nothing.

    Each instance is extracted with padding via the SAME extract_padded_copies as the BED
    path, so the copies are indistinguishable to clustering/MSA downstream.
    """
    gap = max(1, len(rec.get('sequence', '')))
    qlen = len(rec.get('sequence', ''))
    min_cov = getattr(cfg, 'copy_extract_min_cov', 0.5)
    hsps = _rmblastn_hsps(rec['sequence'], cfg)
    if not hsps:
        return []
    instances, loci = hsps_to_instances_loci(hsps, gap, qlen)
    if not instances:
        return []

    def _cov_filtered(insts, lc):
        """Keep only instances whose HSPs cover >= min_cov of the consensus — drops
        partial gene-paralog matches from the copy set fed to the MSA (N7 keeps all)."""
        kept = [ins for ins, l in zip(insts, lc) if l.get('cov', 1.0) >= min_cov]
        return kept if kept else insts        # never starve to empty -> keep all

    def _copies(insts, lc):
        return extract_padded_copies(_cov_filtered(insts, lc), cfg.genome_file, fai,
                                     lambda L: compute_pad(L, cfg), cfg)

    # ── N7 seed-chimera resolution (off by default; uses ALL loci, unfiltered) ────
    if getattr(cfg, 'enable_seed_chimera', False):
        bps = detect_seed_breakpoints(loci, qlen, cfg)
        if bps:
            comps = split_seed(rec, bps, cfg)
            if comps:
                # Re-recruit each component independently; post-split improvement gate.
                units, ok = [], True
                for c in comps:
                    chsps = _rmblastn_hsps(c['sequence'], cfg)
                    cinst, clc = hsps_to_instances_loci(chsps, max(1, len(c['sequence'])),
                                                        len(c['sequence']))
                    if (getattr(cfg, 'chimera_seed_require_standalone', True)
                            and len(cinst) < cfg.min_copies_for_msa):
                        ok = False; break          # a component can't stand alone -> revert
                    units.append((c, _copies(cinst, clc)))
                if ok:
                    logger.info("shard: seed %s de-nested into %d components at %s",
                                rec.get('id', '?'), len(comps), bps)
                    return units
    return [(rec, _copies(instances, loci))]

logger = logging.getLogger(__name__)


# ═══════════════════════════════════════════════════════════════════════
# JSONL helpers (records must be _aln-free and numpy-free on disk)
# ═══════════════════════════════════════════════════════════════════════

def _json_default(o):
    """Coerce numpy scalars (entropy / dust_score may be np.float64) to native."""
    if isinstance(o, np.generic):
        return o.item()
    return str(o)


def _write_jsonl(path: str, records: List[Dict]) -> None:
    with open(path, 'w') as fh:
        for rec in records:
            fh.write(json.dumps(rec, default=_json_default) + '\n')


def _read_jsonl(path: str) -> List[Dict]:
    out: List[Dict] = []
    if not os.path.exists(path):
        return out
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line:
                out.append(json.loads(line))
    return out


def _shard_paths(shard_dir: str, i: int) -> Dict[str, str]:
    return {
        'records': os.path.join(shard_dir, f'shard_{i}.records.jsonl'),
        'instances': os.path.join(shard_dir, f'shard_{i}.instances.tsv'),
        'refined': os.path.join(shard_dir, f'shard_{i}.refined.jsonl'),
        'pruned': os.path.join(shard_dir, f'shard_{i}.pruned.jsonl'),
        'consensus': os.path.join(shard_dir, f'shard_{i}.consensus.fa'),
        'stats': os.path.join(shard_dir, f'shard_{i}.stats.json'),
        'done': os.path.join(shard_dir, f'shard_{i}.done'),
    }


# ═══════════════════════════════════════════════════════════════════════
# 1. Shard plan
# ═══════════════════════════════════════════════════════════════════════

def make_shards(records: List[Dict], shard_size: int) -> List[List[Dict]]:
    """Partition records into shards of ~shard_size, spreading heavy families.

    Records are sorted by len*copies descending and striped round-robin across
    n_shards = ceil(len/shard_size) bins, so the most expensive families do not all
    land in one shard (load balance). Count is conserved exactly: every record lands
    in exactly one shard."""
    n = len(records)
    if n == 0:
        return []
    shard_size = max(1, shard_size)
    n_shards = (n + shard_size - 1) // shard_size
    ordered = sorted(
        records,
        key=lambda r: len(r.get('sequence', '')) * r.get('copies', 1),
        reverse=True,
    )
    shards: List[List[Dict]] = [[] for _ in range(n_shards)]
    for idx, rec in enumerate(ordered):
        shards[idx % n_shards].append(rec)
    return shards


# ═══════════════════════════════════════════════════════════════════════
# 2. Instance routing (streaming group-by over the sorted index)
# ═══════════════════════════════════════════════════════════════════════

def route_instances(index_path: str, record_to_shard: Dict[str, int],
                    shard_dir: str) -> Dict[int, int]:
    """Stream the record_id-sorted instance index once → per-shard instance files.

    The index is pre-sorted by record_id (Phase 0 `sort -k1,1`), so all rows of a
    given family are contiguous: we accumulate a family's rows, then on record_id
    change append them to the owner shard's instances file with a single open().
    This bounds open() calls to the number of distinct families and memory to one
    family's rows — no full in-RAM index (10-30 Gb tractable).

    Returns {shard_idx: n_rows_routed} for logging. Rows whose record_id is not in
    record_to_shard (e.g. a family dropped between index build and sharding) are
    counted and skipped — never misrouted."""
    routed_counts: Dict[int, int] = {}
    n_orphan = 0

    def _flush(rec_id: Optional[str], buf: List[str]) -> None:
        nonlocal n_orphan
        if rec_id is None or not buf:
            return
        shard_idx = record_to_shard.get(rec_id)
        if shard_idx is None:
            n_orphan += len(buf)
            return
        path = _shard_paths(shard_dir, shard_idx)['instances']
        with open(path, 'a') as fh:
            fh.writelines(buf)
        routed_counts[shard_idx] = routed_counts.get(shard_idx, 0) + len(buf)

    cur_id: Optional[str] = None
    buf: List[str] = []
    with open(index_path) as fh:
        for line in fh:
            if not line.strip():
                continue
            rec_id = line.split('\t', 1)[0]
            if rec_id != cur_id:
                _flush(cur_id, buf)
                cur_id = rec_id
                buf = []
            buf.append(line if line.endswith('\n') else line + '\n')
    _flush(cur_id, buf)

    if n_orphan:
        logger.info("route_instances: %d index rows had no owning shard (skipped)",
                    n_orphan)
    return routed_counts


def _count_under_instanced(index_path: Optional[str],
                           record_to_shard: Dict[str, int],
                           pending_shards: set, min_copies: int) -> int:
    """Count families (in pending shards) with fewer than `min_copies` BED rows.

    Streams the record_id-sorted index ONCE (group-by, no RAM index) to get a per
    record row count, then counts records routed to a pending shard whose count is
    below min_copies — these are exactly the families that would trigger the bounded
    blastn fallback, so a non-zero count is the lazy-DB build trigger. Families absent
    from the index (count 0) are under-instanced too. Raw-row counts over-estimate
    scaffold copies slightly (member rows collapse to fewer spanning loci); the only
    cost of an over-count here is building the DB when a borderline scaffold would not
    have used it — never a missed build for a normal under-instanced family."""
    pending_ids = {rid for rid, s in record_to_shard.items()
                   if s in pending_shards}
    if not pending_ids:
        return 0
    # No index at all -> every pending family has 0 instances (all under-instanced).
    if not index_path or not os.path.exists(index_path):
        return len(pending_ids)

    counts: Dict[str, int] = {}
    cur_id: Optional[str] = None
    cur_n = 0
    with open(index_path) as fh:
        for line in fh:
            if not line.strip():
                continue
            rid = line.split('\t', 1)[0]
            if rid != cur_id:
                if cur_id is not None:
                    counts[cur_id] = cur_n
                cur_id = rid
                cur_n = 0
            cur_n += 1
    if cur_id is not None:
        counts[cur_id] = cur_n

    return sum(1 for rid in pending_ids if counts.get(rid, 0) < min_copies)


def _load_shard_instances(path: str) -> Dict[str, List[Instance]]:
    """Read a shard instances TSV → {record_id: [Instance, ...]}."""
    by_rec: Dict[str, List[Instance]] = {}
    if not os.path.exists(path):
        return by_rec
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 7:
                continue
            rec_id = parts[0]
            try:
                inst = Instance(chrom=parts[2], start=int(parts[3]),
                                end=int(parts[4]), strand=parts[5],
                                divergence=float(parts[6]))
            except ValueError:
                continue
            by_rec.setdefault(rec_id, []).append(inst)
    return by_rec


# ═══════════════════════════════════════════════════════════════════════
# 1.7 carried-over behaviors: chimera eligibility + scaffold instances
# ═══════════════════════════════════════════════════════════════════════

def chimera_eligible(rec: Dict, config) -> bool:
    """Tier / length / topology gate deciding *whether* to call detect_chimera.

    Preserves the old policy: topology=='complex' always eligible (overrides the
    length guard on any tier); otherwise T1 needs len>=300, T2 needs len>=150, and
    T3 is eligible only via the topology==complex path. Length is read from the
    record's current (refined) sequence."""
    if rec.get('topology') == 'complex':
        return True
    tier = rec.get('tier', 'T2')
    slen = len(rec.get('sequence', ''))
    if tier == 'T1':
        return slen >= 300
    if tier == 'T2':
        return slen >= 150
    return False  # T3 (non-complex)


def build_scaffold_instances(rec: Dict, routed_instances: List[Instance],
                             config) -> List[Instance]:
    """Reconstruct spanning loci for a scaffold record from its member instances.

    A scaffold record (Phase 0 chain assembly) owns the BED rows of every member
    family. On each chromosome+strand, member instances that sit within
    `fragment_gap` of each other in coordinate order belong to one assembled element
    copy: they are merged into a single spanning Instance
    (chrom, first_start, last_end, strand, mean_member_div). If at least
    `min_copies_for_msa` spanning loci are found the scaffold is refined like a
    normal family; otherwise an empty list is returned and refine_family keeps the
    concatenated seed (consensus_source='original')."""
    if not routed_instances:
        return []
    ordered = sorted(routed_instances, key=lambda x: (x.chrom, x.strand, x.start))
    spanning: List[Instance] = []
    cur_chrom = cur_strand = None
    cur_start = cur_end = None
    cur_divs: List[float] = []

    def _emit():
        if cur_chrom is not None and cur_start is not None:
            spanning.append(Instance(
                chrom=cur_chrom, start=cur_start, end=cur_end,
                strand=cur_strand,
                divergence=mean(cur_divs) if cur_divs else 0.0))

    for it in ordered:
        if (it.chrom == cur_chrom and it.strand == cur_strand
                and cur_end is not None
                and it.start - cur_end <= config.fragment_gap):
            cur_end = max(cur_end, it.end)
            cur_divs.append(it.divergence)
        else:
            _emit()
            cur_chrom, cur_strand = it.chrom, it.strand
            cur_start, cur_end = it.start, it.end
            cur_divs = [it.divergence]
    _emit()

    if len(spanning) < config.min_copies_for_msa:
        return []
    return spanning


# ═══════════════════════════════════════════════════════════════════════
# 2.5 Per-family body (N2): cluster → per-cluster M2 → per-subfamily M3
# ═══════════════════════════════════════════════════════════════════════

def _refine_one_family(rec: Dict, copies: List, cfg):
    """Refine ONE family into one-or-more subfamily-resolved consensus records.

    The N2 per-family body (REFINE_STRATEGY_DESIGN_v2.md §3.4, §3.6): the assembled
    copy set is ALWAYS clustered; every qualifying cluster is handed to M2's
    `refine_family` (the same extend-align-trim, invoked per cluster), then M3's
    chimera check runs per emitted subfamily consensus (unchanged from M4, just
    relocated after clustering). The number of consensi is emergent — one cluster for
    a homogeneous family (output byte-identical to today: original id, all copies, one
    `refine_family` call), k clusters for a divergent one.

    Naming / lineage (§3.5): a single-cluster family keeps its original `mdl_R<N>` id;
    a multi-cluster family suffixes `_sf<k>`. Each subfamily record carries:
      * id           — original or `_sf<k>`
      * copies       — the SUBFAMILY member count (not the parent's), so Phase 2's
                       per-subfamily copy gates act correctly (§7)
      * subfamily    — the biological lineage dict (provenance, never a FASTA header)
    A subfamily's `refine_family` may itself fall back to the seed (never-regress);
    that still emits exactly one record per cluster, so the count is well-defined.

    N5 (REFINE_STRATEGY_DESIGN_v2.md §4): the qualifying-cluster `audit` and the copy
    set drive the conservative pseudo-family discriminant. When the family forms NO
    qualifying cluster AND its copies are positively incoherent, the verdict's
    `is_pseudo` is True and this returns an EMPTY record list — the caller drops the
    family into the audit list (it is never silently lost; the drop is an explicit
    counted term). When N5 is disabled, or the family is coherent / low-copy / not
    recall-verified, `is_pseudo` is False and the family is refined exactly as N4 did.
    Building consensus is SKIPPED for a confirmed pseudo-family (no wasted M2 work).

    Returns (family_out, verdict) where `family_out` is the flat list of output
    records (>= 1 when kept, [] when dropped) and `verdict` is the PseudoFamilyVerdict
    (so the caller can record id + reason + metrics for a dropped family).
    """
    # ── Sweet-spot guard (cost↔value alignment) ──────────────────────────────
    # Subfamily MSA refinement of large / long-element families is intrinsically a
    # minutes-scale cost (MAFFT or abPOA on long divergent sequences) AND, measured,
    # those families fall back to the original seed anyway (R=1567: 927×8.5 kb → 3×
    # original, post-merged away) — pure wasted compute. Long elements are also
    # Look4LTRs' remit in the pipeline. So when the element is long or the total MSA
    # work is large, short-circuit to the mdl-repeat seed: no clustering, no MSA, no
    # subfamily. Value is unaffected (it lives in shorter, moderately-sized families);
    # cost is bounded so the full genome is tractable.
    _slen = len(rec.get('sequence', ''))
    _max_len = getattr(cfg, 'subfamily_max_consensus_len', 5000)
    _work_budget = getattr(cfg, 'subfamily_work_budget', 1_500_000)
    if _slen > _max_len or len(copies) * _slen > _work_budget:
        out = dict(rec)
        out['sequence'] = rec.get('sequence', '')
        out['actual_length'] = _slen
        out['consensus_source'] = 'original'
        return [out], PseudoFamilyVerdict(
            reason='oversized family: kept mdl-repeat seed, subfamily refinement skipped')

    # If there are not even enough copies to MSA, clustering is moot: one cluster of
    # whatever copies exist, refined (refine_family will itself keep the seed). This
    # path reproduces today's single-consensus behavior for sparse families exactly.
    #
    # Bound the O(copies²) clustering cost ONCE here so cluster_copies and
    # qualifying_clusters share one index space: drop length-outlier instances and
    # count-cap ultra-abundant families to a divergence-stratified representative
    # subsample (measured necessary on R=1567: 927 copies incl. a 671 kb mis-annotated
    # instance took 210 s before this cap).
    # Per-copy structural (TSD) boundary trim — cut shared host-gene flanks before the
    # MSA so the consensus does not over-extend into CDS. No-op unless enabled / TSD found.
    copies = trim_copies_structural(copies, cfg)
    copies = prefilter_copies_for_clustering(copies, cfg,
                                             ref_len=len(rec.get('sequence', '')))
    raw_clusters = cluster_copies(copies, cfg)
    quals, audit = qualifying_clusters(raw_clusters, copies, cfg)

    # ── N5 conservative FP prune (§4): decide BEFORE building any consensus so a
    #    confirmed pseudo-family costs no M2 work. Defaults to KEEP; needs positive
    #    incoherence evidence on top of zero qualifying clusters to flip is_pseudo. ──
    verdict = pseudo_family_verdict(rec, copies, audit, cfg)
    if verdict.is_pseudo:
        return [], verdict

    n_sf = len(quals)
    family_out: List[Dict] = []
    for k, cluster in enumerate(quals, start=1):
        sf_id = subfamily_id(rec['id'], k, n_sf)
        # Per-cluster record view: original record fields, but id + copies scoped to
        # this subfamily so M2/M3/Phase 2 see a coherent per-subfamily family.
        sf_rec = dict(rec)
        sf_rec['id'] = sf_id
        if n_sf > 1:
            sf_rec['copies'] = cluster.size
            sf_rec['subfamily'] = subfamily_lineage(rec, k, n_sf, cluster, cfg)

        refined = refine_family(sf_rec, cluster.members, cfg)

        if chimera_eligible(refined, cfg):
            frags = detect_chimera(refined, cfg)
        else:
            frags = [strip_aln(refined)]
        family_out.extend(frags)

    return family_out, verdict


# ═══════════════════════════════════════════════════════════════════════
# 3. Per-shard worker
# ═══════════════════════════════════════════════════════════════════════

def _process_shard(shard_idx: int, shard_dir: str, config_dict: Dict) -> Dict:
    """Worker: refine every family in one shard. Writes refined.jsonl + consensus.fa
    + stats.json, then touches the .done marker. Returns the shard stats dict.

    Per-family try/except guarantees a family is NEVER lost: on any exception the
    original mdl-repeat seed (a single record) is emitted. The worker-internal
    invariant n_out == n_in + extra (extra = sum of fragments beyond the first per
    family) is asserted before the shard is marked done."""
    cfg = SimpleNamespace(**config_dict)
    paths = _shard_paths(shard_dir, shard_idx)

    records = _read_jsonl(paths['records'])
    inst_by_rec = _load_shard_instances(paths['instances'])
    fai = load_fai_lengths(cfg.genome_file + '.fai')

    out_records: List[Dict] = []
    pruned_records: List[Dict] = []   # N5: dropped pseudo-families (id + reason + metrics)
    n_in = len(records)
    extra = 0          # sum over families of (n_fragments - 1)
    split_families = 0
    n_pruned = 0       # N5: families dropped as pseudo (explicit counted term, §4/§8.2)

    # N4 per-worker recall budget: the global cap (config.recall_global_family_budget)
    # is split evenly across the worker pool by the parent (cfg.recall_worker_budget);
    # an absent field (e.g. an old config) means "no recall budget here". Counts
    # families that actually TRIGGERED a Tier-2 recall in THIS worker.
    recall_budget = getattr(cfg, 'recall_worker_budget', 0)
    n_recalled = 0
    n_over_budget = 0

    for rec in records:
        rec_added = 0      # output records this input family contributed (across its units)
        rec_pruned = 0     # pseudo-family drops this input family contributed
        # ── Resolve the unit(s) to refine ─────────────────────────────────
        # rmblastn copy source (#1) + N7 seed de-nesting (>=1 unit), ELSE the legacy
        # BED + N4-recall path (one unit). A unit is (record, copies). On any failure
        # keep the seed.
        try:
            units = (_rmblastn_units(rec, cfg, fai)
                     if getattr(cfg, 'enable_rmblastn_copies', False) else [])
            if not units:
                if rec.get('is_scaffold'):
                    insts = build_scaffold_instances(
                        rec, inst_by_rec.get(rec['id'], []), cfg)
                else:
                    insts = inst_by_rec.get(rec['id'], [])
                # BED copies are the base set; N4 recall only ADDS the divergent tail.
                bed_copies = extract_padded_copies(
                    insts, cfg.genome_file, fai, lambda L: compute_pad(L, cfg), cfg)
                recall_v = (recall_eligible(rec, insts, cfg)
                            if cfg.enable_selective_recall else None)
                is_eligible = bool(recall_v and recall_v.incomplete_recall_eligible
                                   and getattr(cfg, 'genome_blast_db', ''))
                if is_eligible and n_recalled < recall_budget:
                    recruited = recruit_by_blastn(rec, cfg, recall_mode=True)
                    existing_ids = {c.id for c in bed_copies}
                    added = [c for c in recruited if c.id not in existing_ids]
                    copies = list(bed_copies) + added
                    rec['completeness_verified'] = True
                    rec['recalled_copies'] = len(added)
                    n_recalled += 1
                    if added:
                        logger.info("shard %d: family %s recalled %d divergent copies "
                                    "(BED had %d) -> %d total for clustering",
                                    shard_idx, rec.get('id', '?'), len(added),
                                    len(existing_ids), len(copies))
                elif is_eligible:
                    copies = bed_copies
                    rec['completeness_verified'] = False
                    n_over_budget += 1
                    logger.warning("shard %d: family %s recall-eligible but over global "
                                   "budget (%d/worker); keeping BED copies, "
                                   "completeness_verified=False",
                                   shard_idx, rec.get('id', '?'), recall_budget)
                elif (len(insts) < cfg.min_copies_for_msa
                        and getattr(cfg, 'genome_blast_db', '')
                        and cfg.enable_fallback_recruitment):
                    copies = recruit_by_blastn(rec, cfg)
                else:
                    copies = bed_copies
                units = [(rec, copies)]
        except Exception as e:  # noqa: BLE001 — never fake success; keep the seed
            logger.warning("shard %d: family %s copy resolution failed (%s); keeping "
                           "seed", shard_idx, rec.get('id', '?'), e)
            seed = dict(rec)
            seed.pop('_aln', None)
            seed.setdefault('consensus_source', 'original')
            out_records.append(seed)        # 1 input -> 1 output, extra contribution 0
            continue

        # ── Refine each unit; N5 prune handled per unit. Per-input-family accounting
        #    keeps the conservation invariant a real cross-check (seed-splits add units). ──
        for urec, ucopies in units:
            try:
                refine_ret = _refine_one_family(urec, ucopies, cfg)
                if isinstance(refine_ret, tuple):
                    family_out, verdict = refine_ret
                else:
                    family_out, verdict = refine_ret, None
            except Exception as e:  # noqa: BLE001
                logger.warning("shard %d: family %s refine failed (%s); keeping seed",
                               shard_idx, urec.get('id', '?'), e)
                seed = dict(urec)
                seed.pop('_aln', None)
                seed.setdefault('consensus_source', 'original')
                family_out = [seed]
                verdict = None

            if verdict is not None and verdict.is_pseudo:
                n_pruned += 1
                rec_pruned += 1
                pruned_records.append({
                    'id': urec.get('id', '?'),
                    'reason': verdict.reason,
                    'metrics': verdict.metrics,
                })
                logger.warning("shard %d: family %s DROPPED as pseudo-family — %s",
                               shard_idx, urec.get('id', '?'), verdict.reason)
                continue

            if len(family_out) > 1:
                split_families += 1
            rec_added += len(family_out)
            out_records.extend(family_out)

        # One input family counts as 1 toward n_in; its surplus output+drops is `extra`.
        extra += rec_added + rec_pruned - 1

    n_out = len(out_records)
    # Worker-internal data-integrity invariant: every input family is accounted for as
    # exactly one of {kept (with its fragments), pruned}. No family silently lost.
    assert n_out + n_pruned == n_in + extra, (
        f"shard {shard_idx} family loss: n_in={n_in} extra={extra} "
        f"n_pruned={n_pruned} n_out={n_out}")

    _write_jsonl(paths['refined'], out_records)
    if pruned_records:
        _write_jsonl(paths['pruned'], pruned_records)
    with open(paths['consensus'], 'w') as fh:
        for rec in out_records:
            fh.write(f">{rec['id']}\n{rec.get('sequence', '')}\n")

    stats = {'shard': shard_idx, 'n_in': n_in, 'n_out': n_out,
             'extra': extra, 'split_families': split_families,
             'n_recalled': n_recalled, 'n_over_budget': n_over_budget,
             'n_pruned': n_pruned}
    with open(paths['stats'], 'w') as fh:
        json.dump(stats, fh)

    # .done marker LAST — a killed worker leaves no marker, so resume re-runs it.
    with open(paths['done'], 'w') as fh:
        fh.write('done\n')
    return stats


def _write_shard_fallback(shard_idx: int, shard_dir: str) -> Dict:
    """Whole-shard crash recovery: emit every input record as its original seed.

    Used only when a worker process dies before writing its outputs. Preserves the
    no-silent-loss invariant (n_out == n_in, extra == 0)."""
    paths = _shard_paths(shard_dir, shard_idx)
    records = _read_jsonl(paths['records'])
    out_records = []
    for rec in records:
        seed = dict(rec)
        seed.pop('_aln', None)
        seed.setdefault('consensus_source', 'original')
        out_records.append(seed)
    _write_jsonl(paths['refined'], out_records)
    with open(paths['consensus'], 'w') as fh:
        for rec in out_records:
            fh.write(f">{rec['id']}\n{rec.get('sequence', '')}\n")
    stats = {'shard': shard_idx, 'n_in': len(records), 'n_out': len(out_records),
             'extra': 0, 'split_families': 0, 'n_pruned': 0, 'shard_crashed': True}
    with open(paths['stats'], 'w') as fh:
        json.dump(stats, fh)
    with open(paths['done'], 'w') as fh:
        fh.write('done\n')
    return stats


# ═══════════════════════════════════════════════════════════════════════
# 5. Post-refine merge (the looser merge DEFERRED from Phase 0)
# ═══════════════════════════════════════════════════════════════════════

def merge_post_refine(records: List[Dict], config) -> Dict:
    """Single CD-HIT-EST at `post_refine_merge_identity` to collapse boundary-
    convergent duplicates newly exposed by extend-align-trim.

    Reuses phase0._run_cdhit (word size 10 for ~0.95). On CD-HIT failure the
    un-merged set is returned (no data loss; downstream Combine dedups again at 90%).
    Returns {'records': [...], 'merged': int}."""
    if len(records) < 2:
        return {'records': records, 'merged': 0}

    import shutil
    import tempfile
    tmp_dir = tempfile.mkdtemp(prefix='post_refine_merge_')
    in_path = os.path.join(tmp_dir, 'input.fa')
    out_path = os.path.join(tmp_dir, 'merged')

    # Stable representative selection: longer + higher-copy first comer wins.
    ordered = sorted(records,
                     key=lambda r: (len(r.get('sequence', '')), r.get('copies', 0)),
                     reverse=True)
    # Index by a positional key to map survivors back unambiguously (ids may now be
    # non-unique after chimera fragmentation, e.g. ..._chimfrag1).
    keyed = {}
    with open(in_path, 'w') as fh:
        for i, rec in enumerate(ordered):
            tag = f"seq{i}"
            keyed[tag] = rec
            fh.write(f">{tag}\n{rec.get('sequence', '')}\n")

    result_path = _run_cdhit(in_path, out_path,
                             config.post_refine_merge_identity, 10, config)
    if not result_path:
        logger.warning("Post-refine merge CD-HIT-EST failed; keeping un-merged set")
        shutil.rmtree(tmp_dir, ignore_errors=True)
        return {'records': records, 'merged': 0}

    surviving = set()
    with open(result_path) as fh:
        for line in fh:
            if line.startswith('>'):
                surviving.add(line[1:].strip().split()[0])
    merged_records = [keyed[t] for t in keyed if t in surviving]
    shutil.rmtree(tmp_dir, ignore_errors=True)

    merged = len(records) - len(merged_records)
    logger.info("Post-refine merge: %d -> %d sequences (@%.2f, merged %d)",
                len(records), len(merged_records),
                config.post_refine_merge_identity, merged)
    return {'records': merged_records, 'merged': merged}


# ═══════════════════════════════════════════════════════════════════════
# Phase 1 orchestrator entry
# ═══════════════════════════════════════════════════════════════════════

def run_phase1(phase0_output: Dict, config) -> Dict:
    """Shard plan → BED route → parallel per-shard refine → concat → post-merge.

    Returns {'records': [...], 'stats': {...}}. The stats dict reconciles the
    data-integrity equation from independently tracked terms:
        output == input + fragments_added − merged
    and logs every term."""
    logger.info("=" * 60)
    logger.info("Phase 1: BED-seeded extend-align-trim (sharded orchestrator)")
    logger.info("=" * 60)

    records = phase0_output['records']
    index_path = phase0_output.get('instance_index')
    input_count = len(records)

    shard_dir = os.path.join(config.temp_dir, 'phase1_shards')
    os.makedirs(shard_dir, exist_ok=True)

    # 1. Shard plan.
    shards = make_shards(records, config.shard_size)
    n_shards = len(shards)
    record_to_shard: Dict[str, int] = {}
    for i, shard in enumerate(shards):
        paths = _shard_paths(shard_dir, i)
        # Only (re)write the input records file if the shard is not already done —
        # preserves a completed shard's frozen inputs across a resume.
        if not os.path.exists(paths['done']):
            _write_jsonl(paths['records'], shard)
        for rec in shard:
            record_to_shard[rec['id']] = i
    logger.info("Sharding: %d records -> %d shards (size ~%d)",
                input_count, n_shards, config.shard_size)

    # 2. Route instances (only for shards not yet done; clear-then-route pending
    #    shards' instance files so a resume does not double-append).
    if index_path and os.path.exists(index_path):
        for i in range(n_shards):
            if not os.path.exists(_shard_paths(shard_dir, i)['done']):
                inst_file = _shard_paths(shard_dir, i)['instances']
                if os.path.exists(inst_file):
                    os.remove(inst_file)
        pending_for_route = {
            rid: s for rid, s in record_to_shard.items()
            if not os.path.exists(_shard_paths(shard_dir, s)['done'])}
        routed = route_instances(index_path, pending_for_route, shard_dir)
        logger.info("Routed instances to %d shards (%d total rows)",
                    len(routed), sum(routed.values()))
    else:
        logger.info("No instance index; all families route to fallback / seed")

    # 3. Parallel per-shard refine with resume.
    pending = [i for i in range(n_shards)
               if not os.path.exists(_shard_paths(shard_dir, i)['done'])]
    done_already = n_shards - len(pending)
    if done_already:
        logger.info("Resume: %d/%d shards already complete (skipped)",
                    done_already, n_shards)

    # 3a. Lazy genome BLAST DB — built ONCE in this (parent) process, before any
    #     worker spawns, when pending shards contain families that may touch it:
    #     under-instanced families (M5 fallback) OR recall-eligible families (N4
    #     Tier-2). The resulting db path is frozen into config_dict so every worker
    #     reads it; workers NEVER build the DB themselves (that would race on the
    #     shared db files). recall_eligible ⊆ under_instanced (it REQUIRES
    #     under_instanced to co-fire), so the under-instanced count is a sufficient
    #     and tight trigger for both consumers — no extra index pass needed.
    # rmblastn copy source / N7 need the genome DB for EVERY family (not just
    # under-instanced ones), so ensure it whenever those are on.
    rmblastn_phase1 = (getattr(config, 'enable_rmblastn_copies', False)
                       or getattr(config, 'enable_seed_chimera', False))
    db_consumer_on = (getattr(config, 'enable_fallback_recruitment', False)
                      or getattr(config, 'enable_selective_recall', False))
    if rmblastn_phase1 and pending:
        logger.info("rmblastn copy source / N7 enabled -> ensuring genome BLAST DB "
                    "(built ONCE in parent; workers never build it)")
        ensure_genome_blast_db(config)
    if db_consumer_on and pending:
        n_under = _count_under_instanced(
            index_path, record_to_shard, set(pending), config.min_copies_for_msa)
        if n_under > 0:
            logger.info("Recall/fallback: %d under-instanced families in pending "
                        "shards (< %d BED instances) -> ensuring genome BLAST DB "
                        "(built ONCE in parent; workers never build it)",
                        n_under, config.min_copies_for_msa)
            ensure_genome_blast_db(config)
        else:
            logger.info("Recall/fallback: no under-instanced families in pending "
                        "shards; genome BLAST DB not built (common path never "
                        "touches it)")

    # N4 global recall budget split across the worker pool. Each worker gets an even
    # share (ceil) of config.recall_global_family_budget; the per-worker counter caps
    # recalls in that worker. This bounds total recall families without a shared
    # cross-process counter (which a ProcessPoolExecutor cannot cheaply provide).
    config_dict = dict(config.__dict__)
    if getattr(config, 'enable_selective_recall', False) and pending:
        n_workers_eff = max(1, min(config.threads, len(pending)))
        budget = int(getattr(config, 'recall_global_family_budget', 0))
        config_dict['recall_worker_budget'] = (
            (budget + n_workers_eff - 1) // n_workers_eff if budget > 0 else 0)
    else:
        config_dict['recall_worker_budget'] = 0
    if pending:
        max_workers = max(1, min(config.threads, len(pending)))
        with ProcessPoolExecutor(max_workers=max_workers) as ex:
            futs = {ex.submit(_process_shard, i, shard_dir, config_dict): i
                    for i in pending}
            for fut in as_completed(futs):
                i = futs[fut]
                try:
                    st = fut.result()
                    logger.info("shard %d done: in=%d out=%d split_families=%d",
                                i, st['n_in'], st['n_out'], st['split_families'])
                except Exception as e:  # noqa: BLE001 — recover the whole shard
                    logger.error("shard %d worker crashed (%s); falling back to "
                                 "original seeds for this shard", i, e)
                    _write_shard_fallback(i, shard_dir)

    # 4. Concatenate refined shards + read per-shard stats (works for resumed too).
    result_records: List[Dict] = []
    pruned_records: List[Dict] = []   # N5: aggregated dropped-pseudo-family audit
    total_in = 0
    fragments_added = 0
    split_families = 0
    n_recalled = 0
    n_over_budget = 0
    n_pruned = 0
    for i in range(n_shards):
        paths = _shard_paths(shard_dir, i)
        result_records.extend(_read_jsonl(paths['refined']))
        pruned_records.extend(_read_jsonl(paths['pruned']))
        if os.path.exists(paths['stats']):
            with open(paths['stats']) as fh:
                st = json.load(fh)
            total_in += st.get('n_in', 0)
            fragments_added += st.get('extra', 0)
            split_families += st.get('split_families', 0)
            n_recalled += st.get('n_recalled', 0)
            n_over_budget += st.get('n_over_budget', 0)
            n_pruned += st.get('n_pruned', 0)

    pre_merge_count = len(result_records)
    # Cross-check: sharded inputs reconcile with the orchestrator input count and the
    # independently tracked fragment additions and pseudo-family drops (catches any
    # silent family loss). Each family is exactly one of {kept-with-fragments, pruned}.
    if total_in != input_count:
        logger.error("DATA INTEGRITY: sharded input %d != phase0 input %d",
                     total_in, input_count)
    if pre_merge_count + n_pruned != total_in + fragments_added:
        logger.error("DATA INTEGRITY: pre-merge %d + pruned %d != input %d + "
                     "fragments_added %d", pre_merge_count, n_pruned, total_in,
                     fragments_added)
    logger.info("Concatenated %d shards: input=%d, fragments added=%d "
                "(subfamily splits + chimera fragments; %d families split), "
                "pseudo-families pruned=%d -> %d records before merge",
                n_shards, total_in, fragments_added, split_families, n_pruned,
                pre_merge_count)
    # N5: persist the aggregated dropped-family audit (id + reason + metrics) — NEVER
    # mixed into the consensus output, a standalone auditable record of every drop.
    if pruned_records:
        prune_audit_path = os.path.join(shard_dir, 'pruned_pseudo_families.jsonl')
        _write_jsonl(prune_audit_path, pruned_records)
        logger.info("N5 conservative FP prune: %d pseudo-families dropped; audit at %s",
                    n_pruned, prune_audit_path)

    # 5. Post-refine merge.
    merge_out = merge_post_refine(result_records, config)
    result_records = merge_out['records']
    merged = merge_out['merged']
    output_count = len(result_records)

    # Final reconciliation (independent terms): output == input + fragments_added
    # − pruned − merged. `pruned` (N5 pseudo-family drops) is tracked from the per-shard
    # stats, an independent source from `fragments_added`, so the equation stays a real
    # audit (a silent loss anywhere still breaks it).
    expected = input_count + fragments_added - n_pruned - merged
    if output_count != expected:
        logger.error("DATA INTEGRITY: output %d != input %d + fragments_added %d "
                     "- pruned %d - merged %d (= %d)", output_count, input_count,
                     fragments_added, n_pruned, merged, expected)
    logger.info("Phase 1 complete: input=%d + fragments_added=%d - pruned=%d "
                "- merged=%d = output=%d", input_count, fragments_added, n_pruned,
                merged, output_count)
    if getattr(config, 'enable_selective_recall', False):
        logger.info("N4 selective recall: %d families recalled (Tier-2 blastn), "
                    "%d recall-eligible but over budget (kept BED copies, flagged "
                    "completeness_verified=False)", n_recalled, n_over_budget)

    stats = {
        'input_count': input_count,
        'pre_merge_count': pre_merge_count,
        'fragments_added': fragments_added,
        'split_families': split_families,
        'merged': merged,
        'n_pruned': n_pruned,
        'output_count': output_count,
        'n_shards': n_shards,
        'n_recalled': n_recalled,
        'n_over_budget': n_over_budget,
        'conservation_ok': output_count == expected and total_in == input_count,
    }
    return {'records': result_records, 'stats': stats, 'pruned': pruned_records}
