"""
N6 — representativeness / quality metric harness (REFINE_STRATEGY_DESIGN_v2.md §6).

REPORTING ONLY. Nothing here is on the critical path: the harness consumes an already-
emitted consensus library plus the run's ground-truth artifacts (the mdl-repeat BED, the
emitted records with their subfamily lineage, the N5 pruned-family audit) and computes
observable invariants on the REAL run. It invents no target numbers — every figure is
measured from the run, consistent with the project's no-fabricated-numbers discipline.

Metrics (§6):

  1. Member recovery rate (PRIMARY).  blastn the emitted consensus library back against
     the genome; per family compute the fraction of that family's known true members (BED
     instances, optionally ∪ recall-confirmed divergent members) that are recovered (a
     genome hit overlaps the instance at a coverage threshold).  The headline is the
     per-family recovery DISTRIBUTION under two library shapes:
         * multi-consensus   — the v2 output as emitted (subfamily splits kept)
         * single-consensus  — each family collapsed to ONE representative consensus
                               (its longest emitted record), forcing the pre-v2 shape
     The design predicts the DIVERGENT-family tail of the distribution lifts under multi-
     consensus.  Reported as a distribution; no target is asserted.

  2. Divergent-member recovery.  The same recovery, restricted to the divergent-tail
     instances (instance divergence in a configurable band), since that is the population
     subfamily splitting targets.

  3. Subfamily count vs family divergence.  Sanity curve — wider intra-family copy
     divergence should emit more subfamilies.  A tuning diagnostic, not a result.

  4. FP-side metric.  From the N5 pruned-family audit: count dropped, per-reason
     breakdown, and a sample for manual spot-check that none had recoverable support.

  5. Boundary metric.  Refined-vs-original consensus length distribution.

The recovery search reuses the project's blastn conventions (ensure_genome_blast_db,
strip_seqid_prefix, BED 0-based half-open coordinates).  blastn temp files are removed
after use.
"""

import json
import logging
import os
import subprocess
import tempfile
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)

# Reuse the project's seqid normalizer + lazy genome DB builder.
import sys
_here = os.path.dirname(__file__)
if _here not in sys.path:
    sys.path.insert(0, _here)
from phase1_extract import strip_seqid_prefix          # noqa: E402
from phase1_fallback import ensure_genome_blast_db      # noqa: E402


# ═══════════════════════════════════════════════════════════════════════
# Ground-truth: BED instances per family
# ═══════════════════════════════════════════════════════════════════════

@dataclass
class TrueMember:
    """One known true member of a family (a BED instance or a recall-confirmed copy)."""
    chrom: str
    start: int          # 0-based half-open
    end: int
    divergence: float   # 1 - score/1000 for BED; 1 - pident/100 for recall


def parse_bed_truth(bed_path: str) -> Dict[str, List[TrueMember]]:
    """Parse the mdl-repeat instance BED into {family_id: [TrueMember, ...]}.

    BED columns: chrom  start  end  R=N  score  strand.  Family id is normalized to the
    `mdl_R<N>` form used by the emitted records (phase0 ids families `mdl_R<N>` from the
    `R=N` header), so the truth keys line up with the consensus family keys.  Per-instance
    divergence is 1 - score/1000 (the mdl-repeat score convention)."""
    truth: Dict[str, List[TrueMember]] = defaultdict(list)
    if not bed_path or not os.path.exists(bed_path):
        return truth
    with open(bed_path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line or line.startswith('#'):
                continue
            f = line.split('\t')
            if len(f) < 4:
                continue
            try:
                chrom = f[0]
                start = int(f[1])
                end = int(f[2])
            except ValueError:
                continue
            fam = f[3]                       # e.g. "R=42"
            if fam.startswith('R='):
                fam = 'mdl_R' + fam[2:]
            score = 0.0
            if len(f) >= 5:
                try:
                    score = float(f[4])
                except ValueError:
                    score = 0.0
            div = max(0.0, 1.0 - score / 1000.0)
            truth[fam].append(TrueMember(chrom=chrom, start=start, end=end,
                                         divergence=div))
    return dict(truth)


def family_key(record_id: str) -> str:
    """Map an emitted consensus id back to its parent family id.

    Strips the v2 `_sf<k>` subfamily suffix and the `_chimfrag<j>` chimera suffix so a
    subfamily / chimera-fragment consensus groups under its parent `mdl_R<N>` (the BED
    truth key).  A single-cluster family keeps its id verbatim, so this is a no-op for it.
    Suffixes compose (`mdl_R8_sf2_chimfrag1`), so both are stripped."""
    fid = record_id
    for marker in ('_chimfrag', '_sf'):
        idx = fid.find(marker)
        while idx != -1:
            # strip from the FIRST occurrence of the leftmost-applied suffix
            fid = fid[:idx]
            idx = fid.find(marker)
    return fid


# ═══════════════════════════════════════════════════════════════════════
# blastn the consensus library against the genome
# ═══════════════════════════════════════════════════════════════════════

@dataclass
class Hit:
    chrom: str
    start: int          # 0-based half-open
    end: int
    pident: float


def blast_consensi_to_genome(consensi: List[Tuple[str, str]], config,
                             pident_min: float = 70.0,
                             wall_s: int = 3600) -> Dict[str, List[Hit]]:
    """blastn every consensus (id, sequence) against the whole-genome DB.

    Returns {consensus_id: [Hit, ...]}.  Builds the genome DB once via the project's
    lazy builder (ensure_genome_blast_db).  One blastn call over the whole library (the
    query is multi-FASTA), so the genome DB is touched ONLY here — the harness is off the
    common refine path entirely.  Returns {} (and logs) on any blastn failure; never
    fabricates hits.  Temp files are removed."""
    if not consensi:
        return {}
    db = ensure_genome_blast_db(config)
    if not db:
        logger.error("metrics: genome BLAST DB unavailable; cannot compute recovery")
        return {}

    tmp_dir = tempfile.mkdtemp(prefix='metrics_blast_')
    query_path = os.path.join(tmp_dir, 'consensi.fa')
    try:
        with open(query_path, 'w') as fh:
            for cid, seq in consensi:
                if seq:
                    fh.write(f">{cid}\n{seq}\n")
        cmd = [
            config.blastn_exe,
            '-query', query_path,
            '-db', db,
            '-outfmt', '6 qseqid sseqid sstart send pident length',
            '-evalue', str(getattr(config, 'blastn_evalue', 1e-5)),
            '-num_threads', str(max(1, int(getattr(config, 'threads', 1)))),
            # The processed genome is SOFT-masked (lowercase low-complexity/repeat
            # regions).  TE consensi map almost entirely onto those very regions, so
            # blastn's default subject soft-masking (`-dust yes` + soft-masking of
            # lowercase) would zero out legitimate recovery hits.  `-dust no` +
            # `-soft_masking false` make recovery measure TRUE homology against the
            # masked genome rather than an artifact of the mask.
            '-dust', 'no',
            '-soft_masking', 'false',
        ]
        try:
            result = subprocess.run(cmd, capture_output=True, text=True,
                                    timeout=wall_s)
        except Exception as e:  # noqa: BLE001 — report, never fabricate
            logger.error("metrics: blastn failed/timed-out (%s)", e)
            return {}
        if result.returncode != 0:
            logger.error("metrics: blastn exit %d: %s",
                         result.returncode, (result.stderr or '')[:300])
            return {}

        hits: Dict[str, List[Hit]] = defaultdict(list)
        for line in result.stdout.strip().split('\n'):
            if not line:
                continue
            f = line.split('\t')
            if len(f) < 5:
                continue
            try:
                qid = f[0]
                sseqid = strip_seqid_prefix(f[1])
                sstart = int(f[2])
                send = int(f[3])
                pident = float(f[4])
            except (ValueError, IndexError):
                continue
            if pident < pident_min:
                continue
            lo, hi = (sstart, send) if sstart <= send else (send, sstart)
            hits[qid].append(Hit(chrom=sseqid, start=lo - 1, end=hi, pident=pident))
        return dict(hits)
    finally:
        import shutil
        shutil.rmtree(tmp_dir, ignore_errors=True)


# ═══════════════════════════════════════════════════════════════════════
# Recovery computation
# ═══════════════════════════════════════════════════════════════════════

def _overlap_bp(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    return max(0, min(a_end, b_end) - max(a_start, b_start))


def _member_recovered(member: TrueMember, hits: List[Hit],
                     min_cov: float) -> bool:
    """A true member is recovered if any genome hit (on its chrom) overlaps it by
    >= min_cov of the member's length.  Coverage is of the MEMBER (the thing we want to
    recover), so a long consensus hit fully covering a short instance counts."""
    mlen = max(1, member.end - member.start)
    for h in hits:
        if h.chrom != member.chrom:
            continue
        ov = _overlap_bp(member.start, member.end, h.start, h.end)
        if ov / mlen >= min_cov:
            return True
    return False


@dataclass
class FamilyRecovery:
    family: str
    n_members: int
    n_recovered_multi: int
    n_recovered_single: int
    n_subfamilies: int
    mean_member_divergence: float
    # divergent-tail subset
    n_divergent_members: int = 0
    n_divergent_recovered_multi: int = 0
    n_divergent_recovered_single: int = 0

    @property
    def recovery_multi(self) -> float:
        return self.n_recovered_multi / self.n_members if self.n_members else 0.0

    @property
    def recovery_single(self) -> float:
        return self.n_recovered_single / self.n_members if self.n_members else 0.0


def compute_member_recovery(records: List[Dict], truth: Dict[str, List[TrueMember]],
                           config, min_cov: float = 0.5, pident_min: float = 70.0,
                           divergent_band: Tuple[float, float] = (0.12, 0.50)
                           ) -> Dict:
    """Primary metric (§6): per-family member recovery, MULTI vs SINGLE consensus.

    `records` are the emitted Phase 1 records (each a dict with 'id', 'sequence', and,
    for split families, an '_sf'-suffixed id + a 'subfamily' lineage dict).  `truth` is
    parse_bed_truth output.

    multi  = blast ALL emitted consensi; a family's recovered members = those overlapped
             by ANY of its (sub)consensi's hits (the v2 library as shipped).
    single = for each family, keep ONLY its longest emitted consensus (forcing the
             pre-v2 one-per-family shape); recompute recovery from that single record's
             hits.  Same genome, same thresholds — an apples-to-apples ablation.

    The divergent band restricts a second recovery figure to the instance-divergence tail
    that single-consensus is expected to miss.

    Returns a dict with the per-family table and aggregate distributions.  All numbers are
    measured from the blastn run; none invented."""
    # Group emitted records by parent family.
    by_family: Dict[str, List[Dict]] = defaultdict(list)
    for rec in records:
        by_family[family_key(rec['id'])].append(rec)

    # One blastn over the whole library (multi shape) — all (sub)consensi.
    all_consensi = [(rec['id'], rec.get('sequence', '')) for rec in records]
    multi_hits = blast_consensi_to_genome(all_consensi, config,
                                          pident_min=pident_min)
    if not multi_hits and all_consensi:
        # DB or blastn failed — surface, do not fabricate.
        return {'error': 'blastn_failed_or_no_hits', 'families': [], 'n_families': 0}

    # Single-consensus shape: longest record per family (the forced representative).
    single_reps: Dict[str, Dict] = {}
    for fam, recs in by_family.items():
        rep = max(recs, key=lambda r: len(r.get('sequence', '')))
        single_reps[fam] = rep

    dband_lo, dband_hi = divergent_band
    fam_rows: List[FamilyRecovery] = []

    for fam, members in truth.items():
        if fam not in by_family:
            # Family present in BED truth but absent from the emitted library (e.g. it
            # failed Phase 0 hard filters or was N5-pruned). Recorded separately by the
            # caller via FP-side metric; not a recovery row.
            continue
        recs = by_family[fam]
        n_sf = len(recs)

        # multi: union of hits across this family's (sub)consensi.
        fam_multi_hits: List[Hit] = []
        for rec in recs:
            fam_multi_hits.extend(multi_hits.get(rec['id'], []))
        # single: hits of the longest representative only.
        rep_id = single_reps[fam]['id']
        fam_single_hits = list(multi_hits.get(rep_id, []))

        n_rec_multi = n_rec_single = 0
        n_div = n_div_multi = n_div_single = 0
        for m in members:
            rec_m = _member_recovered(m, fam_multi_hits, min_cov)
            rec_s = _member_recovered(m, fam_single_hits, min_cov)
            n_rec_multi += int(rec_m)
            n_rec_single += int(rec_s)
            if dband_lo <= m.divergence <= dband_hi:
                n_div += 1
                n_div_multi += int(rec_m)
                n_div_single += int(rec_s)

        mean_div = (sum(m.divergence for m in members) / len(members)
                    if members else 0.0)
        fam_rows.append(FamilyRecovery(
            family=fam, n_members=len(members),
            n_recovered_multi=n_rec_multi, n_recovered_single=n_rec_single,
            n_subfamilies=n_sf, mean_member_divergence=round(mean_div, 4),
            n_divergent_members=n_div,
            n_divergent_recovered_multi=n_div_multi,
            n_divergent_recovered_single=n_div_single))

    # ── Aggregate distributions (measured) ────────────────────────────────
    def _dist(vals: List[float]) -> Dict:
        if not vals:
            return {'n': 0}
        s = sorted(vals)
        n = len(s)
        def q(p):
            return s[min(n - 1, int(p * n))]
        return {
            'n': n,
            'mean': round(sum(s) / n, 4),
            'median': round(q(0.5), 4),
            'p10': round(q(0.10), 4),
            'p25': round(q(0.25), 4),
            'p75': round(q(0.75), 4),
            'p90': round(q(0.90), 4),
            'min': round(s[0], 4),
            'max': round(s[-1], 4),
        }

    multi_vals = [r.recovery_multi for r in fam_rows]
    single_vals = [r.recovery_single for r in fam_rows]

    # Divergent-family tail: families that actually split (n_subfamilies > 1) are where
    # the design predicts the lift; report their recovery distribution separately.
    split_rows = [r for r in fam_rows if r.n_subfamilies > 1]
    split_multi = [r.recovery_multi for r in split_rows]
    split_single = [r.recovery_single for r in split_rows]

    # Divergent-member recovery (§6 bullet 2): pooled over the divergent band.
    tot_div = sum(r.n_divergent_members for r in fam_rows)
    tot_div_multi = sum(r.n_divergent_recovered_multi for r in fam_rows)
    tot_div_single = sum(r.n_divergent_recovered_single for r in fam_rows)

    return {
        'n_families': len(fam_rows),
        'n_split_families': len(split_rows),
        'min_coverage': min_cov,
        'pident_min': pident_min,
        'divergent_band': list(divergent_band),
        'recovery_multi_dist': _dist(multi_vals),
        'recovery_single_dist': _dist(single_vals),
        'split_family_recovery_multi_dist': _dist(split_multi),
        'split_family_recovery_single_dist': _dist(split_single),
        'divergent_member_recovery': {
            'n_divergent_members': tot_div,
            'recovered_multi': tot_div_multi,
            'recovered_single': tot_div_single,
            'rate_multi': round(tot_div_multi / tot_div, 4) if tot_div else None,
            'rate_single': round(tot_div_single / tot_div, 4) if tot_div else None,
        },
        'families': [
            {
                'family': r.family,
                'n_members': r.n_members,
                'n_subfamilies': r.n_subfamilies,
                'mean_member_divergence': r.mean_member_divergence,
                'recovery_multi': round(r.recovery_multi, 4),
                'recovery_single': round(r.recovery_single, 4),
                'n_divergent_members': r.n_divergent_members,
                'divergent_recovery_multi': (
                    round(r.n_divergent_recovered_multi / r.n_divergent_members, 4)
                    if r.n_divergent_members else None),
                'divergent_recovery_single': (
                    round(r.n_divergent_recovered_single / r.n_divergent_members, 4)
                    if r.n_divergent_members else None),
            }
            for r in sorted(fam_rows, key=lambda x: x.recovery_multi)
        ],
    }


# ═══════════════════════════════════════════════════════════════════════
# Subfamily-count vs divergence (sanity curve)
# ═══════════════════════════════════════════════════════════════════════

def subfamily_count_vs_divergence(records: List[Dict]) -> List[Dict]:
    """Sanity diagnostic (§6 bullet 3): per family, (mean intra-family divergence proxy,
    n_subfamilies).  The intra-family divergence proxy is taken from the subfamily
    lineage's `intra_identity` when present (1 - mean intra-identity), else 0 for a
    single-cluster family.  A wider-divergence family should emit more subfamilies; this
    table lets the caller eyeball the monotonicity (a flat curve = cut too coarse, a
    runaway curve = too fine)."""
    by_family: Dict[str, List[Dict]] = defaultdict(list)
    for rec in records:
        by_family[family_key(rec['id'])].append(rec)

    rows = []
    for fam, recs in by_family.items():
        n_sf = len(recs)
        intra_ids = [r['subfamily']['intra_identity'] for r in recs
                     if isinstance(r.get('subfamily'), dict)
                     and 'intra_identity' in r['subfamily']]
        if intra_ids:
            div_proxy = round(1.0 - sum(intra_ids) / len(intra_ids), 4)
        else:
            div_proxy = 0.0
        rows.append({'family': fam, 'n_subfamilies': n_sf,
                     'intra_divergence_proxy': div_proxy})
    return sorted(rows, key=lambda r: r['intra_divergence_proxy'])


# ═══════════════════════════════════════════════════════════════════════
# FP-side metric (from the N5 pruned-family audit)
# ═══════════════════════════════════════════════════════════════════════

def fp_side_metric(pruned_records: List[Dict], n_input_families: int,
                  sample_n: int = 10) -> Dict:
    """FP-side metric (§6 bullet 4): drop rate + per-reason breakdown + a spot-check
    sample.  `pruned_records` is the aggregated N5 audit (each: id, reason, metrics).
    Reports the fraction of families dropped and a per-reason count; the sample is for a
    human/audit pass to confirm none had recoverable support (the conservative-pruning
    verification — the harness surfaces it, the audit confirms it)."""
    n_dropped = len(pruned_records)
    by_reason: Dict[str, int] = defaultdict(int)
    for pr in pruned_records:
        reason = str(pr.get('reason', 'unspecified')).split(':')[0].strip()
        by_reason[reason] += 1
    return {
        'n_input_families': n_input_families,
        'n_dropped': n_dropped,
        'drop_rate': round(n_dropped / n_input_families, 6) if n_input_families else 0.0,
        'by_reason': dict(by_reason),
        'spot_check_sample': pruned_records[:sample_n],
    }


# ═══════════════════════════════════════════════════════════════════════
# Boundary metric (refined vs original length distribution)
# ═══════════════════════════════════════════════════════════════════════

def boundary_metric(records: List[Dict],
                   original_lengths: Dict[str, int]) -> Dict:
    """Boundary metric (§6 bullet 5): refined-vs-original consensus length distribution.

    `original_lengths` maps family id -> original mdl-repeat seed length (from Phase 0
    input).  For each emitted record we compare its refined length to the parent family's
    original length.  consensus_source=='original' records (never-regress fallbacks) are
    counted separately — they are by construction unchanged in length, so the refined-vs-
    original delta is meaningful only on the truly-refined subset."""
    refined_lens = []
    original_lens = []
    deltas = []
    n_refined = n_original_source = 0
    for rec in records:
        rlen = len(rec.get('sequence', ''))
        fam = family_key(rec['id'])
        olen = original_lengths.get(fam)
        src = rec.get('consensus_source', 'original')
        if src == 'original':
            n_original_source += 1
        else:
            n_refined += 1
        if olen:
            refined_lens.append(rlen)
            original_lens.append(olen)
            deltas.append(rlen - olen)

    def _summary(vals):
        if not vals:
            return {'n': 0}
        s = sorted(vals)
        n = len(s)
        return {'n': n, 'mean': round(sum(s) / n, 1),
                'median': s[n // 2], 'min': s[0], 'max': s[-1]}

    return {
        'n_records': len(records),
        'n_refined_source': n_refined,
        'n_original_source': n_original_source,
        'refined_length': _summary(refined_lens),
        'original_length': _summary(original_lens),
        'length_delta': _summary(deltas),
    }


# ═══════════════════════════════════════════════════════════════════════
# Top-level harness
# ═══════════════════════════════════════════════════════════════════════

def run_metrics(records: List[Dict], config,
               bed_path: Optional[str] = None,
               pruned_records: Optional[List[Dict]] = None,
               original_lengths: Optional[Dict[str, int]] = None,
               n_input_families: Optional[int] = None,
               min_cov: float = 0.5, pident_min: float = 70.0) -> Dict:
    """Run the full §6 metric harness on an emitted library and write a JSON report.

    All inputs come from the completed run (records = Phase 1 output, bed_path = the
    mdl-repeat instance BED used as ground truth, pruned_records = the N5 audit,
    original_lengths = Phase 0 seed lengths).  Returns the metric dict; the caller writes
    it (or pass through to write_metrics_report).  No critical-path side effects."""
    bed_path = bed_path or getattr(config, 'bed_file', '') or ''
    truth = parse_bed_truth(bed_path)

    report: Dict = {
        'n_emitted_records': len(records),
        'n_truth_families': len(truth),
        'bed_path': bed_path,
    }
    report['member_recovery'] = compute_member_recovery(
        records, truth, config, min_cov=min_cov, pident_min=pident_min)
    report['subfamily_count_vs_divergence'] = subfamily_count_vs_divergence(records)
    if pruned_records is not None:
        n_in = (n_input_families if n_input_families is not None
                else len(truth) + len(pruned_records))
        report['fp_side'] = fp_side_metric(pruned_records, n_in)
    if original_lengths is not None:
        report['boundary'] = boundary_metric(records, original_lengths)
    return report


def write_metrics_report(report: Dict, out_path: str) -> str:
    """Write the metric report as JSON.  Returns the path."""
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    with open(out_path, 'w') as fh:
        json.dump(report, fh, indent=2, default=str)
    logger.info("Metrics report written: %s", out_path)
    return out_path
