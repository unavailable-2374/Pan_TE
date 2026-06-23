"""
Phase 1 — N5: conservative false-positive pruning (the pseudo-family discriminant).

This module operationalizes goal-2 of the v2 redesign (REFINE_STRATEGY_DESIGN_v2.md
§4) as a SINGLE, sensitivity-weighted drop decision wired into the per-family body
AFTER recall (N4) and clustering (N2/N3) have produced the family's assembled copy
set and its qualifying-cluster verdict.  The governing principle, made concrete in
§4.1/§4.2:

    A family is PSEUDO  iff  recall cannot assemble a coherent, connected copy set —
    i.e. its purported copies form NO qualifying cluster AND are mutually incoherent
    (pairwise un-homologous), so no cut produces a real subfamily.

The decision rests on POSITIVE evidence of incoherence, NEVER on low copy count.
This distinction is the whole point on chr4, where 84% of *real* families are
low-copy (median 3 copies): low copy count is the NORM, not a noise signal
(§4.1 pseudo-family row, §4.2).  So the discriminant is a conjunction that a real
low-copy family CANNOT satisfy:

    is_pseudo  ⇔  enable_conservative_fp_prune                       (master switch)
              AND  fallback_single_cluster                            (§4.2: ZERO
                       qualifying clusters — no cut produced a real subfamily)
              AND  n_copies >= pseudo_family_min_copies_for_call      (enough copies
                       to TRUST an incoherence call; a 1-2-copy family is never judged
                       incoherent — too little evidence, kept by default)
              AND  median pairwise identity < pseudo_family_max_intra_identity
                       (the copies are mutually un-related — a grab-bag, NOT a tight
                        low-copy family)
              AND  homologous-pair fraction < pseudo_family_min_homologous_pair_frac
                       (most copy PAIRS share no significant local alignment at all —
                        the all-vs-all BLASTN finds almost no connecting HSP)
              AND  n_seed_coherent < seed_coherence_min                (§4.3 seed-coherence
                       protection: the family's OWN consensus/seed recruits FEWER than
                       seed_coherence_min of its copies at high identity+coverage.  If the
                       seed strongly matches even ONE copy the family has real signal and
                       is KEPT despite copy×copy incoherence — the star-topology divergent
                       family and the "real consensus + one anomalous instance" case, e.g.
                       chr4 R=286, are exactly this shape and must NOT be误删.  Computed
                       only when conjuncts 1-5 already hold, so it adds no cost to KEEP.)

Why a real low-copy family is SAFE (the core invariant, tested in tests/test_n5_fp.py):
  * A coherent 2-3-copy family has high mutual BLASTN local identity → its
    `fallback_single_cluster` cluster has high `mean_identity`, its median pairwise
    identity is high, and its homologous-pair fraction is ~1.0.  Two of the three
    incoherence conditions (low median identity, low homologous-pair fraction) FAIL,
    so it is KEPT.  Low copy count alone never trips the drop — by construction the
    drop needs POSITIVE incoherence on TOP of zero-qualifying-clusters.
  * A 1-copy or sub-`min_copies_for_call` family is KEPT unconditionally: too little
    evidence to declare incoherence (§4.3 "when in doubt, keep").

Why this is the §4.1 "key" discriminant and not the others:
  * Low-complexity / SSR is already removed in Phase 0 hard_filter (DUST + Shannon
    entropy) and Phase 2 QC; N5 does NOT re-implement it (§4.1 row 1).  `audit` here
    only ever sees families that survived those gates, so a residual low-complexity
    family is already rare and is caught by the same incoherence logic (random
    low-complexity copies do not mutually align cleanly) without a second filter.
  * Over-fragmented short pieces are re-joined by Phase 0 BED co-occurrence assembly;
    a residual short fragment that STILL recruits a coherent copy set is KEPT (it may
    be a real solo-LTR / MITE).  N5 never drops on length — only on incoherence
    (§4.1 row 2, §4.3).  A short fragment with no coherent copy set fails the SAME
    incoherence conjunction and is handled by it; no length rule is added.
  * Chimera (A+B) is SPLIT by M3, not dropped (§4.1 row 4); N5 leaves M3 untouched.

Conservative boundary (§4.3), enforced structurally:
  * An over-budget family (recall could not be run, `completeness_verified=False`) is
    NEVER dropped here — absence of verification is not evidence of absence.  Such a
    family is kept and its flag already records the honest "not verified" state.
    N5 checks this explicitly and abstains.
  * The drop needs POSITIVE incoherence evidence (all conjuncts), never a single weak
    signal, never low copy count.
  * Refine deliberately UNDER-prunes: downstream Combine (CD-HIT-EST 90% across all
    three sources) and TE-looker's own gates provide further FP control, so a borderline
    family is kept and handed off in the SENSITIVE set (§4.3 last bullet).

Counting / audit (§4 deliverable 4, §8.2):
  * A dropped family is removed from the normal output and recorded in a SEPARATE
    audit list (id + reason + metrics).  The caller (`_process_shard`) counts drops in
    an explicit `n_pruned` term so the data-integrity reconciliation stays a real
    cross-check: `n_out == n_in + extra − pruned`.  A silent loss still trips it.

Safe-disable (degrade to N4):  `enable_conservative_fp_prune=False` makes
`pseudo_family_verdict` ALWAYS return not-pseudo, so the per-family body emits exactly
what N4 produced.  The N1-N4 suites stay green unchanged.

See REFINE_STRATEGY_DESIGN_v2.md §4 (all), §4.1 pseudo-family row, §4.2, §4.3, §8.2,
§10 N5.
"""

import logging
import os
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field
from typing import Dict, List, Optional

import numpy as np

from phase1_cluster import pairwise_distance

logger = logging.getLogger(__name__)


# ═══════════════════════════════════════════════════════════════════════
# Verdict structure
# ═══════════════════════════════════════════════════════════════════════

@dataclass
class PseudoFamilyVerdict:
    """Per-family pseudo-family (drop) verdict (N5).

    Attributes
    ----------
    is_pseudo : bool
        True iff EVERY incoherence conjunct holds (see module docstring) — only then
        is the family dropped.  Default False: the discriminant defaults to KEEP and
        requires positive evidence to flip (§4.3).
    reason : str
        Human-readable one-line justification (audit only; never a FASTA header).
    metrics : Dict
        The numbers the call rested on — n_copies, median pairwise identity,
        homologous-pair fraction, fallback_single_cluster, the thresholds used — so a
        dropped-family audit can be reviewed against the raw evidence.
    """
    is_pseudo: bool = False
    reason: str = ""
    metrics: Dict = field(default_factory=dict)


# ═══════════════════════════════════════════════════════════════════════
# Incoherence metrics over the copy×copy distance matrix
# ═══════════════════════════════════════════════════════════════════════

def _pairwise_identity_stats(copies: List, config) -> Dict:
    """Median pairwise identity + homologous-pair fraction over the copy set.

    Reuses `phase1_cluster.pairwise_distance` (all-vs-all BLASTN local identity, the
    SAME memoized matrix the clustering pass just computed for this family — so this
    adds NO extra BLAST).  Two derived statistics:

      * median_identity — median over the (n choose 2) off-diagonal entries of
        (1 - distance).  Median (not mean) so one accidental high-identity pair in an
        otherwise unrelated grab-bag cannot pull the score up; a real family has the
        whole distribution high, a grab-bag the whole distribution near 0.
      * homologous_pair_frac — fraction of off-diagonal pairs with ANY gate-passing HSP
        (distance < 1.0; a pair with no HSP is left at distance exactly 1.0 by
        `_blastn_identity_distance`).  This is the "are the copies even CONNECTED?"
        signal §4.1/§4.2 turns on: in a real family almost every pair shares homology;
        in a grab-bag almost no pair does.

    n_copies < 2 → identity 1.0 / frac 1.0 (a single copy is trivially self-coherent;
    no pairwise evidence exists to call it incoherent).  Returns a dict of native
    floats (no numpy scalars leak to JSON).
    """
    n = len(copies)
    if n < 2:
        return {'n_copies': n, 'median_identity': 1.0,
                'homologous_pair_frac': 1.0, 'n_pairs': 0}

    dist = pairwise_distance(copies, config)
    # Off-diagonal upper-triangle entries only (unordered pairs, no double count).
    iu = np.triu_indices(n, k=1)
    pair_dists = dist[iu]
    if pair_dists.size == 0:
        return {'n_copies': n, 'median_identity': 1.0,
                'homologous_pair_frac': 1.0, 'n_pairs': 0}

    identities = 1.0 - pair_dists
    median_identity = float(np.median(identities))
    # A pair with NO gate-passing HSP sits at distance == 1.0 (identity 0.0); count the
    # complement as "homologous".  Strict < 1.0 so the no-HSP sentinel is excluded.
    homologous = int(np.count_nonzero(pair_dists < 1.0))
    n_pairs = int(pair_dists.size)
    return {
        'n_copies': n,
        'median_identity': median_identity,
        'homologous_pair_frac': homologous / n_pairs if n_pairs else 1.0,
        'n_pairs': n_pairs,
    }


# ═══════════════════════════════════════════════════════════════════════
# Seed-coherence: does the family's own consensus recruit any of its copies?
# ═══════════════════════════════════════════════════════════════════════

def _seed_coherent_copies(rec: Dict, copies: List, config) -> Optional[int]:
    """Count how many copies the family's OWN seed/consensus strongly matches.

    The copy×copy incoherence metric answers "are the copies coherent WITH EACH OTHER?".
    That is the wrong question for a star-topology divergent family (every copy matches
    the seed but pairwise-diverged copies do not match each other) or a real family
    carrying one anomalous instance (chr4 R=286: 2579 bp consensus; copies
    2578/2578/2578/41659/7307 bp, mutually un-homologous, yet ONE copy is 99.8% identical
    to the consensus).  This helper asks the RIGHT question for the keep decision: does
    the family's consensus itself recruit ≥1 of its purported copies?  If so, the family
    has positive real signal and must NOT be dropped (§4.3 存疑即留).

    Implementation: blastn the seed (rec['sequence'], used as the query) against a tiny
    per-family DB of the copy sequences.  A copy counts as "seed-coherent" when SOME HSP
    clears both `seed_coherence_min_pident` (%identity) AND `seed_coherence_min_coverage`
    of the SEED length (aln_len / seed_len) — coverage measured against the seed, NOT the
    copy, so a copy that is much LONGER than the consensus (e.g. R=286's 7307 bp copy, or
    the 41659 bp 16×-consensus instance) still counts as long as the consensus span is
    well aligned.  O(copies) HSPs on a copy set of a handful of sequences — never the
    genome.

    Returns the integer count of seed-coherent copies, or None on any tool failure
    (missing binary / non-zero exit) so the caller degrades to KEEP (a failure to verify
    the seed is NOT evidence the family is spurious — §4.3, mirrors pairwise_distance's
    never-fabricate contract).  An empty/absent seed or empty copy set returns None.
    """
    seed = (rec.get('sequence') or '').strip().replace('\n', '')
    if not seed or not copies:
        return None
    seed_len = len(seed)
    if seed_len <= 0:
        return None

    makeblastdb = getattr(config, 'makeblastdb_exe', 'makeblastdb') or 'makeblastdb'
    blastn = getattr(config, 'blastn_exe', 'blastn') or 'blastn'
    evalue = str(getattr(config, 'blastn_evalue', 1e-5))
    min_pident = getattr(config, 'seed_coherence_min_pident', 80.0)
    min_cov = getattr(config, 'seed_coherence_min_coverage', 0.50)

    tmp_dir = tempfile.mkdtemp(prefix='n5_seedcoh_')
    qfa = os.path.join(tmp_dir, 'seed.fa')
    dfa = os.path.join(tmp_dir, 'copies.fa')
    db = os.path.join(tmp_dir, 'db')
    try:
        with open(qfa, 'w') as fh:
            fh.write(f">seed\n{seed}\n")
        with open(dfa, 'w') as fh:
            for i, c in enumerate(copies):
                fh.write(f">c{i}\n{c.sequence}\n")

        mk = subprocess.run([makeblastdb, '-in', dfa, '-dbtype', 'nucl', '-out', db],
                            capture_output=True, text=True, timeout=600)
        if mk.returncode != 0:
            logger.warning("seed_coherence(blastn): makeblastdb exit %d: %s",
                           mk.returncode, (mk.stderr or '')[:200])
            return None

        # word_size 7 + dust off: same sensitive regime as the copy×copy pass, so a
        # divergent-but-real seed→copy match in the ~70-80% band is not missed.
        cmd = [blastn, '-query', qfa, '-db', db,
               '-outfmt', '6 sseqid pident length',
               '-evalue', evalue, '-word_size', '7', '-dust', 'no']
        br = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        if br.returncode != 0:
            logger.warning("seed_coherence(blastn): blastn exit %d: %s",
                           br.returncode, (br.stderr or '')[:200])
            return None

        # A copy is seed-coherent if it has SOME HSP clearing both gates. Best HSP per
        # copy by (pident, coverage); accumulate aligned length is unsafe (overlapping
        # HSPs), so use the single best HSP's coverage of the seed — conservative.
        coherent = set()
        for line in br.stdout.splitlines():
            parts = line.split('\t')
            if len(parts) < 3:
                continue
            sid, pid, ln = parts[0], parts[1], parts[2]
            try:
                pident = float(pid)
                aln = int(ln)
            except ValueError:
                continue
            cov = aln / seed_len if seed_len else 0.0
            if pident >= min_pident and cov >= min_cov:
                coherent.add(sid)
        return len(coherent)
    except FileNotFoundError as e:
        logger.warning("seed_coherence(blastn): binary not found (%s)", e)
        return None
    except Exception as e:  # noqa: BLE001 — report + degrade to KEEP, never fabricate
        logger.warning("seed_coherence(blastn): failed (%s)", e)
        return None
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


# ═══════════════════════════════════════════════════════════════════════
# The pseudo-family discriminant (the N5 drop decision)
# ═══════════════════════════════════════════════════════════════════════

def pseudo_family_verdict(rec: Dict, copies: List, audit: Dict,
                          config) -> PseudoFamilyVerdict:
    """Decide whether a family is a pseudo-family (to be dropped) — conservatively.

    Parameters
    ----------
    rec : Dict
        The family record (carries `id`, and — if N4 ran — `completeness_verified`).
    copies : List[phase1_extract.Copy]
        The assembled copy set handed to clustering (BED ± recalled), the SAME list N2
        clustered.  The incoherence metric is measured over THIS set.
    audit : Dict
        The audit dict returned by `phase1_cluster.qualifying_clusters` for this family.
        Supplies `fallback_single_cluster` (True ⇔ ZERO qualifying clusters, §4.2) and
        `n_qualifying`.
    config : RefinerMdlConfig

    Returns
    -------
    PseudoFamilyVerdict
        `is_pseudo=True` ONLY when every incoherence conjunct holds.  Defaults to KEEP.

    The conjunction (module docstring) — each conjunct independently necessary:
      1. master switch on
      2. fallback_single_cluster (no qualifying cluster formed)
      3. n_copies >= pseudo_family_min_copies_for_call (enough evidence to judge)
      4. median pairwise identity < pseudo_family_max_intra_identity (mutually unrelated)
      5. homologous_pair_frac < pseudo_family_min_homologous_pair_frac (un-connected)
      6. n_seed_coherent < seed_coherence_min (the family's own seed recruits no coherent
         copy — §4.3 seed-coherence protection; computed only after 1-5 hold)
    PLUS the explicit §4.3 abstention: a family flagged completeness_verified=False
    (recall could not run — over budget / failed) is NEVER dropped (absence of
    verification ≠ absence of support).
    """
    rid = rec.get('id', '?')

    # ── Master switch: degrade to N4 (never drop) ────────────────────────────
    if not getattr(config, 'enable_conservative_fp_prune', True):
        return PseudoFamilyVerdict(
            is_pseudo=False, reason="conservative FP prune disabled",
            metrics={'enabled': False})

    fallback = bool(audit.get('fallback_single_cluster', False))
    n_qualifying = int(audit.get('n_qualifying', 0))

    # ── §4.3 abstention: recall could not be run → keep + already-flagged ────
    # completeness_verified is set by N4 ONLY on the over-budget / failed path
    # (==False); a successfully recalled family has it True; a family that never
    # needed recall has the key absent.  We abstain ONLY when it is explicitly False.
    if rec.get('completeness_verified', None) is False:
        return PseudoFamilyVerdict(
            is_pseudo=False,
            reason="recall not verified (over budget / failed) — kept, not dropped (§4.3)",
            metrics={'completeness_verified': False,
                     'fallback_single_cluster': fallback})

    # ── Conjunct 2: a qualifying cluster formed ⇒ coherent ⇒ never pseudo ────
    if not fallback:
        return PseudoFamilyVerdict(
            is_pseudo=False,
            reason=f"{n_qualifying} qualifying cluster(s) formed — coherent, kept",
            metrics={'fallback_single_cluster': False,
                     'n_qualifying': n_qualifying})

    # ── Measure incoherence over the copy set (reuses the memoized BLAST) ────
    stats = _pairwise_identity_stats(copies, config)
    n_copies = stats['n_copies']
    median_identity = stats['median_identity']
    homologous_pair_frac = stats['homologous_pair_frac']

    min_copies_for_call = getattr(config, 'pseudo_family_min_copies_for_call', 5)
    max_intra_identity = getattr(config, 'pseudo_family_max_intra_identity', 0.45)
    min_homologous_frac = getattr(config, 'pseudo_family_min_homologous_pair_frac', 0.25)

    metrics = {
        'fallback_single_cluster': True,
        'n_qualifying': n_qualifying,
        'n_copies': n_copies,
        'median_identity': round(median_identity, 4),
        'homologous_pair_frac': round(homologous_pair_frac, 4),
        'n_pairs': stats['n_pairs'],
        'thresholds': {
            'min_copies_for_call': min_copies_for_call,
            'max_intra_identity': max_intra_identity,
            'min_homologous_pair_frac': min_homologous_frac,
        },
    }

    # ── Conjunct 3: too few copies to TRUST an incoherence call → KEEP ───────
    # This is the explicit "never judge a low-copy family" guard.  A real low-copy
    # family (the chr4 84% norm) lands here and is kept regardless of its (noisy)
    # pairwise identity — we simply refuse to call incoherence on thin evidence.
    if n_copies < min_copies_for_call:
        return PseudoFamilyVerdict(
            is_pseudo=False,
            reason=(f"only {n_copies} copies (< {min_copies_for_call} needed to judge "
                    f"incoherence) — low copy count is NOT a noise signal, kept (§4.1/§4.3)"),
            metrics=metrics)

    # ── Conjuncts 4 & 5: POSITIVE incoherence evidence required ──────────────
    incoherent_identity = median_identity < max_intra_identity
    unconnected = homologous_pair_frac < min_homologous_frac

    if incoherent_identity and unconnected:
        # ── Conjunct 6 (seed-coherence protection, §4.3): a final NECESSARY guard ──
        # The copies are mutually incoherent — but that does NOT prove the family is
        # spurious.  A star-topology divergent family (copies match the SEED, not each
        # other) or a real family with one anomalous instance (chr4 R=286) is exactly
        # this shape.  Before dropping, ask whether the family's OWN consensus recruits
        # any of its copies; if it recruits >= seed_coherence_min, the family has real
        # signal and is KEPT.  Computed ONLY here (short-circuit), so it costs the
        # KEEP-majority nothing.  A tool failure returns None → treated as "cannot rule
        # out signal" → KEEP (never fabricate a drop).
        if getattr(config, 'enable_seed_coherence_protection', True):
            seed_coherence_min = getattr(config, 'seed_coherence_min', 1)
            n_seed_coherent = _seed_coherent_copies(rec, copies, config)
            metrics['seed_coherence_min'] = seed_coherence_min
            metrics['n_seed_coherent'] = n_seed_coherent
            if n_seed_coherent is None:
                return PseudoFamilyVerdict(
                    is_pseudo=False,
                    reason=("zero qualifying clusters and incoherent copies, but the "
                            "seed→copy blastn could not run — cannot rule out real signal, "
                            "kept (§4.3 when in doubt, keep)"),
                    metrics=metrics)
            if n_seed_coherent >= seed_coherence_min:
                return PseudoFamilyVerdict(
                    is_pseudo=False,
                    reason=(f"zero qualifying clusters and copies mutually incoherent, but "
                            f"the family consensus strongly matches {n_seed_coherent} of its "
                            f"copies (>= {seed_coherence_min}) — real family with a divergent / "
                            f"anomalous copy set, kept (§4.3 存疑即留)"),
                    metrics=metrics)

        return PseudoFamilyVerdict(
            is_pseudo=True,
            reason=(f"PSEUDO: zero qualifying clusters, incoherent copy set, and the seed "
                    f"recruits no coherent copy — median pairwise identity "
                    f"{median_identity:.3f} < {max_intra_identity} "
                    f"and only {homologous_pair_frac:.1%} of {stats['n_pairs']} copy pairs "
                    f"share homology (< {min_homologous_frac:.0%}); the {n_copies} copies "
                    f"are a mutually-unrelated grab-bag with no recoverable subfamily"),
            metrics=metrics)

    # ── Coherent enough on at least one axis → KEEP (when in doubt, keep) ────
    held = []
    if not incoherent_identity:
        held.append(f"median identity {median_identity:.3f} >= {max_intra_identity}")
    if not unconnected:
        held.append(f"homologous-pair fraction {homologous_pair_frac:.2f} "
                    f">= {min_homologous_frac}")
    return PseudoFamilyVerdict(
        is_pseudo=False,
        reason=("zero qualifying clusters but copies still coherent ("
                + "; ".join(held) + ") — kept (§4.3 when in doubt, keep)"),
        metrics=metrics)
