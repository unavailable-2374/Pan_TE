"""
Phase 1 — N2: subfamily clustering (the "always cluster" module).

This is the structural inversion the v2 redesign turns on (REFINE_STRATEGY_DESIGN_v2.md
§3, §1.1/§1.2): the per-family body changes from "extract → one consensus → chimera"
to "extract → CLUSTER → one consensus PER qualifying cluster". A homogeneous family
emits exactly one cluster (and, per §3.5 naming, keeps its original `mdl_R<N>` id and
reproduces today's single-consensus output byte-for-byte); a divergent family emits
several subfamily clusters (`mdl_R<N>_sf<k>`). There is no "decide whether to split"
branch — the consensus count is an emergent property of the divergence structure
(§3.6).

What N2 implements (BED copies only; recall stays OFF — N4 turns it on):
  * pairwise_distance — copy×copy distance matrix from all-vs-all BLASTN LOCAL-alignment
    identity among a family's OWN copies (§3.2 "Algorithm", the proven precedent in the
    old Refiner/phase3_consensus_building.py: BLASTN recruit → scipy hierarchical). Each
    pair's distance is 1 − (best-HSP %identity × min(1, coverage-of-shorter)); a pair
    with no significant HSP gets distance 1.0. This is local, so it is robust to the
    rearrangement / internal-repeat content that makes mdl-stage copies mis-register
    under a single global MSA (measured below), and it is the pairwise copy×copy
    substrate §3.2 requires — NOT the 1-D BED divergence-to-seed proxy (two copies
    equidistant from the seed can be far from each other). The all-vs-all BLAST runs on a
    per-family DB of the family's own copies only — O(copies²), never the genome.
  * cluster_copies — scipy average-linkage cut at `subfamily_divergence_cut`, with a size
    route to a cheap k-mer graph (connected-component) regime for very large copy sets,
    mirroring M2's MAFFT↔abPOA size split (§3.2 "Algorithm").
  * qualifying_clusters — the conservative gate (§3.3): member≥min, intra-cluster
    homogeneity≥floor, inter-cluster distinctness≥cut. Sub-threshold clusters'
    members are REASSIGNED to the nearest qualifying cluster (never discarded); a
    family too sparse to form even one qualifying cluster falls back to a single
    all-copy cluster (the M2 never-regress invariant, §3.3 last paragraph).
  * subfamily lineage record (§3.5) — parent R, subfamily index, member count, intra-
    cluster mean identity, the cut used. Written into the record dict only (like
    `consensus_source` / `is_chimera_fragment`); NEVER into a FASTA header.

Why BLASTN local identity and not k-mer / global-MSA — MEASURED on real chr4 (§3.2):
Two cheaper substrates were tried on the real chr4 BED copy set and BOTH fail to resolve
subfamilies, for data reasons intrinsic to mdl-stage copies:
  (1) exact-k-mer (canonical-13-mer) Jaccard SATURATES at ~0.98–1.0 for homogeneous and
      divergent families alike (≈0.8^13 ≈ 5% of 13-mers survive at 80% identity, indels
      shatter more; and M1's per-copy-unique padding flank adds non-homologous content);
  (2) a single global MSA (MAFFT/abPOA) MIS-REGISTERS these copies: a 689 bp homogeneous
      family (R=149) aligned to a 6303-column MSA (≈10× the element) with median pairwise
      identity ~0.50 — because the copies carry rearrangements / internal repeats that a
      colinear MSA cannot stack, so "homogeneous" and "divergent" families looked
      identical (~0.5) under MSA-column identity. Removing the padding (pad=0) did NOT fix
      it — the misregistration is in the copies, not the flank.
BLASTN local alignment sidesteps both: it scores the best LOCAL homologous span per pair,
so rearranged copies still register on their shared blocks. MEASURED on real chr4 at
cut 0.30 (≥70% local identity within a subfamily):
  - genuinely divergent multi-copy families RESOLVE: R=8 (40 copies) → 4 subfamilies
    (intra-identity 0.82–0.92, inter-subfamily identity 0.18); R=17 → 3 (intra 0.80–0.92,
    inter 0.08); R=40 → 2 (intra 0.79–0.82, inter 0.65);
  - homogeneous / low-copy families STAY one cluster (R=113/89/94/149/216 → 1), so the
    degrade-to-today invariant holds and no over-splitting occurs.
See tests/integration_n2_cluster.py for the full real numbers reproduced on every run.

The k-mer path is RETAINED only as the size-routed coarse fallback for copy sets too
large for the O(c²) all-vs-all BLAST (cluster_graph_size_threshold) — such a set is in
practice one young high-copy family that should stay one cluster, and the k-mer
saturation degrades safely to exactly that via the never-regress fallback. No path here
is O(families × genome): every BLAST is a family's own copies against a tiny per-family
DB.

§1.2 respected: BLASTN here only MEASURES pairwise distance. Every qualifying cluster's
CONSENSUS is still built by a SEPARATE per-cluster M2 call (extend-align-trim) — v2
opposes "one consensus from a mixed MSA", not "a local-alignment ruler over copies".

See REFINE_STRATEGY_DESIGN_v2.md §3 (all), §1.1/§1.2, §3.2 "Algorithm", §7, §8.2, §10 N2.
"""

import logging
import os
import random
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np

logger = logging.getLogger(__name__)


# ═══════════════════════════════════════════════════════════════════════
# Cluster structure
# ═══════════════════════════════════════════════════════════════════════

@dataclass
class Cluster:
    """One subfamily cluster: its member copies plus intra-cluster statistics.

    Attributes
    ----------
    members : List[Copy]
        The phase1_extract.Copy objects assigned to this cluster (the input to M2's
        per-cluster `refine_family`).
    mean_identity : float
        Mean pairwise identity (1 - distance) within the cluster, over the BLASTN
        local-identity distance matrix. 1.0 for a singleton (trivially self-identical).
        This is the homogeneity score the qualifying gate reads and the lineage stores.
    member_indices : List[int]
        Indices of the members in the ORIGINAL copies list passed to cluster_copies —
        lets callers/tests verify member conservation without identity-matching the
        Copy objects.
    """
    members: List = field(default_factory=list)
    mean_identity: float = 1.0
    member_indices: List[int] = field(default_factory=list)

    @property
    def size(self) -> int:
        return len(self.members)


# ═══════════════════════════════════════════════════════════════════════
# Coarse k-mer (Jaccard) distance — RETAINED as the size-routed fallback only
# ═══════════════════════════════════════════════════════════════════════
# Primary distance is the BLASTN local-identity substrate below; this k-mer path is
# kept ONLY for copy sets too large for an O(c²) all-vs-all BLAST (the graph regime),
# where its saturation is acceptable because such a set is, in practice, one young
# high-copy family that should stay one cluster anyway. See module docstring.

def _kmer_set(seq: str, k: int) -> frozenset:
    """Set of canonical k-mers (min of fwd / revcomp) over the ACGT alphabet.

    Canonicalizing by strand is belt-and-suspenders here — copies arriving from
    phase1_extract are already strand-co-oriented — but it makes the distance robust
    if an upstream orientation ever slips, and it costs nothing. k-mers containing a
    non-ACGT base (N, gap) are skipped: they carry no comparative signal and would
    only inflate the union."""
    s = seq.upper()
    out = set()
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    n = len(s)
    if n < k:
        return frozenset(out)
    for i in range(n - k + 1):
        kmer = s[i:i + k]
        if not all(b in comp for b in kmer):
            continue
        rc = ''.join(comp[b] for b in reversed(kmer))
        out.add(kmer if kmer <= rc else rc)
    return frozenset(out)


def _jaccard_distance(a: frozenset, b: frozenset) -> float:
    """1 - |a∩b| / |a∪b|. Two empty sets -> distance 0 (both featureless, treat as
    identical so they fold together rather than spuriously splitting)."""
    if not a and not b:
        return 0.0
    if not a or not b:
        return 1.0
    inter = len(a & b)
    union = len(a) + len(b) - inter
    return 1.0 - (inter / union) if union else 0.0


def _kmer_pairwise_distance(copies: List, config) -> np.ndarray:
    """Copy×copy distance matrix via canonical-k-mer Jaccard distance (fallback only).

    O(c² · |kmers|) set operations, no alignment. SATURATES near 1.0 on real mdl-stage
    copies (module docstring), so it is NOT the primary substrate — it is the size-routed
    fallback for copy sets too large for the O(c²) all-vs-all BLAST, where its
    saturation degrades safely to one cluster. Symmetric (n×n), zero diagonal; n<=1 -> a
    trivial zero matrix."""
    n = len(copies)
    dist = np.zeros((n, n), dtype=float)
    if n <= 1:
        return dist
    k = getattr(config, 'subfamily_kmer_size', 13)
    sketches = [_kmer_set(c.sequence, k) for c in copies]
    for i in range(n):
        for j in range(i + 1, n):
            d = _jaccard_distance(sketches[i], sketches[j])
            dist[i, j] = d
            dist[j, i] = d
    return dist


# ═══════════════════════════════════════════════════════════════════════
# Primary distance: all-vs-all BLASTN local-identity (§3.2 "Algorithm")
# ═══════════════════════════════════════════════════════════════════════

def _blastn_identity_distance(copies: List, config) -> Optional[np.ndarray]:
    """Copy×copy distance from all-vs-all BLASTN LOCAL-alignment identity.

    distance(i, j) = 1 − (best-HSP %identity / 100), counting only HSPs that cover at
    least `subfamily_min_hsp_coverage` of the SHORTER copy. The coverage GATE (not a
    coverage WEIGHT) is the key choice for mdl-stage copies: they arrive PADDED from M1
    with per-copy-unique flank, and a coverage-weighted score (pident × coverage) is
    deflated by that flank so even a homogeneous family looks divergent — MEASURED on
    real chr4 (R=8 padded: coverage-weighted gave intra-identity 0.25–0.36, garbage). A
    coverage GATE instead rejects only short spurious HSPs, then uses the LOCAL %identity
    of the real homologous span, which is flank-length-invariant. MEASURED: this resolves
    R=8 → 4 subfamilies at intra-identity 0.82–0.90 and keeps homogeneous families at 1,
    on the same PADDED copies the production path passes in. A pair with no HSP clearing
    the gate is left at distance 1.0 (maximally distant — cannot join a cluster).
    Symmetric, zero diagonal.

    All-vs-all is run on a per-family BLAST DB of the family's OWN copies (O(copies²)
    HSPs, never the genome). Returns None on any tool failure (missing binary, non-zero
    exit) so the caller degrades to the never-regress single-cluster fallback rather than
    fabricating a matrix."""
    n = len(copies)
    if n <= 1:
        return np.zeros((n, n), dtype=float)

    makeblastdb = getattr(config, 'makeblastdb_exe', 'makeblastdb') or 'makeblastdb'
    blastn = getattr(config, 'blastn_exe', 'blastn') or 'blastn'
    evalue = str(getattr(config, 'blastn_evalue', 1e-5))
    min_cov = getattr(config, 'subfamily_min_hsp_coverage', 0.4)

    tmp_dir = tempfile.mkdtemp(prefix='n2_blast_')
    fa = os.path.join(tmp_dir, 'copies.fa')
    db = os.path.join(tmp_dir, 'db')
    try:
        with open(fa, 'w') as fh:
            for i, c in enumerate(copies):
                # Positional names s0..s(n-1): map HSP rows back to copy index by name,
                # independent of any Copy.id collisions.
                fh.write(f">s{i}\n{c.sequence}\n")

        mk = subprocess.run([makeblastdb, '-in', fa, '-dbtype', 'nucl', '-out', db],
                            capture_output=True, text=True, timeout=600)
        if mk.returncode != 0:
            logger.warning("pairwise_distance(blastn): makeblastdb exit %d: %s",
                           mk.returncode, (mk.stderr or '')[:200])
            return None

        # word_size 7 + dust off: sensitive enough for the divergent (~70%) subfamily
        # band; the M1 copies are already low-complexity-filtered upstream so dust-off is
        # safe and avoids masking real TE k-mers.
        cmd = [blastn, '-query', fa, '-db', db,
               '-outfmt', '6 qseqid sseqid pident length qlen slen',
               '-evalue', evalue, '-word_size', '7', '-dust', 'no']
        br = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
        if br.returncode != 0:
            logger.warning("pairwise_distance(blastn): blastn exit %d: %s",
                           br.returncode, (br.stderr or '')[:200])
            return None

        # Best gate-passing HSP %identity (as a fraction) per unordered pair.
        best: Dict[Tuple[int, int], float] = {}
        for line in br.stdout.splitlines():
            parts = line.split('\t')
            if len(parts) < 6:
                continue
            q, s, pid, ln, ql, sl = parts
            if q == s:
                continue
            try:
                i = int(q[1:]); j = int(s[1:])
                pident = float(pid); aln = int(ln)
                qlen = int(ql); slen = int(sl)
            except ValueError:
                continue
            if i >= n or j >= n:
                continue
            denom = min(qlen, slen)
            cov = (aln / denom) if denom else 0.0
            if cov < min_cov:
                continue                       # short spurious HSP -> not evidence
            ident = pident / 100.0
            key = (i, j) if i < j else (j, i)
            if ident > best.get(key, 0.0):
                best[key] = ident

        dist = np.ones((n, n), dtype=float)
        np.fill_diagonal(dist, 0.0)
        for (i, j), ident in best.items():
            d = 1.0 - ident
            dist[i, j] = d
            dist[j, i] = d
        return dist
    except FileNotFoundError as e:
        logger.warning("pairwise_distance(blastn): binary not found (%s)", e)
        return None
    except Exception as e:  # noqa: BLE001 — report + fall back, never fabricate a matrix
        logger.warning("pairwise_distance(blastn): failed (%s)", e)
        return None
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


# Single-entry memo so the repeated pairwise_distance(copies, config) calls within one
# family (cluster_copies once, qualifying_clusters up to 3×) run the all-vs-all BLAST
# ONCE. Keyed on the identity-tuple of the copy set (Copy.id order); a new family's first
# call overwrites the previous entry, so the cache never grows and stays correct across
# shard-parallel workers (each worker is its own process with its own module instance).
_LAST_DIST_KEY: Optional[Tuple] = None
_LAST_DIST_VAL: Optional[np.ndarray] = None


def pairwise_distance(copies: List, config) -> np.ndarray:
    """Copy×copy distance matrix — PRIMARY substrate: all-vs-all BLASTN local identity.

    Memoized on the copy set so the per-family call chain (cluster_copies →
    qualifying_clusters) runs the all-vs-all BLAST exactly once.

    distance(i, j) = 1 − (best gate-passing HSP %identity); a coverage GATE (not weight)
    rejects short spurious HSPs so the score is flank-length-invariant local identity.
    Local alignment is robust to the rearrangement / internal-repeat content that makes
    mdl-stage copies mis-register under a single global MSA (measured: a 689 bp
    homogeneous family aligns to a 6303-column MSA at ~0.5 pairwise identity, no
    separation; BLASTN local identity instead resolves R=8 → 4 subfamilies at intra-
    identity 0.82–0.90 and keeps homogeneous families at 1 — see module docstring +
    tests/integration_n2_cluster.py).

    Size route (mirrors M2's MAFFT↔abPOA split, §3.2 "Algorithm"): copy sets larger than
    `cluster_graph_size_threshold` skip the O(c²) all-vs-all BLAST and fall back to the
    coarse k-mer distance — an explicit, documented degradation (such a set is in practice
    one young high-copy family that should stay one cluster, and its k-mer distance
    saturating to ~1.0 yields exactly that via the never-regress fallback). This keeps the
    cost at one per-family BLAST and never an O(families × genome) sweep.

    Never-regress fallback (degrade to today, §3.3): if BLAST cannot run (binary missing,
    non-zero exit), return an all-zero matrix — every pair at distance 0 ⇒ one cluster ⇒
    today's single all-copy consensus. Safe, never fabricated.

    Returns a symmetric (n×n) float matrix with a zero diagonal. n<=1 -> a trivial
    (n×n) zero matrix.
    """
    global _LAST_DIST_KEY, _LAST_DIST_VAL
    key = (tuple(c.id for c in copies),
           getattr(config, 'subfamily_kmer_size', 13),
           getattr(config, 'cluster_graph_size_threshold', 500),
           getattr(config, 'blastn_evalue', 1e-5),
           getattr(config, 'subfamily_min_hsp_coverage', 0.4),
           getattr(config, 'distance_blast_max_total_bp', 500_000))
    if key == _LAST_DIST_KEY and _LAST_DIST_VAL is not None \
            and _LAST_DIST_VAL.shape == (len(copies), len(copies)):
        return _LAST_DIST_VAL
    val = _compute_pairwise_distance(copies, config)
    _LAST_DIST_KEY, _LAST_DIST_VAL = key, val
    return val


def _bp_bounded_subsample(copies: List, max_bp: int, seed: int) -> List[int]:
    """Deterministic indices of a copy subset whose total bp stays <= `max_bp`.

    All-vs-all BLAST cost grows with total bp, so a few-copy giant-element family (e.g. a
    65-copy 20 kb LTR) would otherwise run a multi-minute BLAST. We bound it the way M2
    bounds its MSA: a divergence-stratified deterministic subset (G10 posture). Returns
    the ORIGINAL indices kept, in ascending order. Always keeps at least 2 (so a pair
    distance is still defined). Copies are taken stratified by divergence ascending so the
    subset spans the family's age range rather than only its youngest members."""
    order = sorted(range(len(copies)), key=lambda i: copies[i].divergence)
    rng = random.Random(seed)
    rng.shuffle(order)  # break ties deterministically without biasing toward file order
    order.sort(key=lambda i: copies[i].divergence)
    kept: List[int] = []
    total = 0
    for i in order:
        L = len(copies[i].sequence)
        if kept and total + L > max_bp:
            continue
        kept.append(i)
        total += L
    if len(kept) < 2:                       # guarantee a defined pairwise comparison
        kept = order[:2]
    return sorted(kept)


def prefilter_copies_for_clustering(copies: List, config, ref_len: int = 0) -> List:
    """Bound clustering cost BEFORE the O(copies²) all-vs-all distance.

    Two scalability hazards measured on real full-Arabidopsis data (family R=1567:
    927 copies of an 8.5 kb element, one mis-annotated 'copy' 671 kb long):
      1. all-vs-all distance is O(count²) — the existing bp cap does not bind when copies
         are individually small, so 927 copies → ~860 k pairs → ~210 s for ONE family
         (and the >threshold k-mer fallback is also O(count²)). A hard COUNT cap is the
         real knob.
      2. a single giant mis-annotated instance (≫ the element length) bloats every
         pairwise comparison and is noise, not a real copy.

    So: drop length-outlier copies (> `cluster_outlier_len_mult` × median copy length),
    then, if still over `cluster_distance_max_copies`, keep a divergence-stratified,
    deterministic representative subsample spanning the family's age range. Subfamily
    structure and per-cluster consensi are well captured by the subsample (mirrors M2's
    MSA subsample cap); held-out copies are NOT records and are not lost — the family
    still emits its subfamily consensi. Called once per family so cluster_copies and
    qualifying_clusters operate on the SAME index space.
    """
    n = len(copies)
    if n <= 2:
        return list(copies)

    # 1. length-outlier exclusion. A copy is anomalous only if it is much LONGER than the
    # ELEMENT (consensus/seed length) — e.g. the 671 kb mis-annotated 'copy' of an 8.5 kb
    # family. Copies SHORTER than the element are legitimate truncated fragments and must
    # be kept (a fragment-heavy family's median copy length is tiny, so a median-relative
    # cutoff would wrongly discard the real full-length copies). ref_len is the family
    # consensus length; fall back to the 90th-percentile copy length if not provided.
    mult = getattr(config, 'cluster_outlier_len_mult', 5.0)
    if mult and mult > 0:
        if ref_len and ref_len > 0:
            ref = ref_len
        else:
            lens_sorted = sorted(len(c.sequence) for c in copies)
            ref = lens_sorted[min(len(lens_sorted) - 1, int(len(lens_sorted) * 0.9))]
        if ref > 0:
            kept = [c for c in copies if len(c.sequence) <= mult * ref]
            if len(kept) >= 2:
                copies = kept

    # 2. count cap (divergence-stratified, deterministic, spans the age range)
    cap = getattr(config, 'cluster_distance_max_copies', 150)
    if cap and len(copies) > cap:
        order = sorted(range(len(copies)), key=lambda i: copies[i].divergence)
        rng = random.Random(getattr(config, 'subsample_seed', 42))
        rng.shuffle(order)                          # deterministic tie-break
        order.sort(key=lambda i: copies[i].divergence)
        step = len(order) / cap
        pick = sorted({order[min(len(order) - 1, int(k * step))] for k in range(cap)})
        copies = [copies[i] for i in pick]

    return copies


def _compute_pairwise_distance(copies: List, config) -> np.ndarray:
    """Uncached body of pairwise_distance (BLASTN local identity + size fallbacks)."""
    n = len(copies)
    if n <= 1:
        return np.zeros((n, n), dtype=float)

    graph_thresh = getattr(config, 'cluster_graph_size_threshold', 500)
    if n > graph_thresh:
        logger.debug("pairwise_distance: %d copies > %d -> coarse k-mer fallback",
                     n, graph_thresh)
        return _kmer_pairwise_distance(copies, config)

    # Bound the all-vs-all BLAST by total bp (M2-style cap): a giant-element family is
    # subsampled for the distance, and unsampled copies are left at distance 0 to every
    # other copy so they fold into the main cluster (never-regress) rather than spuriously
    # splitting off — a divergent giant family is rare and N4's recall will revisit it.
    max_bp = getattr(config, 'distance_blast_max_total_bp', 500_000)
    total_bp = sum(len(c.sequence) for c in copies)
    if total_bp > max_bp:
        seed = getattr(config, 'subsample_seed', 42)
        keep = _bp_bounded_subsample(copies, max_bp, seed)
        logger.debug("pairwise_distance: total %d bp > %d -> BLAST on %d/%d copies",
                     total_bp, max_bp, len(keep), n)
        sub_copies = [copies[i] for i in keep]
        sub = _blastn_identity_distance(sub_copies, config)
        if sub is None:
            return np.zeros((n, n), dtype=float)
        full = np.zeros((n, n), dtype=float)
        for a, ia in enumerate(keep):
            for b in range(a + 1, len(keep)):
                ib = keep[b]
                full[ia, ib] = sub[a, b]
                full[ib, ia] = sub[a, b]
        return full

    dist = _blastn_identity_distance(copies, config)
    if dist is None:
        logger.debug("pairwise_distance: BLAST unavailable (%d copies) -> "
                     "single-cluster fallback", n)
        return np.zeros((n, n), dtype=float)
    return dist


# ═══════════════════════════════════════════════════════════════════════
# Hierarchical / graph clustering
# ═══════════════════════════════════════════════════════════════════════

def _mean_intra_identity(dist: np.ndarray, idxs: List[int]) -> float:
    """Mean pairwise identity (1 - distance) over the members in `idxs`.

    A singleton (or empty) cluster is trivially self-identical -> 1.0."""
    m = len(idxs)
    if m <= 1:
        return 1.0
    tot = 0.0
    cnt = 0
    for a in range(m):
        for b in range(a + 1, m):
            tot += 1.0 - dist[idxs[a], idxs[b]]
            cnt += 1
    return tot / cnt if cnt else 1.0


def _hierarchical_labels(dist: np.ndarray, cut: float) -> np.ndarray:
    """scipy average-linkage labels at a divergence `cut` (condensed-matrix path).

    Average linkage on the BLASTN local-identity distance is the natural fit for
    divergence-cut subfamily clustering and is exactly the pattern the old
    Refiner/phase3_consensus_building.py used (BLASTN recruit → scipy hierarchical →
    tiered consensus); v2 adapts it (§3.2 "Algorithm"). fcluster with
    criterion='distance' cuts the dendrogram so every merged cluster has cophenetic
    distance <= cut — i.e. members within `cut` divergence stay together, sub-lineages
    beyond `cut` separate."""
    from scipy.cluster.hierarchy import fcluster, linkage
    from scipy.spatial.distance import squareform

    n = dist.shape[0]
    if n <= 1:
        return np.zeros(n, dtype=int)
    # squareform needs an exactly-symmetric, zero-diagonal matrix; enforce it
    # (float round-trips through Jaccard are symmetric, but guard against drift).
    sym = (dist + dist.T) / 2.0
    np.fill_diagonal(sym, 0.0)
    condensed = squareform(sym, checks=False)
    Z = linkage(condensed, method='average')
    return fcluster(Z, t=cut, criterion='distance')


def _graph_labels(dist: np.ndarray, cut: float) -> np.ndarray:
    """Connected-component labels on the threshold graph (edge iff distance <= cut).

    The cheap graph regime for very large copy sets (§3.2 "graph community detection …
    the alternative for very large copy sets; the choice is by copy-set size,
    mirroring M2's MAFFT-vs-abPOA regime split"). A full Leiden/CPM pass (as TE-looker
    uses) is the heavier option; at refine scale a single-linkage-style connected-
    component cut on the same `cut` threshold is the size-appropriate, dependency-free
    realization and keeps the cut threshold's meaning identical to the hierarchical
    path. Implemented with a union-find over the thresholded adjacency."""
    n = dist.shape[0]
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[max(ra, rb)] = min(ra, rb)

    for i in range(n):
        for j in range(i + 1, n):
            if dist[i, j] <= cut:
                union(i, j)
    # Relabel roots to contiguous 1..K (fcluster convention: labels start at 1).
    roots = {}
    labels = np.zeros(n, dtype=int)
    nxt = 1
    for i in range(n):
        r = find(i)
        if r not in roots:
            roots[r] = nxt
            nxt += 1
        labels[i] = roots[r]
    return labels


def cluster_copies(copies: List, config) -> List[Cluster]:
    """Partition an assembled copy set into raw subfamily clusters (always runs).

    Size route (§3.2, mirrors M2's MAFFT↔abPOA split): copy sets up to
    `cluster_graph_size_threshold` take the scipy hierarchical (average-linkage)
    path; larger sets take the cheap connected-component graph regime so the O(c²)
    linkage memory does not blow up on ultra-abundant families. Both cut at the same
    `subfamily_divergence_cut`, so the threshold's biological meaning (subfamily-level
    divergence) is identical across the route.

    Returns a list of `Cluster` (raw, pre-gate). A 0- or 1-copy input returns a single
    cluster holding whatever was passed (the gate / caller decides if M2 runs). Cluster
    order is deterministic: sorted by size descending, ties by smallest member index,
    so naming `_sf1.._sfk` is stable across runs (G10 posture).
    """
    n = len(copies)
    if n == 0:
        return []
    if n == 1:
        return [Cluster(members=list(copies), mean_identity=1.0, member_indices=[0])]

    cut = config.subfamily_divergence_cut
    dist = pairwise_distance(copies, config)

    graph_thresh = getattr(config, 'cluster_graph_size_threshold', 500)
    if n <= graph_thresh:
        labels = _hierarchical_labels(dist, cut)
    else:
        logger.debug("cluster_copies: %d copies > %d -> graph regime",
                     n, graph_thresh)
        labels = _graph_labels(dist, cut)

    by_label: Dict[int, List[int]] = {}
    for idx, lab in enumerate(labels):
        by_label.setdefault(int(lab), []).append(idx)

    clusters: List[Cluster] = []
    for _lab, idxs in by_label.items():
        clusters.append(Cluster(
            members=[copies[i] for i in idxs],
            mean_identity=_mean_intra_identity(dist, idxs),
            member_indices=list(idxs)))

    # Deterministic order: larger clusters first; ties by smallest member index.
    clusters.sort(key=lambda c: (-c.size, min(c.member_indices)))
    return clusters


# ═══════════════════════════════════════════════════════════════════════
# Qualifying-cluster gate + member reassignment (never drop)
# ═══════════════════════════════════════════════════════════════════════

def _centroid_distance(dist: np.ndarray, a_idxs: List[int],
                       b_idxs: List[int]) -> float:
    """Mean cross-cluster distance between member sets a and b (average linkage).

    Used both for inter-cluster distinctness and for nearest-cluster reassignment, so
    a sub-threshold member rejoins the cluster it is genuinely closest to."""
    if not a_idxs or not b_idxs:
        return 1.0
    tot = 0.0
    for ai in a_idxs:
        for bi in b_idxs:
            tot += dist[ai, bi]
    return tot / (len(a_idxs) * len(b_idxs))


def qualifying_clusters(clusters: List[Cluster], copies: List,
                        config) -> Tuple[List[Cluster], Dict]:
    """Apply the conservative qualifying gate and reassign sub-threshold members.

    Gate (§3.3 — a cluster emits its own consensus only if ALL hold):
      * member count >= `subfamily_min_members`
      * intra-cluster homogeneity (mean pairwise identity) >= `subfamily_intra_homogeneity_floor`
      * (distinctness from siblings is handled by the linkage cut itself: every
        emitted cluster is already >= cut apart by construction, and a merge of two
        too-close qualifying clusters is prevented up front — see the distinctness
        merge below)

    Never-drop invariant (§3.3, §8.2):
      * members of a sub-threshold cluster are REASSIGNED to the nearest qualifying
        cluster (by mean cross-cluster BLASTN-identity distance) — never discarded.
      * if NO cluster qualifies, the whole family falls back to ONE all-copy cluster
        (today's single-consensus M2 path, the never-regress floor).

    Returns (qualifying_clusters, audit) where audit carries member-conservation
    bookkeeping (total members in vs out, reassigned count, fallback flag) so the M4
    counting cross-check (§8.2) is a real audit, not a tautology.
    """
    total_in = sum(c.size for c in clusters)
    n_copies = len(copies)

    min_members = config.subfamily_min_members
    homog_floor = config.subfamily_intra_homogeneity_floor

    qualifying = [c for c in clusters
                  if c.size >= min_members and c.mean_identity >= homog_floor]
    subthreshold = [c for c in clusters
                    if not (c.size >= min_members and c.mean_identity >= homog_floor)]

    # ── No qualifying cluster: never-regress fallback to one all-copy cluster ──
    if not qualifying:
        all_idx = list(range(n_copies))
        dist = pairwise_distance(copies, config) if n_copies > 1 \
            else np.zeros((n_copies, n_copies))
        fallback = Cluster(members=list(copies),
                           mean_identity=_mean_intra_identity(dist, all_idx),
                           member_indices=all_idx)
        audit = {
            'total_members_in': total_in,
            'total_members_out': fallback.size,
            'n_qualifying': 1,
            'n_reassigned': 0,
            'fallback_single_cluster': True,
        }
        return [fallback], audit

    # ── Reassign each sub-threshold member to the nearest qualifying cluster ──
    dist = pairwise_distance(copies, config)
    n_reassigned = 0
    for sub in subthreshold:
        for mi in sub.member_indices:
            # Nearest qualifying cluster by mean distance from this single member.
            best = min(qualifying,
                       key=lambda q: _centroid_distance(dist, [mi], q.member_indices))
            best.members.append(copies[mi])
            best.member_indices.append(mi)
            n_reassigned += 1

    # Recompute intra-cluster homogeneity after absorbing reassigned members.
    for q in qualifying:
        q.mean_identity = _mean_intra_identity(dist, q.member_indices)

    # Re-sort for stable `_sf<k>` naming after reassignment may have changed sizes.
    qualifying.sort(key=lambda c: (-c.size, min(c.member_indices)))

    total_out = sum(c.size for c in qualifying)
    audit = {
        'total_members_in': total_in,
        'total_members_out': total_out,
        'n_qualifying': len(qualifying),
        'n_reassigned': n_reassigned,
        'fallback_single_cluster': False,
    }
    return qualifying, audit


# ═══════════════════════════════════════════════════════════════════════
# Subfamily lineage record (§3.5) — biological provenance, never a header
# ═══════════════════════════════════════════════════════════════════════

def subfamily_lineage(parent_rec: Dict, sf_index: int, n_subfamilies: int,
                      cluster: Cluster, config) -> Dict:
    """Build the lineage provenance dict for one emitted subfamily (§3.5).

    Recorded on the output record alongside `consensus_source` / `is_chimera_fragment`
    — biological provenance ONLY, never written into a FASTA header (§3.5, §10 N2):

      parent_R           : the parent family's R number / id
      subfamily_index    : 1-based index within this family (1..n_subfamilies)
      n_subfamilies      : how many subfamilies this family emitted (1 == homogeneous)
      member_count       : this cluster's member count (also drives the record's
                           `copies` so Phase 2's per-subfamily copy gates act correctly,
                           §7)
      intra_identity     : mean pairwise k-mer identity within the cluster
      divergence_cut     : the cut threshold used (auditable, tunable)
    """
    return {
        'parent_R': parent_rec.get('R', parent_rec.get('id', '')),
        'subfamily_index': sf_index,
        'n_subfamilies': n_subfamilies,
        'member_count': cluster.size,
        'intra_identity': round(float(cluster.mean_identity), 4),
        'divergence_cut': config.subfamily_divergence_cut,
    }


def subfamily_id(parent_id: str, sf_index: int, n_subfamilies: int) -> str:
    """Naming convention (§3.5, §10 N2): single-cluster family keeps its ORIGINAL id
    verbatim (`mdl_R<N>`, no suffix — output byte-identical to today); a multi-cluster
    family suffixes `_sf<k>` (1-based). The `_chimfrag<j>` suffix composes orthogonally
    downstream (N3) for a subfamily that is also structurally split."""
    if n_subfamilies <= 1:
        return parent_id
    return f"{parent_id}_sf{sf_index}"
