"""
Phase 1 — N1: family-completeness PRE-SIGNALS (cheap trigger layer, no recall).

This module answers one question per family, cheaply: *does this family still
have un-recruited members?* — i.e. is it worth spending a (later, N4) targeted
recall on?  It NEVER drops a family and NEVER touches the genome; every signal
here is O(instances) or O(1) over data already in hand (BED instance count,
per-instance BED divergence, consensus length).  The flag only decides where
expensive recall compute will be spent in N4; "selective" governs *where*, never
*whether* completeness is the target (REFINE_STRATEGY_DESIGN_v2.md §2.3).

What is implemented here (N1 cheap signals):
  * bed_divergence_spread  — high spread ⇒ family already straddles a wide
    divergence band ⇒ a divergent tail below mdl-repeat's ~70% acceptance likely
    exists (the R=150 case: instance divs {0.02, 0.676}).
  * copy_count_vs_length_class — a long element with very few BED instances is
    under-instanced relative to a length-class prior (mdl-repeat under-counted it).
  * under_instanced — BED instances < min_copies_for_msa (the M5 trigger, now
    demoted from the *only* trigger to one signal among several).

N4 resolution of the two signals N1 deferred:
  * sketch pre-screen hit-count beyond BED (§2.3 signal 3) — EMPIRICALLY DISABLED.
    A real-chr4 test (sourmash MinHash sketch of the 472 bp R=150 seed vs the
    18.5 Mbp chr4 genome sketch) returned ZERO matches under both `sourmash gather`
    and `sourmash search --containment` even at scaled=1 (all k-mers kept): the
    surviving shared 21-mers per divergent partial are too sparse (~1-13% of the
    446 seed k-mers) for a genome-scale containment estimate to clear any usable
    threshold.  This is the N2 k-mer-saturation result generalized — at 75-90%
    identity the sketch is blind to exactly the divergent tail it would gate on.
    So the SELECTIVE recall gate is the N1 specific-signal co-fire below, NOT a
    sketch.  See config.enable_genome_sketch_prescreen (default False) for the
    full record.
  * boundary-truncation signal from a first-pass MSA (§2.3 signal 4) — still
    deferred (needs a first-pass MSA occupancy pass; not required for the N4
    recall gate, which is satisfied by the cheap signals).

N4 adds `recall_eligible` (below): a SELECTIVE gate that, unlike
`assess_completeness`'s "any signal fires", requires `under_instanced` to CO-FIRE
with a specific signal (bed_divergence_spread or copy_count_vs_length_class).  On
real chr4 this fires for 20.9% of families vs the 85.2% that `under_instanced`
alone would — `under_instanced` is the norm in an 84%-low-copy genome and cannot
by itself justify the O(genome) Tier-2 blastn.

The signal functions take a `List[Instance]` (phase1_extract.Instance) and a
`record` dict (the Phase-0 record carrying at least `id`, `actual_length`/`length`,
`copies`).  `assess_completeness` combines them into a CompletenessVerdict.
"""

import logging
from dataclasses import dataclass, field
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)

# Import Instance for type clarity; phase1_extract is a sibling module.
try:
    from phase1_extract import Instance  # noqa: F401
except Exception:  # pragma: no cover - import-path robustness for test harnesses
    Instance = object  # type: ignore


# ═══════════════════════════════════════════════════════════════════════
# Verdict structure
# ═══════════════════════════════════════════════════════════════════════

@dataclass
class CompletenessVerdict:
    """Per-family completeness pre-signal verdict (N1).

    Attributes
    ----------
    record_id : str
        The family/record this verdict is for.
    incomplete_recall_eligible : bool
        True iff at least one cheap signal fired ⇒ this family is worth a later
        (N4) targeted recall.  A family that is NOT eligible is declared
        BED-complete for N1 purposes and proceeds to clustering on its BED copies
        alone.  This flag drops no family either way.
    flags : Dict[str, bool]
        Which individual signals fired (keys: 'bed_divergence_spread',
        'copy_count_vs_length_class', 'under_instanced').
    reasons : List[str]
        Human-readable one-line justification per fired signal (for the metric
        output / audit; never enters FASTA headers).
    n_instances : int
        BED instance count used (echoed for downstream bookkeeping).
    bed_divergence_spread : Optional[float]
        Population std of per-instance BED divergence (None when n < 1).
    """
    record_id: str
    incomplete_recall_eligible: bool = False
    flags: Dict[str, bool] = field(default_factory=dict)
    reasons: List[str] = field(default_factory=list)
    n_instances: int = 0
    bed_divergence_spread: Optional[float] = None


# ═══════════════════════════════════════════════════════════════════════
# Cheap signal functions (O(instances) / O(1), zero genome alignment)
# ═══════════════════════════════════════════════════════════════════════

def _divergences(instances) -> List[float]:
    return [getattr(i, 'divergence', 0.0) for i in instances]


def bed_divergence_spread(instances) -> Optional[float]:
    """Population std of per-instance BED divergence.

    High spread means the family already spans a wide divergence band, which
    indicates a divergent tail likely exists below mdl-repeat's acceptance floor
    (§2.3 signal 1).  Returns None for <1 instance, 0.0 for a single instance
    (a single point has zero spread but is trivially defined), and the population
    std otherwise.  O(n), no alignment.
    """
    divs = _divergences(instances)
    n = len(divs)
    if n == 0:
        return None
    if n == 1:
        return 0.0
    mean = sum(divs) / n
    var = sum((d - mean) ** 2 for d in divs) / n  # population variance
    return var ** 0.5


def copy_count_vs_length_class(record: Dict, config) -> bool:
    """True iff the family is under-instanced for its length class (§2.3 signal 2).

    mdl-repeat tends to under-count longer elements (a full-length element leaves
    fewer near-identical full copies than a short high-copy SINE/MITE), so a long
    element with very few BED instances signals missing members.  The expectation
    table (config.length_class_copy_expectation) maps length classes to an
    expected minimum copy count; the family is flagged when its BED instance
    count (preferred, the real evidence) — or `copies` if no instance count is
    available — is below that expectation.  O(1), no alignment.
    """
    table = getattr(config, 'length_class_copy_expectation', ())
    if not table:
        return False
    length = record.get('actual_length') or record.get('length') or 0
    # Prefer the actual BED instance count if the record carries one; fall back
    # to the header `copies` (mdl-repeat's accepted-copy count).
    n = record.get('n_instances')
    if n is None:
        n = record.get('copies', 0)
    expect = _expected_copies(length, table)
    return n < expect


def _expected_copies(length: int, table) -> int:
    """First (descending-floor) length-class entry whose floor <= length."""
    best_floor = -1
    best_expect = 0
    for floor, expect in table:
        if length >= floor and floor > best_floor:
            best_floor = floor
            best_expect = expect
    return best_expect


def under_instanced(record: Dict, config) -> bool:
    """True iff BED instance count < config.min_copies_for_msa (§2.3 signal 5).

    The implemented M5 trigger, demoted to one signal among several.  Prefers the
    real BED instance count; falls back to `copies`.  O(1).
    """
    n = record.get('n_instances')
    if n is None:
        n = record.get('copies', 0)
    return n < config.min_copies_for_msa


# ═══════════════════════════════════════════════════════════════════════
# Combined verdict
# ═══════════════════════════════════════════════════════════════════════

def assess_completeness(record: Dict, instances, config) -> CompletenessVerdict:
    """Combine the N1 cheap signals into a recall-eligibility verdict.

    Parameters
    ----------
    record : Dict
        Phase-0 record dict.  Must carry `id`; `actual_length`/`length`, `copies`,
        and (optionally) `n_instances` are used by the length-class / under-
        instanced signals.  This function will populate `n_instances` from the
        passed `instances` so callers need not pre-fill it.
    instances : List[Instance]
        The family's BED instances (phase1_extract.Instance).  May be empty.
    config : RefinerMdlConfig

    Returns
    -------
    CompletenessVerdict
        `incomplete_recall_eligible` is True iff any cheap signal fired.  Nothing
        is dropped; the verdict only decides whether a later (N4) recall is worth
        running for this family.  Decision cost is O(len(instances)) + O(1), with
        no genome alignment.

    NOTE(N4): two further §2.3 signals — a sourmash/minimizer sketch pre-screen
    hit count beyond the BED, and an M2 first-pass boundary-truncation signal —
    are intentionally NOT consulted here; they require the shared genome sketch
    index and a first-pass MSA that land in N4.  When N4 wires them in they should
    be added as additional `flags`/`reasons` here without changing the N1 logic.
    """
    n = len(instances)
    # Populate n_instances on a shallow copy of the relevant fields so the O(1)
    # signal functions see the real BED count rather than the header `copies`.
    rec_view = dict(record)
    rec_view['n_instances'] = n

    flags: Dict[str, bool] = {}
    reasons: List[str] = []

    # Signal 1: BED divergence spread.
    spread = bed_divergence_spread(instances)
    spread_fired = False
    if (spread is not None
            and n >= getattr(config, 'bed_spread_min_instances', 2)
            and spread >= config.bed_divergence_spread_trigger):
        spread_fired = True
        reasons.append(
            f"bed_divergence_spread={spread:.3f} >= "
            f"{config.bed_divergence_spread_trigger} over {n} instances "
            f"(wide divergence band ⇒ likely divergent tail)")
    flags['bed_divergence_spread'] = spread_fired

    # Signal 2: copy count vs length class.
    lc_fired = copy_count_vs_length_class(rec_view, config)
    if lc_fired:
        length = rec_view.get('actual_length') or rec_view.get('length') or 0
        expect = _expected_copies(length,
                                  getattr(config, 'length_class_copy_expectation', ()))
        reasons.append(
            f"copy_count_vs_length_class: n_instances={n} < expected {expect} "
            f"for length={length}bp (under-instanced for length class)")
    flags['copy_count_vs_length_class'] = lc_fired

    # Signal 3: under-instanced (M5 trigger, demoted).
    ui_fired = under_instanced(rec_view, config)
    if ui_fired:
        reasons.append(
            f"under_instanced: n_instances={n} < min_copies_for_msa="
            f"{config.min_copies_for_msa}")
    flags['under_instanced'] = ui_fired

    # N4: sketch pre-screen signal EMPIRICALLY DISABLED (see module docstring +
    # config.enable_genome_sketch_prescreen) — sourmash is blind to the 75-90%
    # divergent tail, so it cannot gate recall.  The selective gate is the
    # specific-signal co-fire in recall_eligible(), not a sketch.
    # N4: boundary-truncation signal still deferred (needs first-pass MSA occupancy);
    # not required for the recall gate.

    eligible = any(flags.values())
    return CompletenessVerdict(
        record_id=record.get('id', ''),
        incomplete_recall_eligible=eligible,
        flags=flags,
        reasons=reasons,
        n_instances=n,
        bed_divergence_spread=spread,
    )


# ═══════════════════════════════════════════════════════════════════════
# N4 selective recall gate (co-fire — keeps the bottleneck out)
# ═══════════════════════════════════════════════════════════════════════

def recall_eligible(record: Dict, instances, config) -> CompletenessVerdict:
    """SELECTIVE N4 recall gate: should this family pay for a Tier-2 blastn recall?

    Unlike `assess_completeness` (which sets `incomplete_recall_eligible` when ANY
    signal fires — useful as an audit of *all* families with a hint of
    incompleteness), this gate requires the EXPENSIVE-recall trigger to be the
    co-fire of `under_instanced` with a SPECIFIC signal:

        recall_eligible  ⇔  under_instanced  AND  (bed_divergence_spread
                                                    OR copy_count_vs_length_class)

    Why co-fire and not `under_instanced` alone (the implemented M5 trigger):
    `under_instanced` (BED instances < min_copies_for_msa) fires for 85.2% of real
    chr4 families because that genome is 84% low-copy — it is the NORM, not a
    missing-tail signal, so triggering O(genome) blastn on it requests exactly the
    O(families × genome) sweep the v1 design correctly fought
    (REFINE_STRATEGY_DESIGN_v2.md §2.5).  Requiring a SPECIFIC co-fire signal —
    bed_divergence_spread (the family already straddles a wide divergence band ⇒ a
    tail below mdl-repeat's acceptance likely exists; the R=150 case, spread 0.328)
    or copy_count_vs_length_class (a long element with implausibly few copies) —
    drops the recall set to 20.9% on real chr4 while still firing on R=150.

    The sketch pre-screen (§2.5 Tier 1) that would otherwise be a third, cheaper
    confirm is empirically disabled (module docstring); this gate IS the selective
    layer.  It is still O(len(instances)) + O(1) per family — no genome alignment to
    decide *whether* to recall (the alignment is Tier 2, only for the families this
    gate passes).

    Returns a CompletenessVerdict whose `incomplete_recall_eligible` carries the
    CO-FIRE decision (not the any-signal one).  `flags`/`reasons` still record every
    individual signal so the metric output can audit why a family did or did not
    recall.  When `config.enable_selective_recall` is False, NOTHING is eligible —
    the pipeline degrades to N3 (cluster BED copies only), the safe-disable invariant.
    """
    verdict = assess_completeness(record, instances, config)

    if not getattr(config, 'enable_selective_recall', True):
        verdict.incomplete_recall_eligible = False
        verdict.reasons.append("selective recall disabled (enable_selective_recall=False)")
        return verdict

    ui = verdict.flags.get('under_instanced', False)
    specific = (verdict.flags.get('bed_divergence_spread', False)
                or verdict.flags.get('copy_count_vs_length_class', False))
    cofire = bool(ui and specific)

    verdict.incomplete_recall_eligible = cofire
    if cofire:
        verdict.reasons.append(
            "RECALL-ELIGIBLE: under_instanced co-fires with a specific signal "
            "(divergent tail likely; spend Tier-2 blastn here)")
    elif ui and not specific:
        verdict.reasons.append(
            "not recall-eligible: under_instanced alone (the 84%-low-copy norm) — "
            "no specific divergent-tail signal, recall would be the O(genome) bottleneck")
    return verdict


# ═══════════════════════════════════════════════════════════════════════
# Pre-signal TSV reader (consume what Phase 0 emitted)
# ═══════════════════════════════════════════════════════════════════════

def load_presignals(presignal_path: str) -> Dict[str, Dict]:
    """Parse completeness_presignals.tsv → {record_id: presignal-row dict}.

    Row dict keys: n_instances (int), bed_divergence_mean (float|None),
    bed_divergence_spread (float|None), length (int),
    length_class_expectation (int), length_class_flag (bool).
    """
    out: Dict[str, Dict] = {}
    with open(presignal_path) as fh:
        header = fh.readline()  # skip header
        if not header:
            return out
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 7:
                continue
            rid = parts[0]
            out[rid] = {
                'n_instances': int(parts[1]),
                'bed_divergence_mean': float(parts[2]) if parts[2] else None,
                'bed_divergence_spread': float(parts[3]) if parts[3] else None,
                'length': int(parts[4]),
                'length_class_expectation': int(parts[5]),
                'length_class_flag': parts[6] == '1',
            }
    return out
