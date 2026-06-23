"""
Copy gate (between Phase 0 and Phase 1): recruit real genomic copies, apply the
recurrence hard-floor, then rescue dropped-but-TE-structured consensi.

Runs BEFORE the expensive representativeness work (Phase 1 subfamily clustering /
consensus / recall) so that:
  * the floor thresholds on the WHOLE family's genomic copy count (rmblastn), not the
    per-subfamily count that splitting would dilute below the floor;
  * Phase 1 only refines the survivors, not families about to be dropped.

Flow:  recruit_genomic_copies -> filter_short_lowcopy(hard floor) -> rescue dropped
        -> survivors = kept ∪ rescued  (records dict otherwise preserved for Phase 1).
"""

import logging
from typing import Dict

logger = logging.getLogger(__name__)


def run_copy_gate(phase0_output: Dict, config) -> Dict:
    from phase2_copy_recruit import recruit_genomic_copies
    from phase2_lowcopy_filter import filter_short_lowcopy
    from phase_rescue import rescue_te_structured

    records = phase0_output.get('records', [])
    if not records or not getattr(config, 'enable_copy_recruit', False):
        return phase0_output

    logger.info("=" * 60)
    logger.info(f"Copy gate: recruit + hard-floor + rescue ({len(records)} families)")
    logger.info("=" * 60)

    # 1) real genomic copy number (writes rec['genomic_copies'])
    recruit_stats = recruit_genomic_copies(records, config)

    # 2) recurrence hard-floor / short-low-copy filter (reads genomic_copies)
    kept, dropped, lc_stats = filter_short_lowcopy(records, config)

    # 3) rescue dropped consensi that carry TE evidence (HMM + structure)
    rescued, rescue_stats = rescue_te_structured(dropped, config)

    survivors = kept + rescued
    gate_stats = {'input': len(records), 'kept_by_floor': len(kept),
                  'dropped': len(dropped), 'rescued': len(rescued),
                  'survivors': len(survivors),
                  'recruit': recruit_stats, 'lowcopy': lc_stats, 'rescue': rescue_stats}
    logger.info(f"Copy gate: {len(records)} -> {len(survivors)} "
                f"(floor-kept {len(kept)} + rescued {len(rescued)} of {len(dropped)} dropped)")

    out = dict(phase0_output)
    out['records'] = survivors
    out['copy_gate_stats'] = gate_stats
    return out
