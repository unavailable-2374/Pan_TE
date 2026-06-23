"""
Phase 2: Short low-copy noise filter (genome-size adaptive).

mdl-repeat's MDL criterion reports ANY sequence whose family encoding beats literal
encoding, so its raw output is dominated (~80%) by SHORT, LOW-COPY fragments (median
~120 bp, 2-4 copies) — the bottom of the MDL barrel: marginal recent duplications /
segmental fragments, not confident TE families. They are ~50% of the refined library by
count but contribute <10% of genome masking, and are too short + too few-copy to be
caught by TE-looker (which targets divergent families at >=5 copies) either.

This filter drops only the SHORT-AND-LOW-COPY corner:

    drop  ⇔  copies < min_copies(genome_size)  AND  length < lowcopy_max_len

The JOINT condition protects the two legitimate edges seen in the data:
  * long low-copy  (e.g. a 2-copy 8 kb full-length element)      — length >= max_len keeps it
  * short high-copy (e.g. a 120 bp MITE at 50 copies)            — copies >= min_copies keeps it

min_copies scales with genome size: the number of chance / segmental occurrences of a
short sequence grows with genome length, so the copy floor to call a real family rises.
(Tiers mirror the genome-size tiering already used for window-stride.)
"""

import logging
from typing import Dict, List, Tuple

logger = logging.getLogger(__name__)


def lowcopy_min_copies(config) -> int:
    """Genome-size-adaptive minimum copy count to keep a short sequence."""
    if getattr(config, 'lowcopy_min_copies_override', 0) > 0:
        return config.lowcopy_min_copies_override
    g = getattr(config, 'genome_size_bp', 0) or 0
    if g <= 0:
        return 3                       # unknown genome -> conservative base
    mb = g / 1e6
    if mb < 100:
        return 3                       # XS
    if mb < 500:
        return 4                       # S-M  (e.g. Arabidopsis 134 Mb)
    if mb < 1000:
        return 5                       # L
    if mb < 3000:
        return 6                       # XL
    return 8                           # XXL


def filter_short_lowcopy(records: List[Dict], config) -> Tuple[List[Dict], List[Dict], Dict]:
    """Drop short-AND-low-copy consensi (marginal MDL noise).

    Returns (kept, dropped, stats). Dropped records get a 'lowcopy_noise' annotation.
    """
    stats = {'input': len(records),
             'enabled': bool(getattr(config, 'enable_lowcopy_filter', False))}
    if not records or not getattr(config, 'enable_lowcopy_filter', False):
        stats['skipped_reason'] = 'disabled' if records else 'no_records'
        return records, [], stats

    # Large-genome SAMPLED mode: copy counts here come from the SAMPLE (undercount), so
    # the real floor is deferred to phase2b, which recruits genome-wide counts on the
    # final consensi. Keep every record (the copy gate then only recruits + rescues).
    if getattr(config, 'defer_copy_floor', False):
        stats['skipped_reason'] = 'deferred_to_phase2b (sampled mode)'
        logger.info("Low-copy filter DEFERRED to phase2b (sampled mode); keeping all "
                    f"{len(records)} families for now")
        return records, [], stats

    min_c = lowcopy_min_copies(config)
    max_l = config.lowcopy_max_len
    hard_c = getattr(config, 'hard_min_copies', 0)
    kept, dropped = [], []
    n_hard = n_joint = 0
    use_genomic = any('genomic_copies' in rec for rec in records)
    for rec in records:
        # Prefer the rmblastn-recruited genomic copy count over mdl's unreliable estimate.
        c = rec['genomic_copies'] if 'genomic_copies' in rec else rec.get('copies', 0)
        l = len(rec.get('sequence', ''))
        hard_hit = hard_c > 0 and c < hard_c                 # hard recurrence floor
        joint_hit = c < min_c and l < max_l                  # short AND low-copy gate
        if hard_hit or joint_hit:
            rec['lowcopy_noise'] = {'copies': c, 'length': l,
                                    'reason': 'hard_floor' if hard_hit else 'short_lowcopy'}
            dropped.append(rec)
            n_hard += hard_hit
            n_joint += (joint_hit and not hard_hit)
        else:
            kept.append(rec)

    stats.update({'min_copies': min_c, 'max_len': max_l, 'hard_min_copies': hard_c,
                  'genome_size_bp': getattr(config, 'genome_size_bp', 0),
                  'copy_source': 'genomic_copies(rmblastn)' if use_genomic else 'mdl_copies',
                  'dropped': len(dropped), 'dropped_hard_floor': n_hard,
                  'dropped_short_lowcopy': n_joint, 'kept': len(kept)})
    logger.info(f"Low-copy filter (genome {stats['genome_size_bp']/1e6:.0f} Mb -> "
                f"min_copies={min_c}, max_len={max_l}, hard_floor={hard_c or 'off'}): "
                f"{len(records)} -> kept {len(kept)}, dropped {len(dropped)} "
                f"(hard_floor={n_hard}, short_lowcopy={n_joint})")
    return kept, dropped, stats
