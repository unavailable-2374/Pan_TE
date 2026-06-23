"""
Phase 1 (rewrite) — Module M5: bounded per-family blastn recruitment + lazy
whole-genome BLAST DB construction.

For under-instanced families (fewer than `min_copies_for_msa` BED instances) the
BED-seeded extractor has too few copies to build a trustworthy MSA. This module
runs ONE bounded blastn of the family seed against a lazily-built whole-genome DB
and extracts the (padded, strand-oriented) subject spans as extra copies. The cost
is O(few under-instanced families × genome) — never O(all families × genome) —
because the common path (families with enough BED instances) never touches the DB.

Two strict invariants (plan §1.5, §7, §11):
  * The DB is built ONCE in the PARENT process (run_phase1, before workers spawn) and
    its path is propagated into the worker config_dict. Workers NEVER build it — that
    would be a race on shared files. `ensure_genome_blast_db` is idempotent.
  * Any failure (no DB, blastn non-zero exit, no hits, extraction failure) returns an
    empty copy list → the caller (`refine_family`) sees < min_copies_for_msa copies
    and keeps the original mdl-repeat seed. No copy is ever fabricated, no faked
    success.

N4 (v2 §2.5, §7 M5-promote row) promotes `recruit_by_blastn` to **Tier-2 of
selective recall**: the SAME bounded recruiter, but invoked whenever the cheap
co-fire completeness gate flags a family (phase1_completeness.recall_eligible),
not only when BED < min_copies.  A `recall_mode=True` call swaps the pident floor
to `config.recall_identity_floor` (the ~60-70% refine↔TE-looker handoff band,
§2.6) and the subprocess wall to `config.recall_per_family_wall_s`; everything else
(one bounded blastn, cap at blastn_max_targets by score, same M1 padded extractor →
identical `Copy` objects, any-failure→[] never-fabricate) is unchanged.  The M5
default-arg path (recall_mode=False) is byte-for-byte the original behaviour so the
M5 suite stays green.

See REFINE_IMPLEMENTATION_PLAN.md §1.5, §7, §10, §12.M5 and
REFINE_STRATEGY_DESIGN_v2.md §2.5, §2.6.
"""

import logging
import os
import subprocess
import tempfile
from typing import List, Tuple

from phase1_extract import (
    Copy,
    Instance,
    compute_pad,
    extract_padded_copies,
    load_fai_lengths,
    strip_seqid_prefix,
)

logger = logging.getLogger(__name__)

# Divergent-recruitment pident floor: with enable_divergent_blast_recruitment the
# 45-75% identity band is opened (otherwise that band is TE-looker's remit, §7).
_DIVERGENT_PIDENT_FLOOR = 45.0


def ensure_genome_blast_db(config) -> str:
    """Build the whole-genome BLAST DB once (lazily) and return its path.

    Idempotent: if the DB already exists (or config.genome_blast_db is already set to
    an existing DB) it is returned without rebuilding. Built in the parent process
    only. On any failure config.genome_blast_db is set to "" and an ERROR is logged —
    under-instanced families then keep their seeds; the common path is unaffected
    because it never reads the DB.

    Deliberately NO short timeout (the old 600 s trap turned a slow large-genome
    makeblastdb into a silent failure). A generous ceiling guards against a truly
    hung process without misclassifying a long-but-progressing build as a failure."""
    existing = getattr(config, 'genome_blast_db', '') or ''
    if existing and (os.path.exists(existing + '.nsq')
                     or os.path.exists(existing + '.nhr')):
        logger.info("Genome BLAST DB already available: %s", existing)
        return existing

    db_dir = os.path.join(config.temp_dir, 'genome_blastdb')
    db_path = os.path.join(db_dir, 'genome')
    if os.path.exists(db_path + '.nsq') or os.path.exists(db_path + '.nhr'):
        logger.info("Genome BLAST DB found on disk: %s", db_path)
        config.genome_blast_db = db_path
        return db_path

    os.makedirs(db_dir, exist_ok=True)
    cmd = [
        config.makeblastdb_exe,
        '-in', config.genome_file,
        '-dbtype', 'nucl',
        '-out', db_path,
        '-parse_seqids',
    ]
    logger.info("Building genome BLAST DB (lazy, fallback-triggered): %s",
                ' '.join(cmd))
    try:
        # Generous ceiling, not the old 600 s trap. A real failure logs ERROR below.
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=7200)
    except Exception as e:  # noqa: BLE001 — never fake success; disable fallback
        logger.error("Genome BLAST DB build failed (%s); fallback recruitment "
                     "disabled, under-instanced families keep their seeds", e)
        config.genome_blast_db = ""
        return ""
    if result.returncode != 0:
        logger.error("makeblastdb returned %d (%s); fallback recruitment disabled, "
                     "under-instanced families keep their seeds",
                     result.returncode, (result.stderr or '')[:300])
        config.genome_blast_db = ""
        return ""
    config.genome_blast_db = db_path
    logger.info("Genome BLAST DB built successfully: %s", db_path)
    return db_path


def _blast_hits_to_instances(stdout: str, pident_min: float,
                             max_keep: int = None) -> List[Instance]:
    """Parse `-outfmt '6 sseqid sstart send pident length [score]'` into Instances.

    Each HSP becomes one Instance: chrom is the seqid with any NCBI `ref|...|`
    wrapper stripped, BED 0-based half-open coords are derived from the 1-based
    inclusive blast span, and strand comes from the sign of (sstart vs send). The
    per-copy divergence proxy is 1 - pident/100. Hits below pident_min are dropped.

    When `max_keep` is set, the hits are ranked by alignment score (column 6 if
    present, else 0) descending and only the top `max_keep` are returned — this is
    how the recruiter bounds total recruited copies to `blastn_max_targets`
    independent of genome contiguity (a single chromosome with many TE copies yields
    many HSPs; -max_target_seqs alone would not cap them)."""
    scored: List[Tuple[float, Instance]] = []
    for line in stdout.strip().split('\n'):
        if not line:
            continue
        f = line.split('\t')
        if len(f) < 4:
            continue
        try:
            sseqid = strip_seqid_prefix(f[0])
            sstart = int(f[1])
            send = int(f[2])
            pident = float(f[3])
            score = float(f[5]) if len(f) >= 6 else 0.0
        except (ValueError, IndexError):
            continue
        if pident < pident_min:
            continue
        strand = '+' if sstart <= send else '-'
        lo, hi = (sstart, send) if sstart <= send else (send, sstart)
        # 1-based inclusive [lo, hi] -> 0-based half-open [lo-1, hi).
        scored.append((score, Instance(chrom=sseqid, start=lo - 1, end=hi,
                                       strand=strand,
                                       divergence=1.0 - pident / 100.0)))
    if max_keep is not None and len(scored) > max_keep:
        scored.sort(key=lambda x: x[0], reverse=True)
        scored = scored[:max_keep]
    return [inst for _score, inst in scored]


def _recall_pident_floor(config) -> float:
    """Tier-2 recall pident floor (%): config.recall_identity_floor (fraction) → %.

    The N4 recall band (§2.6): recruit seed-homologous members down to the
    productive ~60-70% identity floor.  Stored as a fraction in config
    (recall_identity_floor=0.65) and converted to the %-scale blastn pident here.
    Falls back to the M5 divergent/standard floors if the field is absent."""
    floor = getattr(config, 'recall_identity_floor', None)
    if floor is None:
        return (_DIVERGENT_PIDENT_FLOOR
                if getattr(config, 'enable_divergent_blast_recruitment', False)
                else config.min_recruit_pident)
    return float(floor) * 100.0


def recruit_by_blastn(rec, config, recall_mode: bool = False) -> List[Copy]:
    """Bounded per-family copy recruitment (§1.5, §7; N4 Tier-2 §2.5, §2.6).

    Runs ONE blastn of rec['sequence'] against the lazily-built whole-genome DB,
    parses the HSP spans into Instances, and reuses the M1 padded extractor to return
    strand-co-oriented Copy objects. Any failure → [] (caller keeps the seed).

    recall_mode (N4): when True this is Tier-2 of selective recall — the pident floor
    becomes config.recall_identity_floor (the ~60-70% refine↔TE-looker handoff band,
    §2.6) and the subprocess wall becomes config.recall_per_family_wall_s.  When False
    (the M5 default) the floor is min_recruit_pident / the divergent floor and the wall
    is the original 600 s — byte-for-byte the pre-N4 behaviour."""
    db = getattr(config, 'genome_blast_db', '') or ''
    if not db:
        logger.warning("recruit_by_blastn: no genome DB for family %s; keeping seed",
                       rec.get('id', '?'))
        return []

    seq = rec.get('sequence', '')
    if not seq:
        return []

    if recall_mode:
        pident_min = _recall_pident_floor(config)
        wall_s = int(getattr(config, 'recall_per_family_wall_s', 300))
    else:
        pident_min = (_DIVERGENT_PIDENT_FLOOR
                      if getattr(config, 'enable_divergent_blast_recruitment', False)
                      else config.min_recruit_pident)
        wall_s = 600

    tmp_dir = tempfile.mkdtemp(prefix='recruit_blastn_')
    query_path = os.path.join(tmp_dir, 'query.fa')
    try:
        with open(query_path, 'w') as fh:
            fh.write(f">{rec.get('id', 'query')}\n{seq}\n")

        # NOTE (deviation from plan §1.5/§7, documented): the plan listed
        # `-max_hsps 1`, but that caps BLAST to ONE HSP per *subject sequence*
        # (chromosome). A TE family typically has many copies on the same
        # chromosome, so `-max_hsps 1` would recruit at most one copy per chromosome
        # and defeat the recruiter's purpose on contiguous genomes (verified: 1 of 10
        # planted copies recovered on a single-chromosome synthetic genome). The
        # plan's clear intent (`-max_target_seqs 50` → up to ~50 recruited copies) is
        # honored instead by allowing up to blastn_max_targets HSPs per subject and
        # capping the TOTAL recruited copies to blastn_max_targets by alignment score
        # (`_blast_hits_to_instances(max_keep=...)`). This is a root-cause fix, not a
        # masking workaround.
        max_targets = config.blastn_max_targets
        cmd = [
            config.blastn_exe,
            '-query', query_path,
            '-db', db,
            '-outfmt', '6 sseqid sstart send pident length score',
            '-max_target_seqs', str(max_targets),
            '-max_hsps', str(max_targets),
            '-evalue', str(config.blastn_evalue),
            '-dust', 'no',
            '-num_threads', '1',
        ]
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=wall_s)
        except Exception as e:  # noqa: BLE001 — report + fall back, never fake a copy
            # Includes subprocess.TimeoutExpired when a recall exceeds wall_s: the
            # family keeps its BED copies (caller flags completeness_verified=False);
            # never a fabricated copy, never a faked success.
            logger.warning("recruit_by_blastn: blastn failed/timed-out (wall=%ds) for "
                           "family %s (%s)", wall_s, rec.get('id', '?'), e)
            return []
        if result.returncode != 0:
            logger.warning("recruit_by_blastn: blastn exit %d for family %s: %s",
                           result.returncode, rec.get('id', '?'),
                           (result.stderr or '')[:200])
            return []

        instances = _blast_hits_to_instances(result.stdout, pident_min,
                                             max_keep=config.blastn_max_targets)
        if not instances:
            logger.info("recruit_by_blastn: no qualifying hits (pident>=%.0f) for "
                        "family %s; keeping seed", pident_min, rec.get('id', '?'))
            return []

        fai = load_fai_lengths(config.genome_file + '.fai')
        copies = extract_padded_copies(
            instances, config.genome_file, fai,
            lambda L: compute_pad(L, config), config)
        logger.info("recruit_by_blastn: family %s recruited %d copies from %d hits "
                    "(pident>=%.0f)", rec.get('id', '?'), len(copies),
                    len(instances), pident_min)
        return copies
    finally:
        import shutil
        shutil.rmtree(tmp_dir, ignore_errors=True)
