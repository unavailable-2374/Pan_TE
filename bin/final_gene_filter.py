#!/usr/bin/env python3
"""Final library step: remove host-gene sequence from the classified TE library.

Multi-copy host gene families (F-box, LRR-RLK, PPR, lectins, ...) are recurrent and get
collected as repeats; they mask coding sequence without being TEs. This step runs on the
finished library, against a cellular protein DB with TE proteins removed:

1. Drop by the single-HSP rule (Refiner_mdl.phase2_gene_exclusion): coverage >= 0.5 with no
   TE evidence, or coverage >= 0.7 regardless of evidence.
2. Drop by merged coverage: HSPs with identity >= 35% are merged per protein, because a
   multi-exon gene yields several short HSPs that the single-HSP rule misses. Drop at
   coverage >= 0.1 with no TE evidence, or >= 0.5 with TE evidence.
3. Excise from every remaining entry the regions hit by any cellular protein, extended by
   30 bp on each side. Segments of at least 100 bp stay as separate entries named
   <id>_s<k>#<class>; the parent is recoverable from the name.

The step has no skip path: a missing DB or a failed BLAST aborts the run, so a library is
never reported as filtered when it was not.

Usage: final_gene_filter.py <library.fa> <cellular_db> <output_dir> <threads>
"""
import collections
import concurrent.futures
import logging
import os
import shutil
import subprocess
import sys
import tempfile
from typing import Dict, List, Tuple

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from Refiner_mdl.config import RefinerMdlConfig  # noqa: E402
from Refiner_mdl.phase2_gene_exclusion import _cellular_db_ready, exclude_gene_derived  # noqa: E402
from Refiner_mdl.te_structure_filter import annotate_te_structure  # noqa: E402

logger = logging.getLogger(__name__)

NO_EVIDENCE_MIN_COV = 0.1
EVIDENCE_MIN_COV = 0.5
MIN_PIDENT = 35.0
MARGIN = 30
MIN_SEGMENT = 100
HSP_EVALUE = "1e-5"
HSP_MAX_TARGETS = "20"

Interval = Tuple[int, int]
# {query: {protein: [query intervals]}}
Hsps = Dict[str, Dict[str, List[Interval]]]


def resolve_cellular_db() -> str:
    """Locate the cellular protein DB; returns '' when no usable DB exists."""
    env = os.environ.get("PAN_TE_CELLULAR_DB")
    if env:
        return env if _cellular_db_ready(env) else ""
    candidates = []
    prefix = os.environ.get("CONDA_PREFIX")
    if prefix:
        candidates.append(os.path.join(prefix, "share", "pan_te_data", "sprot_cell"))
    candidates.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "share",
                                   "pan_te_data", "sprot_cell"))
    return next((c for c in candidates if _cellular_db_ready(c)), "")


def merge_intervals(intervals: List[Interval]) -> List[List[int]]:
    merged: List[List[int]] = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return merged


def merged_coverage(by_protein: Dict[str, List[Interval]], qlen: int) -> float:
    """Largest fraction of the query covered by the merged HSPs of a single protein."""
    if not by_protein or qlen <= 0:
        return 0.0
    best = max(sum(e - s + 1 for s, e in merge_intervals(v)) for v in by_protein.values())
    return min(1.0, best / qlen)


def excise_segments(seq_len: int, by_protein: Dict[str, List[Interval]]) -> List[Interval]:
    """1-based segments left after removing protein-hit regions plus the margin."""
    hits = [(s - MARGIN, e + MARGIN) for v in by_protein.values() for s, e in v]
    segments, pos = [], 1
    for start, end in merge_intervals(hits):
        if start > pos:
            segments.append((pos, start - 1))
        pos = max(pos, end + 1)
    if pos <= seq_len:
        segments.append((pos, seq_len))
    clipped = [(s, min(e, seq_len)) for s, e in segments]
    return [(s, e) for s, e in clipped if e - s + 1 >= MIN_SEGMENT]


def segment_id(entry_id: str, k: int) -> str:
    name, sep, cls = entry_id.partition("#")
    return f"{name}_s{k}{sep}{cls}"


def read_fasta(path: str) -> List[Dict]:
    records: List[Dict] = []
    for line in open(path):
        if line.startswith(">"):
            records.append({"id": line[1:].split()[0], "sequence": []})
        elif records:
            records[-1]["sequence"].append(line.strip())
    for record in records:
        record["sequence"] = "".join(record["sequence"])
    if len({r["id"] for r in records}) != len(records):
        raise ValueError(f"duplicate sequence ids in {path}")
    return records


def run_hsp_blastx(records: List[Dict], db: str, threads: int, work_dir: str,
                   blastx_exe: str = "blastx") -> Hsps:
    """All HSPs at the identity floor. Any BLAST failure raises."""
    n_jobs = max(1, min(threads, len(records)))
    chunks = [os.path.join(work_dir, f"hsp_q{j}.fa") for j in range(n_jobs)]
    handles = [open(p, "w") for p in chunks]
    for i, record in enumerate(records):
        handles[i % n_jobs].write(f">{record['id']}\n{record['sequence']}\n")
    for handle in handles:
        handle.close()

    def blast(query: str) -> str:
        out = query + ".tsv"
        result = subprocess.run(
            [blastx_exe, "-query", query, "-db", db, "-evalue", HSP_EVALUE,
             "-max_target_seqs", HSP_MAX_TARGETS, "-num_threads", "1",
             "-outfmt", "6 qseqid qstart qend sseqid pident", "-out", out],
            capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"blastx failed ({result.returncode}): {result.stderr[:300]}")
        return out

    hsps: Hsps = collections.defaultdict(lambda: collections.defaultdict(list))
    with concurrent.futures.ThreadPoolExecutor(max_workers=n_jobs) as pool:
        for out in pool.map(blast, chunks):
            for line in open(out):
                query, start, end, protein, pident = line.rstrip("\n").split("\t")
                if float(pident) >= MIN_PIDENT:
                    a, b = sorted((int(start), int(end)))
                    hsps[query][protein].append((a, b))
    return hsps


def apply_rules(records: List[Dict], signals: Dict[str, List[str]], single_hsp_dropped: set,
                hsps: Hsps) -> Tuple[List[Tuple[str, str]], List[Tuple[str, float, str]],
                                     List[Tuple[str, int, int, int]]]:
    """Returns (library entries, dropped rows, excised rows)."""
    library, dropped, excised = [], [], []
    for record in records:
        rid, seq = record["id"], record["sequence"]
        by_protein = hsps.get(rid, {})
        cov = merged_coverage(by_protein, len(seq))
        threshold = EVIDENCE_MIN_COV if signals.get(rid) else NO_EVIDENCE_MIN_COV
        if rid in single_hsp_dropped:
            dropped.append((rid, cov, "single_hsp"))
        elif cov >= threshold:
            dropped.append((rid, cov, "merged_coverage"))
        elif not by_protein:
            library.append((rid, seq))
        else:
            segments = excise_segments(len(seq), by_protein)
            for k, (start, end) in enumerate(segments, 1):
                library.append((segment_id(rid, k), seq[start - 1:end]))
            excised.append((rid, len(seq), len(segments), sum(e - s + 1 for s, e in segments)))
    return library, dropped, excised


def filter_library(library_path: str, db: str, out_dir: str, threads: int) -> Dict:
    if not _cellular_db_ready(db):
        raise FileNotFoundError(f"cellular protein DB not found: {db or '(unset)'}; "
                                f"run: pan_te-setup-data --cellular-db <DIR>")
    records = read_fasta(library_path)
    if not records:
        raise ValueError(f"empty TE library: {library_path}")
    os.makedirs(out_dir, exist_ok=True)
    config = RefinerMdlConfig(threads=threads, cellular_protein_db=db, enable_gene_exclusion=True)
    signals, _ = annotate_te_structure(records, config)
    _, single_dropped, gene_stats = exclude_gene_derived(records, signals, config)
    if gene_stats.get("skipped_reason"):
        raise RuntimeError(f"gene exclusion did not run: {gene_stats['skipped_reason']}")
    work_dir = tempfile.mkdtemp(prefix="final_gene_filter_", dir=out_dir)
    try:
        hsps = run_hsp_blastx(records, db, threads, work_dir)
    finally:
        shutil.rmtree(work_dir, ignore_errors=True)
    library, dropped, excised = apply_rules(records, signals, {r["id"] for r in single_dropped}, hsps)
    if not library:
        raise RuntimeError("gene filter removed every library entry")

    tmp = os.path.join(out_dir, "TEs.fa.tmp")
    with open(tmp, "w") as out:
        for rid, seq in library:
            out.write(f">{rid}\n{seq}\n")
    os.replace(tmp, os.path.join(out_dir, "TEs.fa"))
    with open(os.path.join(out_dir, "gene_filter_dropped.tsv"), "w") as out:
        out.write("id\tmerged_coverage\trule\tte_signals\n")
        for rid, cov, rule in dropped:
            out.write(f"{rid}\t{cov:.3f}\t{rule}\t{','.join(signals.get(rid, []))}\n")
    with open(os.path.join(out_dir, "gene_filter_excised.tsv"), "w") as out:
        out.write("id\tlength\tsegments_kept\tbp_kept\n")
        for row in excised:
            out.write("\t".join(map(str, row)) + "\n")
    stats = {"input": len(records), "input_bp": sum(len(r["sequence"]) for r in records),
             "dropped": len(dropped), "excised": len(excised), "output": len(library),
             "output_bp": sum(len(s) for _, s in library)}
    with open(os.path.join(out_dir, "gene_filter_stats.tsv"), "w") as out:
        for key, value in stats.items():
            out.write(f"{key}\t{value}\n")
    logger.info("Gene filter: %s", stats)
    return stats


if __name__ == "__main__":
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    print(filter_library(sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])))
