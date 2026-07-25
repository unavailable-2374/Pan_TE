#!/usr/bin/env python3
"""Compare GFA-node RepeatMasker annotations with a direct reference run.

All internal coordinates are zero-based, half-open.  Only repeats present in
the supplied custom library are retained from the direct RepeatMasker output;
this excludes RepeatMasker's default simple/low-complexity annotations.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import json
import statistics
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple


@dataclass(frozen=True)
class Hit:
    chrom: str
    start: int
    end: int
    family: str
    te_class: str
    score: float
    source_id: str
    node_id: str = ""


def reverse_complement(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def read_fasta(path: Path) -> Dict[str, str]:
    records: Dict[str, List[str]] = {}
    name = None
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].split()[0]
                if name in records:
                    raise ValueError(f"Duplicate FASTA record: {name}")
                records[name] = []
            elif name is None:
                raise ValueError(f"Sequence before FASTA header in {path}")
            else:
                records[name].append(line)
    return {name: "".join(parts).upper() for name, parts in records.items()}


def read_library(path: Path) -> Dict[str, str]:
    families: Dict[str, str] = {}
    with path.open() as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            token = line[1:].strip().split()[0]
            family, sep, te_class = token.partition("#")
            families[family] = te_class if sep else "Unknown"
    if not families:
        raise ValueError(f"No FASTA headers found in repeat library: {path}")
    return families


def parse_tags(fields: Sequence[str]) -> Dict[str, str]:
    tags = {}
    for field in fields:
        parts = field.split(":", 2)
        if len(parts) == 3:
            tags[parts[0]] = parts[2]
    return tags


def reference_nodes(gfa: Path, reference: Dict[str, str]):
    nodes = {}
    orientation_counts = defaultdict(int)
    intervals = defaultdict(list)
    with gfa.open() as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.startswith("S\t"):
                continue
            fields = raw.rstrip("\n").split("\t")
            tags = parse_tags(fields[3:])
            if tags.get("SR") != "0":
                continue
            node_id, sequence = fields[1], fields[2].upper()
            if sequence == "*":
                raise ValueError(f"Reference node {node_id} has no sequence")
            try:
                chrom, offset = tags["SN"].split("#")[-1], int(tags["SO"])
            except (KeyError, ValueError) as exc:
                raise ValueError(f"Invalid SN/SO at GFA line {line_number}") from exc
            if chrom not in reference:
                raise ValueError(f"GFA chromosome {chrom!r} is absent from reference FASTA")
            target = reference[chrom][offset:offset + len(sequence)]
            if sequence == target:
                orientation = "+"
            elif reverse_complement(sequence) == target:
                orientation = "-"
            else:
                orientation = "mismatch"
            orientation_counts[orientation] += 1
            nodes[node_id] = (chrom, offset, len(sequence), orientation)
            intervals[chrom].append((offset, offset + len(sequence)))

    if orientation_counts["mismatch"]:
        raise ValueError(
            f"{orientation_counts['mismatch']} SR=0 nodes do not match the reference "
            "at SN/SO in either orientation"
        )
    if not nodes:
        raise ValueError("No SR:i:0 reference nodes found in GFA")
    return nodes, dict(orientation_counts), intervals


def parse_projected_hits(path: Path, nodes) -> List[Hit]:
    hits = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"node_id", "node_start", "node_end", "te_family", "te_class", "score"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(f"Missing required columns in {path}")
        for row_number, row in enumerate(reader, 2):
            node = nodes.get(row["node_id"])
            if node is None:
                continue
            chrom, offset, length, orientation = node
            start, end = int(row["node_start"]), int(row["node_end"])
            if not (0 <= start < end <= length):
                raise ValueError(f"Invalid node interval at {path}:{row_number}")
            if orientation == "+":
                ref_start, ref_end = offset + start, offset + end
            else:
                ref_start, ref_end = offset + length - end, offset + length - start
            hits.append(Hit(
                chrom, ref_start, ref_end, row["te_family"], row["te_class"],
                float(row["score"]), f"gfa:{row_number}", row["node_id"]
            ))
    return sorted(hits, key=lambda h: (h.chrom, h.start, h.end, h.family))


def normalize_query(query: str, reference: Dict[str, str]) -> str:
    if query in reference:
        return query
    suffix = query.split("#")[-1]
    if suffix in reference:
        return suffix
    raise ValueError(f"RepeatMasker query {query!r} cannot be mapped to reference FASTA")


def parse_repeatmasker_out(path: Path, library, reference) -> Tuple[List[Hit], dict]:
    hits = []
    counts = defaultdict(int)
    with path.open() as handle:
        for line_number, raw in enumerate(handle, 1):
            fields = raw.split()
            if len(fields) < 14 or not fields[0].isdigit():
                continue
            counts["all_records"] += 1
            family = fields[9]
            if family not in library:
                counts["not_in_custom_library"] += 1
                continue
            start, end = int(fields[5]) - 1, int(fields[6])
            if start >= end:
                counts["invalid_intervals"] += 1
                continue
            chrom = normalize_query(fields[4], reference)
            if end > len(reference[chrom]):
                raise ValueError(f"Interval exceeds {chrom} at {path}:{line_number}")
            hits.append(Hit(
                chrom, start, end, family, fields[10], float(fields[0]),
                f"direct:{line_number}"
            ))
            counts["custom_library_records"] += 1
    return sorted(hits, key=lambda h: (h.chrom, h.start, h.end, h.family)), dict(counts)


def merge_intervals(intervals: Iterable[Tuple[int, int]]) -> List[Tuple[int, int]]:
    merged: List[List[int]] = []
    for start, end in sorted(intervals):
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        elif end > merged[-1][1]:
            merged[-1][1] = end
    return [(start, end) for start, end in merged]


def interval_bp(intervals: Iterable[Tuple[int, int]]) -> int:
    return sum(end - start for start, end in intervals)


def intersect_bp(a: Sequence[Tuple[int, int]], b: Sequence[Tuple[int, int]]) -> int:
    i = j = total = 0
    while i < len(a) and j < len(b):
        total += max(0, min(a[i][1], b[j][1]) - max(a[i][0], b[j][0]))
        if a[i][1] <= b[j][1]:
            i += 1
        else:
            j += 1
    return total


def subtract_intervals(a, b):
    result = []
    j = 0
    for start, end in a:
        cursor = start
        while j < len(b) and b[j][1] <= start:
            j += 1
        k = j
        while k < len(b) and b[k][0] < end:
            if b[k][0] > cursor:
                result.append((cursor, min(end, b[k][0])))
            cursor = max(cursor, b[k][1])
            if cursor >= end:
                break
            k += 1
        if cursor < end:
            result.append((cursor, end))
    return result


def flatten_best(hits: Sequence[Hit]) -> List[Hit]:
    """Resolve overlapping RepeatMasker records using the highest score."""
    result = []
    by_chrom = defaultdict(list)
    for hit in hits:
        by_chrom[hit.chrom].append(hit)
    for chrom, chrom_hits in by_chrom.items():
        events = defaultdict(lambda: [[], []])
        for index, hit in enumerate(chrom_hits):
            events[hit.start][1].append(index)
            events[hit.end][0].append(index)
        active = set()
        previous = None
        for position in sorted(events):
            if previous is not None and previous < position and active:
                winner = max(active, key=lambda i: (chrom_hits[i].score, chrom_hits[i].end - chrom_hits[i].start))
                hit = chrom_hits[winner]
                piece = Hit(chrom, previous, position, hit.family, hit.te_class,
                            hit.score, hit.source_id, hit.node_id)
                if (result and result[-1].chrom == chrom and result[-1].end == previous
                        and result[-1].family == piece.family
                        and result[-1].te_class == piece.te_class
                        and result[-1].source_id == piece.source_id):
                    old = result[-1]
                    result[-1] = Hit(chrom, old.start, position, old.family, old.te_class,
                                     old.score, old.source_id, old.node_id)
                else:
                    result.append(piece)
            removals, additions = events[position]
            active.difference_update(removals)
            active.update(additions)
            previous = position
    return sorted(result, key=lambda h: (h.chrom, h.start, h.end))


def labeled_agreement(a: Sequence[Hit], b: Sequence[Hit]):
    by_a, by_b = defaultdict(list), defaultdict(list)
    for hit in a:
        by_a[hit.chrom].append(hit)
    for hit in b:
        by_b[hit.chrom].append(hit)
    totals = defaultdict(int)
    class_pairs, family_pairs = defaultdict(int), defaultdict(int)
    for chrom in sorted(set(by_a) | set(by_b)):
        aa, bb = by_a[chrom], by_b[chrom]
        i = j = 0
        while i < len(aa) and j < len(bb):
            overlap = max(0, min(aa[i].end, bb[j].end) - max(aa[i].start, bb[j].start))
            if overlap:
                totals["overlap_bp"] += overlap
                if aa[i].te_class == bb[j].te_class:
                    totals["same_class_bp"] += overlap
                if aa[i].family == bb[j].family:
                    totals["same_family_bp"] += overlap
                class_pairs[(aa[i].te_class, bb[j].te_class)] += overlap
                family_pairs[(aa[i].family, bb[j].family)] += overlap
            if aa[i].end <= bb[j].end:
                i += 1
            else:
                j += 1
    return dict(totals), class_pairs, family_pairs


def feature_comparison(gfa_hits: Sequence[Hit], direct_hits: Sequence[Hit],
                       node_boundaries=None):
    node_boundaries = node_boundaries or {}
    by_chrom = defaultdict(list)
    direct_by_chrom = defaultdict(list)
    for hit in gfa_hits:
        by_chrom[hit.chrom].append(hit)
    for hit in direct_hits:
        direct_by_chrom[hit.chrom].append(hit)
    rows = []
    for chrom in sorted(direct_by_chrom):
        candidates = by_chrom.get(chrom, [])
        boundaries = node_boundaries.get(chrom, [])
        active = []
        next_candidate = 0
        for direct in direct_by_chrom[chrom]:
            while (next_candidate < len(candidates)
                   and candidates[next_candidate].start < direct.end):
                active.append(candidates[next_candidate])
                next_candidate += 1
            active = [hit for hit in active if hit.end > direct.start]
            overlaps = [hit for hit in active if hit.start < direct.end]
            overlap_bp = interval_bp(merge_intervals(
                (max(hit.start, direct.start), min(hit.end, direct.end)) for hit in overlaps
            ))
            same_class = [hit for hit in overlaps if hit.te_class == direct.te_class]
            same_family = [hit for hit in overlaps if hit.family == direct.family]
            same_class_bp = interval_bp(merge_intervals(
                (max(hit.start, direct.start), min(hit.end, direct.end)) for hit in same_class
            ))
            same_family_bp = interval_bp(merge_intervals(
                (max(hit.start, direct.start), min(hit.end, direct.end)) for hit in same_family
            ))
            best = max(same_family, key=lambda hit: max(
                0, min(hit.end, direct.end) - max(hit.start, direct.start)
            ), default=None)
            nodes_crossed = 1 + (
                bisect.bisect_left(boundaries, direct.end)
                - bisect.bisect_right(boundaries, direct.start))
            rows.append({
                "chrom": direct.chrom, "start": direct.start, "end": direct.end,
                "family": direct.family, "class": direct.te_class,
                "direct_length": direct.end - direct.start,
                "any_overlap_bp": overlap_bp, "same_class_bp": same_class_bp,
                "same_family_bp": same_family_bp, "overlapping_gfa_fragments": len(overlaps),
                "reference_nodes_crossed": nodes_crossed,
                "same_family_gfa_fragments": len(same_family),
                "best_gfa_start": "" if best is None else best.start,
                "best_gfa_end": "" if best is None else best.end,
                "best_gfa_node": "" if best is None else best.node_id,
            })
    return rows


def write_hits(path: Path, hits: Sequence[Hit]):
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["chrom", "start", "end", "family", "class", "score", "source_id", "node_id"])
        for hit in hits:
            writer.writerow([hit.chrom, hit.start, hit.end, hit.family, hit.te_class,
                             hit.score, hit.source_id, hit.node_id])


def write_bed(path: Path, intervals_by_chrom, label: str):
    with path.open("w") as handle:
        for chrom in sorted(intervals_by_chrom):
            for start, end in intervals_by_chrom[chrom]:
                handle.write(f"{chrom}\t{start}\t{end}\t{label}\t{end-start}\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gfa", required=True, type=Path)
    parser.add_argument("--reference-fasta", required=True, type=Path)
    parser.add_argument("--node-hits", required=True, type=Path)
    parser.add_argument("--direct-out", required=True, type=Path)
    parser.add_argument("--repeat-library", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    reference = read_fasta(args.reference_fasta)
    library = read_library(args.repeat_library)
    nodes, orientations, node_intervals = reference_nodes(args.gfa, reference)
    gfa_hits = parse_projected_hits(args.node_hits, nodes)
    direct_hits, direct_parse = parse_repeatmasker_out(
        args.direct_out, library, reference
    )
    write_hits(args.out_dir / "gfa_projected_reference_hits.tsv", gfa_hits)
    write_hits(args.out_dir / "direct_reference_hits.tsv", direct_hits)

    gfa_union, direct_union = {}, {}
    for chrom in reference:
        gfa_union[chrom] = merge_intervals(
            (hit.start, hit.end) for hit in gfa_hits if hit.chrom == chrom
        )
        direct_union[chrom] = merge_intervals(
            (hit.start, hit.end) for hit in direct_hits if hit.chrom == chrom
        )
    coverage_rows = []
    for chrom in list(reference) + ["ALL"]:
        chromosomes = list(reference) if chrom == "ALL" else [chrom]
        genome_bp = sum(len(reference[name]) for name in chromosomes)
        gfa_bp = sum(interval_bp(gfa_union[name]) for name in chromosomes)
        direct_bp = sum(interval_bp(direct_union[name]) for name in chromosomes)
        common_bp = sum(intersect_bp(gfa_union[name], direct_union[name]) for name in chromosomes)
        union_bp = gfa_bp + direct_bp - common_bp
        coverage_rows.append({
            "chrom": chrom, "genome_bp": genome_bp, "gfa_te_bp": gfa_bp,
            "direct_te_bp": direct_bp, "intersection_bp": common_bp,
            "gfa_only_bp": gfa_bp - common_bp,
            "direct_only_bp": direct_bp - common_bp,
            "gfa_vs_direct_recall": common_bp / direct_bp if direct_bp else 1.0,
            "direct_vs_gfa_recall": common_bp / gfa_bp if gfa_bp else 1.0,
            "jaccard": common_bp / union_bp if union_bp else 1.0,
        })
    with (args.out_dir / "coverage_summary.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", lineterminator="\n",
                                fieldnames=list(coverage_rows[0]))
        writer.writeheader(); writer.writerows(coverage_rows)

    direct_only, gfa_only = {}, {}
    for chrom in reference:
        direct_only[chrom] = subtract_intervals(direct_union[chrom], gfa_union[chrom])
        gfa_only[chrom] = subtract_intervals(gfa_union[chrom], direct_union[chrom])
    write_bed(args.out_dir / "direct_only.bed", direct_only, "direct_only")
    write_bed(args.out_dir / "gfa_only.bed", gfa_only, "gfa_only")

    label_totals, class_pairs, family_pairs = labeled_agreement(
        flatten_best(gfa_hits), flatten_best(direct_hits)
    )
    for name, pairs in (("class_confusion.tsv", class_pairs), ("family_confusion.tsv", family_pairs)):
        with (args.out_dir / name).open("w") as handle:
            handle.write("gfa_label\tdirect_label\toverlap_bp\n")
            for (left, right), bp in sorted(pairs.items(), key=lambda item: -item[1]):
                handle.write(f"{left}\t{right}\t{bp}\n")

    node_boundaries = defaultdict(list)
    for chrom, offset, _length, _orientation in nodes.values():
        if offset:
            node_boundaries[chrom].append(offset)
    for boundaries in node_boundaries.values():
        boundaries.sort()
    feature_rows = feature_comparison(gfa_hits, direct_hits, node_boundaries)
    with (args.out_dir / "direct_hit_comparison.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", lineterminator="\n",
                                fieldnames=list(feature_rows[0]) if feature_rows else [])
        if feature_rows:
            writer.writeheader(); writer.writerows(feature_rows)

    fragment_counts = [row["same_family_gfa_fragments"] for row in feature_rows]
    direct_lengths = [row["direct_length"] for row in feature_rows]
    crossing_rows = [row for row in feature_rows if row["reference_nodes_crossed"] > 1]
    contained_rows = [row for row in feature_rows if row["reference_nodes_crossed"] == 1]
    same_family_half = lambda row: row["same_family_bp"] * 2 >= row["direct_length"]
    feature_summary = {
        "direct_hits": len(feature_rows),
        "gfa_projected_hits": len(gfa_hits),
        "direct_hits_with_any_overlap": sum(row["any_overlap_bp"] > 0 for row in feature_rows),
        "direct_hits_with_same_class_overlap": sum(row["same_class_bp"] > 0 for row in feature_rows),
        "direct_hits_with_same_family_overlap": sum(row["same_family_bp"] > 0 for row in feature_rows),
        "direct_hits_with_same_family_coverage_ge_50pct": sum(
            row["same_family_bp"] * 2 >= row["direct_length"] for row in feature_rows
        ),
        "direct_hits_split_across_multiple_same_family_gfa_fragments": sum(
            count > 1 for count in fragment_counts
        ),
        "direct_hits_crossing_node_boundary": len(crossing_rows),
        "direct_hits_contained_in_one_node": len(contained_rows),
        "crossing_hits_same_family_coverage_ge_50pct": sum(
            same_family_half(row) for row in crossing_rows),
        "contained_hits_same_family_coverage_ge_50pct": sum(
            same_family_half(row) for row in contained_rows),
        "crossing_hits_same_family_coverage_ge_50pct_fraction": (
            sum(same_family_half(row) for row in crossing_rows) / len(crossing_rows)
            if crossing_rows else 1.0),
        "contained_hits_same_family_coverage_ge_50pct_fraction": (
            sum(same_family_half(row) for row in contained_rows) / len(contained_rows)
            if contained_rows else 1.0),
        "median_same_family_gfa_fragments_per_direct_hit": statistics.median(fragment_counts) if fragment_counts else 0,
        "median_direct_hit_length": statistics.median(direct_lengths) if direct_lengths else 0,
    }
    all_coverage = coverage_rows[-1]
    overlap_bp = label_totals.get("overlap_bp", 0)
    summary = {
        "inputs": {key: str(value) for key, value in vars(args).items()},
        "reference": {
            "sequences": {name: len(seq) for name, seq in reference.items()},
            "sr0_nodes": len(nodes), "node_orientations": orientations,
            "sr0_interval_bp": sum(interval_bp(merge_intervals(value)) for value in node_intervals.values()),
        },
        "direct_repeatmasker_parse": direct_parse,
        "coverage": all_coverage,
        "labels": {
            **label_totals,
            "same_class_fraction_of_overlap": label_totals.get("same_class_bp", 0) / overlap_bp if overlap_bp else 1.0,
            "same_family_fraction_of_overlap": label_totals.get("same_family_bp", 0) / overlap_bp if overlap_bp else 1.0,
        },
        "features": feature_summary,
        "parameter_caveat": (
            "RepeatMasker sensitivity/search options must match between runs; filtering "
            "to custom-library families removes default simple repeats but cannot correct "
            "a sensitivity-mode mismatch."
        ),
    }
    with (args.out_dir / "comparison_summary.json").open("w") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
