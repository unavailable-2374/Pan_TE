#!/usr/bin/env python3
"""Strict rGFA + GAF support for locus-aware Pan_TE enrichment.

rGFA does not encode complete sample paths.  This module recovers only the
complete GFA nodes visited by primary GAF records, uses ordered SR:i:0 nodes as
reference anchors, and emits graph-path allele candidates plus conservative
0/1/NA locus calls.  It deliberately never turns a missing GAF traversal into
an absence call.
"""
from __future__ import annotations

from bisect import bisect_left
from collections import Counter, defaultdict
from dataclasses import dataclass
import csv
import hashlib
import os
import re
import sqlite3
from typing import Callable, Dict, Iterable, Iterator, List, Optional, Sequence, Set, Tuple

from gfa_alleles import AlleleCandidate, _classify_allele, build_sequence, window_starts


_WALK_RE = re.compile(r"([<>])([^<>]+)")
_COORD_RE = re.compile(r"^(.*):(\d+)-(\d+)$")
_CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
_COMPLEMENT = str.maketrans(
    "ACGTacgtRYSWKMBDHVNryswkmbdhvn",
    "TGCAtgcaYRSWMKVHDBNryswmkvhdbn",
)


class RgfaGafError(ValueError):
    """Raised when an rGFA/GAF input cannot support safe path recovery."""


@dataclass(frozen=True)
class Segment:
    node_id: str
    sequence: str
    length: int
    rank: int
    source_name: str
    source_offset: int


@dataclass(frozen=True)
class GafRecord:
    record_id: int
    query: str
    query_length: int
    query_start: int
    query_end: int
    strand: str
    path: str
    path_length: int
    path_start: int
    path_end: int
    matches: int
    block_length: int
    mapq: int
    tags: Dict[str, str]

    @property
    def identity(self) -> float:
        return self.matches / self.block_length if self.block_length else 0.0


@dataclass(frozen=True)
class GafStep:
    node_id: str
    orientation: str


@dataclass(frozen=True)
class ReferenceContig:
    name: str
    node_ids: Tuple[str, ...]
    prefix: Tuple[int, ...]
    rank_by_node: Dict[str, int]


@dataclass(frozen=True)
class LocusCall:
    locus_key: str
    contig: str
    left_anchor: str
    right_anchor: str
    ref_span_bp: int
    sample: str
    sequence_sha256: str
    is_candidate: bool


def reverse_complement(sequence: str) -> str:
    return sequence.translate(_COMPLEMENT)[::-1]


def sample_from_query(name: str) -> str:
    return name.split("#", 1)[0] if "#" in name else name


def contig_from_query(name: str) -> str:
    return name.split("#")[-1] if "#" in name else name


def _parse_tags(fields: Sequence[str]) -> Dict[str, str]:
    result = {}
    for field in fields:
        pieces = field.split(":", 2)
        if len(pieces) == 3:
            result[pieces[0]] = pieces[2]
    return result


def read_fasta_lengths(path: str) -> Dict[str, int]:
    lengths: Dict[str, int] = {}
    name: Optional[str] = None
    length = 0
    with open(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    if name in lengths:
                        raise RgfaGafError(f"duplicate FASTA query name: {name!r}")
                    lengths[name] = length
                name = line[1:].split()[0]
                length = 0
            elif name is None:
                raise RgfaGafError(f"sequence before FASTA header in {path}")
            else:
                length += len(line)
    if name is not None:
        if name in lengths:
            raise RgfaGafError(f"duplicate FASTA query name: {name!r}")
        lengths[name] = length
    if not lengths:
        raise RgfaGafError(f"no FASTA records found in {path}")
    return lengths


def load_rgfa_graph(gfa_file: str) -> Tuple[Dict[str, Segment], Dict[Tuple[str, str, str, str], str], Dict[Tuple[str, int, int], str], Dict[str, ReferenceContig], Dict[str, List[Tuple[int, int, str]]]]:
    """Load rGFA nodes, zero-overlap links, coordinate tokens and SR0 anchors."""
    segments: Dict[str, Segment] = {}
    links: Dict[Tuple[str, str, str, str], str] = {}
    coordinates: Dict[Tuple[str, int, int], str] = {}
    reference_by_contig: Dict[str, List[Tuple[int, str]]] = defaultdict(list)
    source_intervals: Dict[str, List[Tuple[int, int, str]]] = defaultdict(list)

    with open(gfa_file) as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if fields[0] == "S":
                if len(fields) < 3:
                    raise RgfaGafError(f"malformed S record at GFA line {line_number}")
                node_id, sequence = fields[1], fields[2]
                if node_id in segments:
                    raise RgfaGafError(f"duplicate GFA segment ID: {node_id!r}")
                if not sequence or sequence == "*":
                    raise RgfaGafError(
                        f"rGFA+GAF requires embedded sequences; node {node_id!r} has none")
                tags = _parse_tags(fields[3:])
                if "SR" not in tags or "SN" not in tags or "SO" not in tags:
                    raise RgfaGafError(
                        f"rGFA node {node_id!r} lacks required SR/SN/SO tags")
                try:
                    rank = int(tags["SR"])
                    offset = int(tags["SO"])
                    length = int(tags.get("LN", len(sequence)))
                except ValueError as exc:
                    raise RgfaGafError(f"invalid rGFA tags at line {line_number}") from exc
                if length != len(sequence):
                    raise RgfaGafError(
                        f"LN differs from embedded sequence length for node {node_id!r}")
                source_name = tags["SN"]
                segment = Segment(node_id, sequence, length, rank, source_name, offset)
                segments[node_id] = segment
                key = (source_name, offset, offset + length)
                if key in coordinates:
                    raise RgfaGafError(f"duplicate rGFA coordinate interval: {key!r}")
                coordinates[key] = node_id
                source_intervals[source_name].append((offset, offset + length, node_id))
                if rank == 0:
                    reference_by_contig[contig_from_query(source_name)].append((offset, node_id))
            elif fields[0] == "L":
                if len(fields) < 6:
                    raise RgfaGafError(f"malformed L record at GFA line {line_number}")
                left, left_orient, right, right_orient, overlap = fields[1:6]
                if overlap != "0M":
                    raise RgfaGafError(
                        "rGFA+GAF currently requires 0M links for exact node "
                        f"coordinates; found {overlap!r} at GFA line {line_number}")
                links[(left, left_orient, right, right_orient)] = overlap
                reverse = (
                    right, "+" if right_orient == "-" else "-",
                    left, "+" if left_orient == "-" else "-",
                )
                links[reverse] = overlap

    if not segments:
        raise RgfaGafError("GFA contains no S records")
    if not reference_by_contig:
        raise RgfaGafError("rGFA contains no SR:i:0 reference nodes")

    references: Dict[str, ReferenceContig] = {}
    for contig, entries in reference_by_contig.items():
        entries.sort()
        previous_end = -1
        ids: List[str] = []
        prefix = [0]
        for offset, node_id in entries:
            node = segments[node_id]
            if offset < previous_end:
                raise RgfaGafError(
                    f"overlapping SR0 reference coordinates on {contig}: {node_id!r}")
            if offset != previous_end and previous_end >= 0:
                raise RgfaGafError(
                    f"gapped SR0 reference coordinates on {contig}: "
                    f"{previous_end}..{offset}")
            previous_end = offset + node.length
            ids.append(node_id)
            prefix.append(prefix[-1] + node.length)
        references[contig] = ReferenceContig(
            contig, tuple(ids), tuple(prefix), {node: i for i, node in enumerate(ids)})
    for source_name in source_intervals:
        source_intervals[source_name].sort()
    return segments, links, coordinates, references, dict(source_intervals)


def _parse_cigar_spans(cigar: str) -> Tuple[int, int]:
    matches = list(_CIGAR_RE.finditer(cigar))
    if not matches or "".join(item.group(0) for item in matches) != cigar:
        raise RgfaGafError(f"malformed cg tag: {cigar!r}")
    query = target = 0
    for match in matches:
        length, op = int(match.group(1)), match.group(2)
        if op in "M=XIS":
            query += length
        if op in "M=XDN":
            target += length
    return query, target


def _resolve_steps(token: str, orientation: str,
                   coordinates: Dict[Tuple[str, int, int], str],
                   source_intervals: Dict[str, List[Tuple[int, int, str]]],
                   segments: Dict[str, Segment]) -> List[GafStep]:
    """Expand one minigraph coordinate token to its exact consecutive S nodes."""
    if token in segments:
        return [GafStep(token, orientation)]
    match = _COORD_RE.match(token)
    if match is None:
        raise RgfaGafError(f"GAF path token cannot be resolved to an rGFA node: {token!r}")
    source, start, end = match.group(1), int(match.group(2)), int(match.group(3))
    exact = coordinates.get((source, start, end))
    if exact is not None:
        return [GafStep(exact, orientation)]
    intervals = source_intervals.get(source, [])
    expanded: List[str] = []
    cursor = start
    index = bisect_left(intervals, (start, -1, ''))
    for node_start, node_end, node_id in intervals[index:]:
        if node_start != cursor or node_end > end:
            break
        expanded.append(node_id)
        cursor = node_end
        if cursor == end:
            break
    if cursor != end or not expanded:
        raise RgfaGafError(
            f"GAF coordinate token does not expand to complete consecutive rGFA nodes: {token!r}")
    if orientation == "+":
        return [GafStep(node_id, "+") for node_id in expanded]
    return [GafStep(node_id, "-") for node_id in reversed(expanded)]


def _parse_walk(path: str, coordinates: Dict[Tuple[str, int, int], str],
                source_intervals: Dict[str, List[Tuple[int, int, str]]],
                segments: Dict[str, Segment]) -> List[GafStep]:
    if path.startswith((">", "<")):
        raw_steps = _WALK_RE.findall(path)
        if not raw_steps or "".join(direction + token for direction, token in raw_steps) != path:
            raise RgfaGafError(f"malformed GAF walk: {path!r}")
        result: List[GafStep] = []
        for direction, token in raw_steps:
            result.extend(_resolve_steps(
                token, "+" if direction == ">" else "-", coordinates,
                source_intervals, segments))
        return result
    if path in source_intervals:
        # minigraph may encode a source contig as a named path and use the
        # GAF path_start/path_end fields to select its aligned subinterval.
        return [GafStep(node_id, "+")
                for _start, _end, node_id in source_intervals[path]]
    return _resolve_steps(path, "+", coordinates, source_intervals, segments)


def iter_gaf_records(
        gaf_file: str,
        segments: Dict[str, Segment],
        links: Dict[Tuple[str, str, str, str], str],
        coordinates: Dict[Tuple[str, int, int], str],
        source_intervals: Dict[str, List[Tuple[int, int, str]]],
        query_lengths: Optional[Dict[str, int]] = None,
        require_primary: bool = True,
        require_cigar: bool = True,
        sample_sink: Optional[Callable[[str], None]] = None,
        stats: Optional[Counter] = None,
) -> Iterator[Tuple[GafRecord, List[GafStep]]]:
    """Yield structurally valid GAF records and exact graph-node walks."""
    counters = stats if stats is not None else Counter()
    with open(gaf_file) as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 12:
                raise RgfaGafError(f"malformed GAF record at line {line_number}")
            try:
                record = GafRecord(
                    line_number, fields[0], int(fields[1]), int(fields[2]), int(fields[3]),
                    fields[4], fields[5], int(fields[6]), int(fields[7]), int(fields[8]),
                    int(fields[9]), int(fields[10]), int(fields[11]), _parse_tags(fields[12:]))
            except ValueError as exc:
                raise RgfaGafError(f"invalid numeric GAF field at line {line_number}") from exc
            counters["records_total"] += 1
            if record.strand not in {"+", "-"}:
                raise RgfaGafError(f"invalid GAF strand at line {line_number}: {record.strand!r}")
            if not (0 <= record.query_start <= record.query_end <= record.query_length):
                raise RgfaGafError(f"invalid GAF query interval at line {line_number}")
            if sample_sink is not None:
                sample_sink(sample_from_query(record.query))
            if not (0 <= record.path_start <= record.path_end <= record.path_length):
                raise RgfaGafError(f"invalid GAF path interval at line {line_number}")
            if query_lengths is not None:
                length = query_lengths.get(record.query)
                if length is None:
                    raise RgfaGafError(f"GAF query absent from --gfa-gaf-fasta: {record.query!r}")
                if length != record.query_length:
                    raise RgfaGafError(
                        f"GAF query length differs from FASTA for {record.query!r}: "
                        f"{record.query_length} != {length}")
            if require_primary and record.tags.get("tp") != "P":
                counters["records_skipped_nonprimary"] += 1
                continue
            cigar = record.tags.get("cg")
            if require_cigar and not cigar:
                raise RgfaGafError(f"primary GAF record lacks cg tag at line {line_number}")
            if cigar:
                qspan, pspan = _parse_cigar_spans(cigar)
                if qspan != record.query_end - record.query_start:
                    raise RgfaGafError(f"cg query span disagrees with GAF fields at line {line_number}")
                if pspan != record.path_end - record.path_start:
                    raise RgfaGafError(f"cg path span disagrees with GAF fields at line {line_number}")
            if record.block_length <= 0 or record.matches < 0 or record.matches > record.block_length:
                raise RgfaGafError(f"invalid GAF match/block lengths at line {line_number}")
            steps = _parse_walk(record.path, coordinates, source_intervals, segments)
            node_length = sum(segments[step.node_id].length for step in steps)
            if node_length != record.path_length:
                raise RgfaGafError(
                    f"GAF path length differs from resolved nodes at line {line_number}: "
                    f"{record.path_length} != {node_length}")
            for left, right in zip(steps, steps[1:]):
                if (left.node_id, left.orientation, right.node_id, right.orientation) not in links:
                    raise RgfaGafError(
                        f"GAF walk edge is absent from GFA at line {line_number}: "
                        f"{left.node_id}{left.orientation}->{right.node_id}{right.orientation}")
            counters["records_valid"] += 1
            counters[f"mapq_{record.mapq}"] += 1
            yield record, steps


def _complete_steps(record: GafRecord, steps: Sequence[GafStep], segments: Dict[str, Segment]) -> Tuple[List[GafStep], int]:
    """Return only nodes fully covered by the aligned path interval."""
    offset = 0
    complete: List[GafStep] = []
    partial = 0
    for step in steps:
        end = offset + segments[step.node_id].length
        if offset >= record.path_start and end <= record.path_end:
            complete.append(step)
        elif end > record.path_start and offset < record.path_end:
            partial += 1
        offset = end
    return complete, partial


def _normalize_to_query_forward(record: GafRecord, steps: Sequence[GafStep]) -> List[GafStep]:
    if record.strand == "+":
        return list(steps)
    return [GafStep(step.node_id, "+" if step.orientation == "-" else "-")
            for step in reversed(steps)]


def _maximum_collinear_anchors(
        steps: Sequence[GafStep], reference: ReferenceContig, segments: Dict[str, Segment]
) -> List[Tuple[int, int, str]]:
    """Return an increasing chain of complete, forward SR0 anchors."""
    tails: List[int] = []
    tail_trace: List[int] = []
    trace: List[Tuple[int, int, str, int]] = []
    for path_index, step in enumerate(steps):
        if step.orientation != "+":
            continue
        segment = segments[step.node_id]
        if segment.rank != 0 or contig_from_query(segment.source_name) != reference.name:
            continue
        rank = reference.rank_by_node.get(step.node_id)
        if rank is None:
            continue
        chain_rank = bisect_left(tails, rank)
        previous = tail_trace[chain_rank - 1] if chain_rank else -1
        trace_index = len(trace)
        trace.append((path_index, rank, step.node_id, previous))
        if chain_rank == len(tails):
            tails.append(rank)
            tail_trace.append(trace_index)
        else:
            tails[chain_rank] = rank
            tail_trace[chain_rank] = trace_index
    if not tail_trace:
        return []
    result: List[Tuple[int, int, str]] = []
    trace_index = tail_trace[-1]
    while trace_index >= 0:
        path_index, rank, node_id, previous = trace[trace_index]
        result.append((path_index, rank, node_id))
        trace_index = previous
    return list(reversed(result))


def _reference_internal_sequence(reference: ReferenceContig, left_rank: int, right_rank: int,
                                 segments: Dict[str, Segment]) -> Tuple[str, int]:
    if right_rank <= left_rank:
        raise RgfaGafError("non-increasing reference anchor ranks after collinear chaining")
    ids = reference.node_ids[left_rank + 1:right_rank]
    return "".join(segments[node_id].sequence for node_id in ids), (
        reference.prefix[right_rank] - reference.prefix[left_rank + 1])


def _sha256(sequence: str) -> str:
    return hashlib.sha256(sequence.upper().encode()).hexdigest()


def _reference_samples(segments: Dict[str, Segment]) -> Set[str]:
    """Return source-sample names that define the SR:i:0 backbone."""
    return {sample_from_query(segment.source_name) for segment in segments.values()
            if segment.rank == 0}


def extract_rgfa_gaf_alleles(
        gfa_file: str,
        gaf_file: str,
        min_size: int = 100,
        window_size: int = 30000,
        window_overlap: int = 1000,
        min_mapq: int = 20,
        min_identity: float = 0.90,
        query_fasta: Optional[str] = None,
        candidate_sink: Optional[Callable[[AlleleCandidate], None]] = None,
        call_sink: Optional[Callable[[LocusCall], None]] = None,
        sample_sink: Optional[Callable[[str], None]] = None,
        stats: Optional[dict] = None,
) -> List[AlleleCandidate]:
    """Stream primary GAF records into complete-node rGFA allele candidates.

    The emitted sequence is deliberately labelled ``graph_complete_nodes``.
    It is not a CIGAR replay of the sample query and terminal partial nodes are
    excluded rather than approximated.
    """
    if min_size < 1:
        raise ValueError("min_size must be positive")
    if not 0 <= min_identity <= 1:
        raise ValueError("min_identity must be between 0 and 1")
    segments, links, coordinates, references, source_intervals = load_rgfa_graph(gfa_file)
    query_lengths = read_fasta_lengths(query_fasta) if query_fasta else None
    counters: Counter = Counter()
    records: List[AlleleCandidate] = []
    sequences = {node_id: segment.sequence for node_id, segment in segments.items()}
    reference_samples = _reference_samples(segments)
    counters["reference_samples"] = len(reference_samples)
    candidate_number = 0

    for record, raw_steps in iter_gaf_records(
            gaf_file, segments, links, coordinates, source_intervals, query_lengths,
            sample_sink=sample_sink, stats=counters):
        sample = sample_from_query(record.query)
        # SR:i:0 is the graph's declared reference truth. A mapper can choose
        # an equivalent non-reference route for this query, but that must not
        # manufacture a reference-sample alternate allele.
        if sample in reference_samples:
            counters["records_skipped_reference_query"] += 1
            continue
        contig = contig_from_query(record.query)
        reference = references.get(contig)
        if reference is None:
            counters["records_skipped_unmatched_contig"] += 1
            continue
        complete, partial = _complete_steps(record, raw_steps, segments)
        counters["partial_terminal_nodes_skipped"] += partial
        if record.mapq < min_mapq or record.identity < min_identity:
            counters["records_skipped_low_quality"] += 1
            continue
        steps = _normalize_to_query_forward(record, complete)
        anchors = _maximum_collinear_anchors(steps, reference, segments)
        counters["alternate_records_processed"] += 1
        counters["collinear_anchors"] += len(anchors)
        if len(anchors) < 2:
            counters["fragments_below_two_anchors"] += 1
            continue
        counters["fragments_with_two_anchors"] += 1

        for pair_index, (left, right) in enumerate(zip(anchors, anchors[1:])):
            left_path, left_rank, left_node = left
            right_path, right_rank, right_node = right
            if right_path <= left_path or right_rank <= left_rank:
                continue
            alt_steps = steps[left_path + 1:right_path]
            alt_sequence = build_sequence(
                [(step.node_id, step.orientation) for step in alt_steps],
                sequences)
            ref_sequence, ref_span = _reference_internal_sequence(
                reference, left_rank, right_rank, segments)
            locus_key = f"{reference.name}:{left_rank}:{right_rank}"
            changed = alt_sequence.upper() != ref_sequence.upper()
            is_candidate = changed and len(alt_sequence) >= min_size
            if call_sink is not None:
                call_sink(LocusCall(
                    locus_key, reference.name, left_node, right_node, ref_span,
                    sample, _sha256(alt_sequence), is_candidate))
            if not changed:
                counters["unchanged_sequence_intervals"] += 1
                continue
            if len(alt_sequence) < min_size:
                counters["alleles_below_min"] += 1
                counters["alleles_below_min_bp"] += len(alt_sequence)
                continue
            allele_type, delta = _classify_allele(ref_span, len(alt_sequence), min_size)
            parent_id = f"rgfa_gaf_parent_{record.record_id}_{pair_index + 1}"
            full_digest = _sha256(alt_sequence)
            starts = window_starts(len(alt_sequence), window_size, window_overlap, min_size)
            if len(starts) > 1:
                counters["long_alleles_windowed"] += 1
            for chunk_index, start in enumerate(starts, 1):
                chunk = alt_sequence[start:start + window_size]
                if len(chunk) < min_size:
                    raise RgfaGafError("windowing produced a sub-minimum candidate tail")
                candidate_number += 1
                candidate = AlleleCandidate(
                    candidate_id=f"gfa_rgfa_gaf_{candidate_number:08d}",
                    parent_allele_id=parent_id,
                    sequence=chunk,
                    sample=sample,
                    path=record.query,
                    contig=contig,
                    ref_path=reference.name,
                    path_step_start=left_path + 1,
                    path_step_end=right_path,
                    left_anchor=left_node,
                    right_anchor=right_node,
                    left_ref_rank=left_rank,
                    right_ref_rank=right_rank,
                    ref_span_bp=ref_span,
                    alt_span_bp=len(alt_sequence),
                    net_change_bp=delta,
                    allele_type=allele_type,
                    chunk_index=chunk_index,
                    chunk_start=start,
                    chunk_end=start + len(chunk),
                    locus_key=locus_key,
                    sequence_sha256=_sha256(chunk),
                    sequence_source="graph_complete_nodes",
                    source_query=record.query,
                    query_start=record.query_start,
                    query_end=record.query_end,
                    gaf_record_id=record.record_id,
                    mapq=record.mapq,
                    alignment_identity=record.identity,
                )
                if candidate_sink is None:
                    records.append(candidate)
                else:
                    candidate_sink(candidate)
                counters["candidate_occurrences"] += 1
                counters["candidate_bp"] += len(chunk)
                counters[f"allele_type_{allele_type}"] += 1
                counters["parent_sequence_hashes"] += 1 if chunk_index == 1 else 0
                counters["parent_full_digest_records"] += 1 if full_digest else 0

    if stats is not None:
        stats.update(counters)
    return records


def populate_rgfa_gaf_locus_calls(
        gfa_file: str, gaf_file: str, representatives: Iterable[AlleleCandidate],
        call_sink: Callable[[LocusCall], None], min_size: int = 100,
        min_mapq: int = 20, min_identity: float = 0.90,
        query_fasta: Optional[str] = None,
        sample_sink: Optional[Callable[[str], None]] = None,
        stats: Optional[dict] = None,
) -> None:
    """Populate calls for candidate loci crossed by each quality GAF record.

    Candidate extraction uses adjacent maximum-collinear anchors.  A different
    sample can legitimately have additional reference anchors between those
    flanks; this second streaming pass therefore calls a locus whenever the
    *outer* candidate anchors are both present in one complete, quality record.
    """
    segments, links, coordinates, references, source_intervals = load_rgfa_graph(gfa_file)
    query_lengths = read_fasta_lengths(query_fasta) if query_fasta else None
    targets: Dict[str, Dict[int, Set[int]]] = defaultdict(lambda: defaultdict(set))
    for record in representatives:
        if record.left_ref_rank is None or record.right_ref_rank is None:
            continue
        targets[record.ref_path][record.left_ref_rank].add(record.right_ref_rank)
    counters: Counter = Counter()
    sequences = {node_id: segment.sequence for node_id, segment in segments.items()}
    reference_samples = _reference_samples(segments)
    counters["reference_samples"] = len(reference_samples)
    for record, raw_steps in iter_gaf_records(
            gaf_file, segments, links, coordinates, source_intervals, query_lengths,
            sample_sink=sample_sink, stats=counters):
        if sample_from_query(record.query) in reference_samples:
            counters['records_skipped_reference_query'] += 1
            continue
        if record.mapq < min_mapq or record.identity < min_identity:
            counters['records_skipped_low_quality'] += 1
            continue
        contig = contig_from_query(record.query)
        reference = references.get(contig)
        if reference is None or reference.name not in targets:
            continue
        complete, partial = _complete_steps(record, raw_steps, segments)
        counters['partial_terminal_nodes_skipped'] += partial
        steps = _normalize_to_query_forward(record, complete)
        anchors = _maximum_collinear_anchors(steps, reference, segments)
        rank_to_path = {rank: path_index for path_index, rank, _node in anchors}
        sample = sample_from_query(record.query)
        for left_rank, right_ranks in targets[reference.name].items():
            left_path = rank_to_path.get(left_rank)
            if left_path is None:
                continue
            left_node = reference.node_ids[left_rank]
            for right_rank in right_ranks:
                right_path = rank_to_path.get(right_rank)
                if right_path is None or right_path <= left_path:
                    continue
                right_node = reference.node_ids[right_rank]
                alt_steps = steps[left_path + 1:right_path]
                alt_sequence = build_sequence(
                    [(step.node_id, step.orientation) for step in alt_steps],
                    sequences)
                ref_sequence, ref_span = _reference_internal_sequence(
                    reference, left_rank, right_rank, segments)
                changed = alt_sequence.upper() != ref_sequence.upper()
                call_sink(LocusCall(
                    f'{reference.name}:{left_rank}:{right_rank}', reference.name,
                    left_node, right_node, ref_span, sample,
                    _sha256(alt_sequence), changed and len(alt_sequence) >= min_size))
                counters['locus_calls'] += 1
    if stats is not None:
        stats.update(counters)


class RgfaGafCallStore:
    """Disk-backed conservative 0/1/NA locus calls for rGFA GAF records."""

    def __init__(self, path: str):
        self.path = path
        if os.path.exists(path):
            os.remove(path)
        self.connection = sqlite3.connect(path)
        self.connection.execute("PRAGMA journal_mode=OFF")
        self.connection.execute("PRAGMA synchronous=OFF")
        self.connection.executescript("""
            CREATE TABLE calls (
                locus_key TEXT NOT NULL,
                sample TEXT NOT NULL,
                sequence_sha256 TEXT NOT NULL,
                is_candidate INTEGER NOT NULL,
                PRIMARY KEY (locus_key, sample, sequence_sha256, is_candidate)
            );
            CREATE TABLE loci (
                locus_key TEXT PRIMARY KEY,
                contig TEXT NOT NULL,
                left_anchor TEXT NOT NULL,
                right_anchor TEXT NOT NULL,
                ref_span_bp INTEGER NOT NULL
            );
            CREATE TABLE samples (sample TEXT PRIMARY KEY);
        """)

    def add(self, call: LocusCall) -> None:
        self.register_sample(call.sample)
        self.connection.execute(
            "INSERT OR IGNORE INTO loci VALUES (?,?,?,?,?)",
            (call.locus_key, call.contig, call.left_anchor, call.right_anchor,
             call.ref_span_bp))
        self.connection.execute(
            "INSERT OR IGNORE INTO calls VALUES (?,?,?,?)",
            (call.locus_key, call.sample, call.sequence_sha256, int(call.is_candidate)))

    def register_sample(self, sample: str) -> None:
        """Record every GAF sample even if it produces no locus call."""
        self.connection.execute(
            "INSERT OR IGNORE INTO samples(sample) VALUES (?)", (sample,))

    def write_presence_matrix(self, path: str, locus_keys: Iterable[str],
                              reference_samples: Iterable[str] = ()) -> List[str]:
        self.connection.commit()
        reference_samples = set(reference_samples)
        samples = sorted(reference_samples | {row[0] for row in self.connection.execute(
            "SELECT sample FROM samples")})
        keys = sorted(set(locus_keys))
        with open(path, "w", newline="") as handle:
            fields = ["locus_key", "contig", "left_anchor", "right_anchor", "ref_span_bp"] + [
                f"sample:{sample}" for sample in samples]
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for key in keys:
                meta = self.connection.execute(
                    "SELECT contig,left_anchor,right_anchor,ref_span_bp FROM loci WHERE locus_key=?", (key,)).fetchone()
                if meta is None:
                    continue
                row = {"locus_key": key, "contig": meta[0], "left_anchor": meta[1],
                       "right_anchor": meta[2], "ref_span_bp": meta[3]}
                values = defaultdict(set)
                for sample, digest, candidate in self.connection.execute(
                        "SELECT sample,sequence_sha256,is_candidate FROM calls WHERE locus_key=?", (key,)):
                    values[sample].add((digest, candidate))
                for sample in samples:
                    if sample in reference_samples:
                        row[f"sample:{sample}"] = "0"
                        continue
                    observed = values.get(sample, set())
                    if not observed:
                        row[f"sample:{sample}"] = "NA"
                    elif len({digest for digest, _ in observed}) != 1:
                        row[f"sample:{sample}"] = "NA"
                    elif any(candidate for _, candidate in observed):
                        row[f"sample:{sample}"] = "1"
                    else:
                        row[f"sample:{sample}"] = "0"
                writer.writerow(row)
        return samples

    def close(self, remove: bool = False) -> None:
        if self.connection is not None:
            self.connection.close()
            self.connection = None
        if remove and os.path.exists(self.path):
            os.remove(self.path)


def derive_rgfa_node_presence(
        gfa_file: str, gaf_file: str, min_mapq: int = 20,
        min_identity: float = 0.90, query_fasta: Optional[str] = None,
) -> Tuple[Set[str], Dict[Tuple[str, str], str], dict]:
    """Return conservative node 1/NA evidence from complete high-quality walks."""
    segments, links, coordinates, _references, source_intervals = load_rgfa_graph(gfa_file)
    query_lengths = read_fasta_lengths(query_fasta) if query_fasta else None
    counters: Counter = Counter()
    samples: Set[str] = set()
    states: Dict[Tuple[str, str], str] = {}
    for record, steps in iter_gaf_records(
            gaf_file, segments, links, coordinates, source_intervals, query_lengths,
            sample_sink=samples.add, stats=counters):
        sample = sample_from_query(record.query)
        samples.add(sample)
        if record.mapq < min_mapq or record.identity < min_identity:
            counters["records_skipped_low_quality"] += 1
            continue
        complete, partial = _complete_steps(record, steps, segments)
        counters["partial_terminal_nodes_skipped"] += partial
        for step in complete:
            states[(step.node_id, sample)] = "1"
    return samples, states, dict(counters)


__all__ = [
    "RgfaGafError", "RgfaGafCallStore", "extract_rgfa_gaf_alleles",
    "populate_rgfa_gaf_locus_calls",
    "derive_rgfa_node_presence", "load_rgfa_graph", "read_fasta_lengths",
    "GafRecord", "GafStep", "Segment", "iter_gaf_records",
]
