#!/usr/bin/env python3
"""Annotate GFA segments with a Pan_TE library and sample presence.

RepeatMasker is run on sequence contexts reconstructed from GFA P/W records.
Every base in a context is tracked back to its source segment, allowing repeat
hits to be projected onto segment-local coordinates.  Segment coverage is
computed as an interval union, so a segment may be partially TE and repeated
occurrences of the same segment do not inflate its annotated fraction.
"""

from __future__ import annotations

import bisect
import logging
import os
import re
import shutil
import subprocess
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

from rgfa_gaf import derive_rgfa_node_presence


logger = logging.getLogger(__name__)

_COMPLEMENT = str.maketrans(
    "ACGTacgtRYSWKMBDHVNryswkmbdhvn",
    "TGCAtgcaYRSWMKVHDBNyrswmkvhdbn",
)
_WALK_STEP_RE = re.compile(r"([><])([^><]+)")
_CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(_COMPLEMENT)[::-1]


@dataclass(frozen=True)
class Node:
    node_id: str
    sequence: Optional[str]
    length: int


@dataclass(frozen=True)
class Traversal:
    traversal_id: str
    sample: str
    steps: Tuple[Tuple[str, str], ...]
    overlaps: Tuple[Optional[str], ...]
    kind: str


@dataclass(frozen=True)
class Occurrence:
    query_id: str
    traversal_id: str
    sample: str
    node_id: str
    query_start: int
    query_end: int
    node_start: int
    node_end: int
    orientation: str


@dataclass(frozen=True)
class RepeatHit:
    score: int
    divergence: float
    deletion: float
    insertion: float
    query_id: str
    query_start: int
    query_end: int
    strand: str
    repeat_name: str
    repeat_class: str


def _sample_from_path_name(path_name: str) -> str:
    """Use the PanSN sample component when present, otherwise the path name."""
    return path_name.split("#", 1)[0] if "#" in path_name else path_name


def _merge_intervals(intervals: Iterable[Tuple[int, int]]) -> List[Tuple[int, int]]:
    merged: List[List[int]] = []
    for start, end in sorted(intervals):
        if end <= start:
            continue
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return [(start, end) for start, end in merged]


def _covered_bp(intervals: Iterable[Tuple[int, int]]) -> int:
    return sum(end - start for start, end in _merge_intervals(intervals))


def _simple_overlap(cigar: Optional[str]) -> Optional[int]:
    """Return exact path overlap for match-only CIGARs, else None.

    Indel-bearing overlaps require a base-level coordinate transform.  Rather
    than silently making an approximate projection, callers split the context
    at such an edge.
    """
    if not cigar or cigar == "*":
        return None
    matches = list(_CIGAR_RE.finditer(cigar))
    if not matches or "".join(match.group(0) for match in matches) != cigar:
        return None
    if any(match.group(2) not in {"M", "=", "X"} for match in matches):
        return None
    return sum(int(match.group(1)) for match in matches)


class GfaAnnotator:
    """Build path contexts, run RepeatMasker, and summarize segment annotation."""

    def __init__(
        self,
        gfa_file: str,
        output_dir: str,
        te_library: str,
        threads: int = 4,
        repeatmasker_exe: str = "RepeatMasker",
        gfa_gaf: Optional[str] = None,
        gfa_gaf_fasta: Optional[str] = None,
        gfa_gaf_min_mapq: int = 20,
        gfa_gaf_min_identity: float = 0.90,
    ):
        self.gfa_file = os.path.abspath(gfa_file)
        self.output_dir = os.path.abspath(output_dir)
        self.te_library = os.path.abspath(te_library)
        self.threads = max(1, int(threads))
        self.repeatmasker_exe = repeatmasker_exe
        self.gfa_gaf = os.path.abspath(gfa_gaf) if gfa_gaf else None
        self.gfa_gaf_fasta = (os.path.abspath(gfa_gaf_fasta)
                              if gfa_gaf_fasta else None)
        self.gfa_gaf_min_mapq = int(gfa_gaf_min_mapq)
        self.gfa_gaf_min_identity = float(gfa_gaf_min_identity)
        if self.gfa_gaf_min_mapq < 0:
            raise ValueError('gfa_gaf_min_mapq must be non-negative')
        if not 0 <= self.gfa_gaf_min_identity <= 1:
            raise ValueError('gfa_gaf_min_identity must be between 0 and 1')
        self.te_families: Set[str] = set()

        self.nodes: Dict[str, Node] = {}
        self.traversals: List[Traversal] = []
        self.samples: Set[str] = set()
        self.presence: Dict[str, Set[str]] = defaultdict(set)
        self.presence_states: Dict[Tuple[str, str], str] = {}
        self._presence_is_tristate = False
        self.rgfa_gaf_presence_stats: Dict[str, int] = {}
        self.occurrences: List[Occurrence] = []
        self.query_sequences: Dict[str, str] = {}
        self._links: Dict[Tuple[str, str, str, str], str] = {}
        self._query_counter = 0
        self._unsupported_edges = 0

        os.makedirs(self.output_dir, exist_ok=True)

    @property
    def context_fasta(self) -> str:
        return os.path.join(self.output_dir, "gfa_annotation_context.fa")

    @property
    def manifest_path(self) -> str:
        return os.path.join(self.output_dir, "gfa_path_node_map.tsv")

    @property
    def node_table_path(self) -> str:
        return os.path.join(self.output_dir, "gfa_node_annotations.tsv")

    @property
    def rgfa_node_presence_path(self) -> str:
        return os.path.join(self.output_dir, "gfa_rgfa_node_presence.tsv")

    @property
    def hit_table_path(self) -> str:
        return os.path.join(self.output_dir, "gfa_node_te_hits.tsv")

    def annotate(self) -> str:
        if not os.path.isfile(self.te_library) or os.path.getsize(self.te_library) == 0:
            raise FileNotFoundError(f"TE library is missing or empty: {self.te_library}")
        with open(self.te_library, "rt") as library:
            self.te_families = {
                line[1:].strip().split()[0].partition("#")[0]
                for line in library if line.startswith(">")
            }

        self.parse_gfa()
        self.build_contexts()
        self.write_context_inputs()

        if self.query_sequences:
            repeat_out = self.run_repeatmasker()
            hits = self.parse_repeatmasker_output(repeat_out)
        else:
            logger.warning("GFA contains no available segment sequences; RepeatMasker skipped")
            hits = []

        projected = self.project_hits(hits)
        self.write_outputs(projected)
        logger.info(
            "GFA annotation complete: %d nodes, %d samples, %d RepeatMasker hits",
            len(self.nodes), len(self.samples), len(hits),
        )
        return self.node_table_path

    def parse_gfa(self) -> None:
        """Parse S/L/P/W records and collect exact sample traversals."""
        raw_paths: List[Tuple[str, List[Tuple[str, str]], List[Optional[str]]]] = []
        raw_walks: List[Traversal] = []

        with open(self.gfa_file, "rt") as handle:
            for line_number, raw_line in enumerate(handle, 1):
                if not raw_line.strip() or raw_line.startswith("#"):
                    continue
                fields = raw_line.rstrip("\n").split("\t")
                record_type = fields[0]

                if record_type == "S":
                    if len(fields) < 3:
                        raise ValueError(f"Malformed S line at {line_number}")
                    node_id, raw_sequence = fields[1], fields[2]
                    if node_id in self.nodes:
                        raise ValueError(f"Duplicate GFA segment ID: {node_id}")
                    tags = {tag[:2]: tag[5:] for tag in fields[3:]
                            if len(tag) >= 6 and tag[2:5] in {":i:", ":Z:"}}
                    sequence = None if raw_sequence == "*" else raw_sequence
                    length = len(sequence) if sequence is not None else int(tags.get("LN", 0))
                    self.nodes[node_id] = Node(node_id, sequence, length)

                elif record_type == "L":
                    if len(fields) < 6:
                        raise ValueError(f"Malformed L line at {line_number}")
                    key = (fields[1], fields[2], fields[3], fields[4])
                    self._links[key] = fields[5]
                    reverse_key = (
                        fields[3], "+" if fields[4] == "-" else "-",
                        fields[1], "+" if fields[2] == "-" else "-",
                    )
                    self._links[reverse_key] = fields[5]

                elif record_type == "P":
                    if len(fields) < 3:
                        raise ValueError(f"Malformed P line at {line_number}")
                    path_name = fields[1]
                    step_tokens = re.split(r"[,;]", fields[2])
                    steps = []
                    for token in step_tokens:
                        if len(token) < 2 or token[-1] not in "+-":
                            raise ValueError(
                                f"Malformed P step {token!r} at line {line_number}")
                        steps.append((token[:-1], token[-1]))
                    overlap_field = fields[3] if len(fields) > 3 else "*"
                    if overlap_field == "*":
                        overlaps: List[Optional[str]] = [None] * max(0, len(steps) - 1)
                    else:
                        overlaps = overlap_field.split(",")
                        if len(overlaps) != max(0, len(steps) - 1):
                            raise ValueError(
                                f"P overlap count does not match steps at line {line_number}")
                    raw_paths.append((path_name, steps, overlaps))

                elif record_type == "W":
                    if len(fields) < 7:
                        raise ValueError(f"Malformed W line at {line_number}")
                    sample, haplotype, seq_id = fields[1], fields[2], fields[3]
                    parsed = _WALK_STEP_RE.findall(fields[6])
                    if not parsed or "".join(direction + node for direction, node in parsed) != fields[6]:
                        raise ValueError(f"Malformed W walk at line {line_number}")
                    steps = tuple((node, "+" if direction == ">" else "-")
                                  for direction, node in parsed)
                    traversal_id = (
                        f"W:{sample}#{haplotype}#{seq_id}:{fields[4]}-{fields[5]}"
                    )
                    raw_walks.append(Traversal(
                        traversal_id, sample, steps,
                        tuple("0M" for _ in range(max(0, len(steps) - 1))), "W",
                    ))

        if not self.nodes:
            raise ValueError("GFA contains no S records")

        for path_name, steps, overlaps in raw_paths:
            resolved = []
            for index, overlap in enumerate(overlaps):
                if overlap is None:
                    left_id, left_orient = steps[index]
                    right_id, right_orient = steps[index + 1]
                    overlap = self._links.get(
                        (left_id, left_orient, right_id, right_orient))
                resolved.append(overlap)
            sample = _sample_from_path_name(path_name)
            self.traversals.append(Traversal(
                f"P:{path_name}", sample, tuple(steps), tuple(resolved), "P"))

        self.traversals.extend(raw_walks)
        for traversal in self.traversals:
            self.samples.add(traversal.sample)
            for node_id, _ in traversal.steps:
                if node_id not in self.nodes:
                    raise ValueError(
                        f"Traversal {traversal.traversal_id} references missing segment {node_id}")
                self.presence[node_id].add(traversal.sample)

        if not self.traversals:
            if self.gfa_gaf:
                samples, states, presence_stats = derive_rgfa_node_presence(
                    self.gfa_file, self.gfa_gaf,
                    min_mapq=self.gfa_gaf_min_mapq,
                    min_identity=self.gfa_gaf_min_identity,
                    query_fasta=self.gfa_gaf_fasta)
                self.samples.update(samples)
                self.presence_states.update(states)
                self._presence_is_tristate = True
                self.rgfa_gaf_presence_stats = presence_stats
                logger.info(
                    "Recovered conservative rGFA node presence from GAF: %d samples, %d complete traversals",
                    len(samples), len(states))
            else:
                logger.warning(
                    "GFA has no P/W records: nodes can be annotated, but sample presence "
                    "columns cannot be inferred; provide --gfa-gaf for rGFA")

    def _new_query_id(self, prefix: str) -> str:
        self._query_counter += 1
        return f"{prefix}_{self._query_counter:08d}"

    def build_contexts(self) -> None:
        represented: Set[str] = set()

        for traversal in self.traversals:
            chunk_sequence = ""
            pending: List[Tuple[str, int, int, int, int, str]] = []

            def flush_chunk() -> None:
                nonlocal chunk_sequence, pending
                if not chunk_sequence:
                    pending = []
                    return
                targets = [item for item in pending if item[0] not in represented]
                if not targets:
                    chunk_sequence = ""
                    pending = []
                    return

                # The first traversal contributes its full sequence. Later
                # traversals contribute only windows around previously unseen
                # nodes, avoiding N-sample copies of the shared backbone while
                # retaining 10 kb of RepeatMasker boundary context.
                windows = _merge_intervals(
                    (max(0, query_start - 10000),
                     min(len(chunk_sequence), query_end + 10000))
                    for _, query_start, query_end, _, _, _ in targets)
                for window_start, window_end in windows:
                    query_id = self._new_query_id("path")
                    self.query_sequences[query_id] = chunk_sequence[window_start:window_end]
                    for (node_id, query_start, query_end,
                         node_start, node_end, orient) in pending:
                        if query_end <= window_start or query_start >= window_end:
                            continue
                        clipped_start = max(query_start, window_start)
                        clipped_end = min(query_end, window_end)
                        left_offset = clipped_start - query_start
                        right_offset = clipped_end - query_start
                        if orient == "+":
                            projected_node_start = node_start + left_offset
                            projected_node_end = node_start + right_offset
                        else:
                            projected_node_start = node_end - right_offset
                            projected_node_end = node_end - left_offset
                        self.occurrences.append(Occurrence(
                            query_id, traversal.traversal_id, traversal.sample, node_id,
                            clipped_start - window_start, clipped_end - window_start,
                            projected_node_start, projected_node_end, orient,
                        ))
                        if clipped_start <= query_start and clipped_end >= query_end:
                            represented.add(node_id)
                chunk_sequence = ""
                pending = []

            for index, (node_id, orientation) in enumerate(traversal.steps):
                node = self.nodes[node_id]
                if node.sequence is None:
                    flush_chunk()
                    continue
                oriented_sequence = (
                    node.sequence if orientation == "+" else reverse_complement(node.sequence)
                )

                if not chunk_sequence:
                    chunk_sequence = oriented_sequence
                    pending.append((node_id, 0, len(oriented_sequence), 0,
                                    len(oriented_sequence), orientation))
                    continue

                cigar = traversal.overlaps[index - 1]
                overlap = _simple_overlap(cigar)
                if (overlap is None
                        or overlap > min(len(chunk_sequence), len(oriented_sequence))
                        or (overlap and chunk_sequence[-overlap:].upper()
                            != oriented_sequence[:overlap].upper())):
                    self._unsupported_edges += 1
                    flush_chunk()
                    chunk_sequence = oriented_sequence
                    pending.append((node_id, 0, len(oriented_sequence), 0,
                                    len(oriented_sequence), orientation))
                    continue

                query_start = len(chunk_sequence) - overlap
                query_end = query_start + len(oriented_sequence)
                chunk_sequence += oriented_sequence[overlap:]
                pending.append((node_id, query_start, query_end, 0,
                                len(oriented_sequence), orientation))

            flush_chunk()

        # A graph segment may be absent from every emitted traversal, or a graph
        # may have no P/W records.  Annotate such sequences directly as a safe
        # fallback.  They lack boundary context but are never silently omitted.
        for node_id, node in self.nodes.items():
            if node_id in represented or node.sequence is None:
                continue
            query_id = self._new_query_id("node")
            self.query_sequences[query_id] = node.sequence
            self.occurrences.append(Occurrence(
                query_id, "unplaced", "", node_id, 0, len(node.sequence),
                0, len(node.sequence), "+",
            ))

        if self._unsupported_edges:
            logger.warning(
                "%d path edges had missing/complex overlap CIGARs; contexts were "
                "split at those edges to preserve exact coordinates",
                self._unsupported_edges,
            )

    def write_context_inputs(self) -> None:
        with open(self.context_fasta, "wt") as fasta:
            for query_id, sequence in self.query_sequences.items():
                fasta.write(f">{query_id}\n")
                for start in range(0, len(sequence), 80):
                    fasta.write(sequence[start:start + 80] + "\n")

        with open(self.manifest_path, "wt") as manifest:
            manifest.write(
                "query_id\tquery_start\tquery_end\tnode_id\tnode_start\t"
                "node_end\torientation\tpath_id\tsample\n")
            for occurrence in self.occurrences:
                manifest.write(
                    f"{occurrence.query_id}\t{occurrence.query_start}\t"
                    f"{occurrence.query_end}\t{occurrence.node_id}\t"
                    f"{occurrence.node_start}\t{occurrence.node_end}\t"
                    f"{occurrence.orientation}\t{occurrence.traversal_id}\t"
                    f"{occurrence.sample}\n")

    def run_repeatmasker(self) -> str:
        executable = shutil.which(self.repeatmasker_exe)
        if executable is None and not os.path.isfile(self.repeatmasker_exe):
            raise FileNotFoundError(
                f"RepeatMasker executable not found: {self.repeatmasker_exe}")
        executable = executable or self.repeatmasker_exe

        repeat_dir = os.path.join(self.output_dir, "repeatmasker")
        os.makedirs(repeat_dir, exist_ok=True)
        command = [
            executable,
            "-pa", str(self.threads),
            "-lib", self.te_library,
            "-dir", repeat_dir,
            "-s", "-gff",
            self.context_fasta,
        ]
        logger.info("Running RepeatMasker on %d GFA sequence contexts", len(self.query_sequences))
        result = subprocess.run(command, capture_output=True, text=True, check=False)
        if result.returncode != 0:
            tail = "\n".join((result.stderr or "").strip().splitlines()[-20:])
            raise RuntimeError(
                f"RepeatMasker failed with exit code {result.returncode}:\n{tail}")

        output = os.path.join(repeat_dir, os.path.basename(self.context_fasta) + ".out")
        if not os.path.isfile(output):
            raise RuntimeError(f"RepeatMasker produced no annotation file: {output}")
        return output

    def parse_repeatmasker_output(self, output_file: str) -> List[RepeatHit]:
        hits: List[RepeatHit] = []
        skipped_nonlibrary = 0
        with open(output_file, "rt") as handle:
            for line_number, line in enumerate(handle, 1):
                fields = line.split()
                if len(fields) < 14:
                    continue
                try:
                    score = int(fields[0])
                    divergence = float(fields[1])
                    deletion = float(fields[2])
                    insertion = float(fields[3])
                    query_start = int(fields[5]) - 1
                    query_end = int(fields[6])
                except ValueError:
                    continue
                if query_end <= query_start:
                    logger.warning("Ignoring invalid RepeatMasker interval at line %d", line_number)
                    continue
                if fields[9] not in self.te_families:
                    skipped_nonlibrary += 1
                    continue
                hits.append(RepeatHit(
                    score, divergence, deletion, insertion,
                    fields[4], query_start, query_end, fields[8], fields[9], fields[10],
                ))
        if skipped_nonlibrary:
            logger.info("Ignored %d RepeatMasker hits absent from the custom TE library",
                        skipped_nonlibrary)
        return hits

    @staticmethod
    def _project_interval(hit: RepeatHit, occurrence: Occurrence) -> Tuple[int, int]:
        overlap_start = max(hit.query_start, occurrence.query_start)
        overlap_end = min(hit.query_end, occurrence.query_end)
        left_offset = overlap_start - occurrence.query_start
        right_offset = overlap_end - occurrence.query_start
        if occurrence.orientation == "+":
            return occurrence.node_start + left_offset, occurrence.node_start + right_offset
        return occurrence.node_end - right_offset, occurrence.node_end - left_offset

    def project_hits(self, hits: Sequence[RepeatHit]) -> List[dict]:
        occurrences_by_query: Dict[str, List[Occurrence]] = defaultdict(list)
        hits_by_query: Dict[str, List[RepeatHit]] = defaultdict(list)
        for occurrence in self.occurrences:
            occurrences_by_query[occurrence.query_id].append(occurrence)
        for hit in hits:
            if hit.query_id not in self.query_sequences:
                logger.warning("RepeatMasker hit references unknown query %s", hit.query_id)
                continue
            hits_by_query[hit.query_id].append(hit)

        projected: List[dict] = []
        for query_id, query_hits in hits_by_query.items():
            occurrences = sorted(
                occurrences_by_query[query_id], key=lambda item: item.query_start)
            starts = [occurrence.query_start for occurrence in occurrences]
            left = 0
            for hit in sorted(query_hits, key=lambda item: item.query_start):
                while left < len(occurrences) and occurrences[left].query_end <= hit.query_start:
                    left += 1
                stop = bisect.bisect_left(starts, hit.query_end)
                for occurrence in occurrences[left:stop]:
                    node_start, node_end = self._project_interval(hit, occurrence)
                    node = self.nodes[occurrence.node_id]
                    node_start = max(0, min(node.length, node_start))
                    node_end = max(0, min(node.length, node_end))
                    if node_end <= node_start:
                        continue
                    projected.append({
                        "node_id": occurrence.node_id,
                        "node_start": node_start,
                        "node_end": node_end,
                        "te_family": hit.repeat_name,
                        "te_class": hit.repeat_class,
                        "query_id": hit.query_id,
                        "query_start": hit.query_start,
                        "query_end": hit.query_end,
                        "path_id": occurrence.traversal_id,
                        "sample": occurrence.sample,
                        "score": hit.score,
                        "divergence": hit.divergence,
                        "deletion": hit.deletion,
                        "insertion": hit.insertion,
                        "strand": hit.strand,
                    })
        return projected

    def write_outputs(self, projected: Sequence[dict]) -> None:
        fields = [
            "node_id", "node_start", "node_end", "overlap_bp", "te_family",
            "te_class", "query_id", "query_start", "query_end", "path_id",
            "sample", "score", "divergence", "deletion", "insertion", "strand",
        ]
        with open(self.hit_table_path, "wt") as output:
            output.write("\t".join(fields) + "\n")
            for row in sorted(projected, key=lambda item: (
                    item["node_id"], item["node_start"], item["node_end"],
                    item["te_family"], item["path_id"])):
                values = dict(row)
                values["overlap_bp"] = row["node_end"] - row["node_start"]
                output.write("\t".join(str(values[field]) for field in fields) + "\n")

        by_node_label: Dict[str, Dict[Tuple[str, str], List[Tuple[int, int]]]] = defaultdict(
            lambda: defaultdict(list))
        hit_counts: Dict[str, int] = defaultdict(int)
        for row in projected:
            label = (row["te_family"], row["te_class"])
            by_node_label[row["node_id"]][label].append(
                (row["node_start"], row["node_end"]))
            hit_counts[row["node_id"]] += 1

        fixed_columns = [
            "node_id", "length", "annotation_status", "te_annotations",
            "primary_te_family", "primary_te_class", "te_covered_bp",
            "te_fraction", "te_uncovered_bp", "te_type_count",
            "projected_hit_count",
        ]
        sample_columns = [f"sample:{sample}" for sample in sorted(self.samples)]
        with open(self.node_table_path, "wt") as output:
            output.write("\t".join(fixed_columns + sample_columns) + "\n")
            for node_id, node in self.nodes.items():
                label_intervals = by_node_label.get(node_id, {})
                label_coverage = {
                    label: _covered_bp(intervals)
                    for label, intervals in label_intervals.items()
                }
                ordered_labels = sorted(
                    label_coverage, key=lambda label: (-label_coverage[label], label))
                all_intervals = [interval
                                 for intervals in label_intervals.values()
                                 for interval in intervals]
                covered = min(node.length, _covered_bp(all_intervals))

                if node.sequence is None:
                    status = "SequenceUnavailable"
                elif covered == 0:
                    status = "None"
                elif covered < node.length:
                    status = "PartialTE"
                else:
                    status = "TE"

                if ordered_labels:
                    annotation = ";".join(
                        f"{family}#{te_class}" for family, te_class in ordered_labels)
                    primary_family, primary_class = ordered_labels[0]
                else:
                    annotation = "None"
                    primary_family = "None"
                    primary_class = "None"

                fraction = covered / node.length if node.length else 0.0
                row = [
                    node_id, str(node.length), status, annotation,
                    primary_family, primary_class, str(covered), f"{fraction:.6f}",
                    str(max(0, node.length - covered)), str(len(ordered_labels)),
                    str(hit_counts.get(node_id, 0)),
                ]
                if self._presence_is_tristate:
                    # A missing rGFA GAF traversal is not evidence of absence.
                    row.extend(self.presence_states.get((node_id, sample), "NA")
                               for sample in sorted(self.samples))
                else:
                    present_samples = self.presence.get(node_id, set())
                    row.extend("1" if sample in present_samples else "0"
                               for sample in sorted(self.samples))
                output.write("\t".join(row) + "\n")

        if self._presence_is_tristate:
            with open(self.rgfa_node_presence_path, "wt") as output:
                output.write("node_id\tsample\tstate\n")
                for (node_id, sample), state in sorted(self.presence_states.items()):
                    output.write(f"{node_id}\t{sample}\t{state}\n")


__all__ = [
    "GfaAnnotator", "Node", "Occurrence", "RepeatHit", "Traversal",
    "_covered_bp", "_merge_intervals", "_simple_overlap",
]
