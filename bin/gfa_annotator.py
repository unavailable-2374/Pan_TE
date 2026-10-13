#!/usr/bin/env python3
"""Project per-assembly GFF3/TE annotations onto GFA node occurrences.

The occurrence-level table is authoritative.  A node-level identity is a
summary because the same graph segment can occur in different genomic contexts.
All internal coordinates are 0-based, half-open; GFF3 is converted on input and
converted back to 1-based inclusive only in the node-local GFF3 output.
"""

from __future__ import annotations

import argparse
import bz2
import csv
import gzip
import hashlib
import json
import logging
import os
import re
from array import array
from bisect import bisect_right
from collections import Counter, defaultdict
from dataclasses import dataclass
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Set, Tuple
from urllib.parse import quote, unquote

import numpy as np

from rgfa_gaf import RgfaGafError, iter_gaf_records, load_rgfa_graph


logger = logging.getLogger(__name__)
_WALK_STEP_RE = re.compile(r"([><])([^><]+)")
_CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
_TE_WORD_RE = re.compile(r"(?:transpos|retrotranspos|retroelement|repeat_region|mobile_element)", re.I)
_NON_TE_RM_CLASS_RE = re.compile(
    r"^(?:Simple_repeat|Low_complexity|Satellite(?:/.*)?|rRNA|tRNA|snRNA|scRNA|srpRNA)$",
    re.I,
)
_GENE_TYPES = {
    "gene", "pseudogene", "mrna", "transcript", "exon", "cds", "intron",
    "five_prime_utr", "three_prime_utr", "utr", "ncrna", "trna", "rrna",
    "snrna", "snorna", "mirna", "lncrna",
}
# PGGB/smoothxg nodes are often only tens of bases long.  A 1 kb feature bin
# keeps node queries selective without building a base-level interval index.
_BIN_SIZE = 1_000
# Shared by every segment without SAM-style tags (all PGGB segments); never mutated.
_NO_TAGS: Dict[str, str] = {}
# Rows handed from numpy columns to Python per chunk while streaming outputs.
_CHUNK = 1 << 20

# RepeatMasker-compatible report rows. Aggregates and subtypes are both
# reported; occupied lengths are independent interval unions and must not be
# summed across rows.
_RM_ROWS = (
    ("retroelements", "Retroelements"),
    ("sines", "    SINEs:"),
    ("penelope", "    Penelope:"),
    ("lines", "    LINEs:"),
    ("line_cre", "        CRE/SLACS:"),
    ("line_l2", "        L2/CR1/Rex:"),
    ("line_r1", "        R1/LOA/Jockey:"),
    ("line_r2", "        R2/R4/NeSL:"),
    ("line_rte", "        RTE/Bov-B:"),
    ("line_l1", "        L1/CIN4:"),
    ("ltr", "    LTR elements:"),
    ("ltr_bel", "        BEL/Pao:"),
    ("ltr_copia", "        Ty1/Copia:"),
    ("ltr_gypsy", "        Gypsy/DIRS1:"),
    ("retroviral", "        Retroviral:"),
    ("dna", "DNA transposons"),
    ("dna_hobo", "    hobo-Activator:"),
    ("dna_tc1", "    Tc1-IS630-Pogo:"),
    ("dna_ensp", "    En-Spm:"),
    ("dna_mule", "    MuDR-IS905:"),
    ("dna_piggybac", "    PiggyBac:"),
    ("dna_tourist", "    Tourist/Harbinger:"),
    ("dna_other", "    Other (Mirage, P-element, Transib):"),
    ("rolling", "Rolling-circles"),
    ("unclassified", "Unclassified:"),
    ("interspersed", "Total interspersed repeats:"),
    ("small_rna", "Small RNA:"),
    ("satellite", "Satellites:"),
    ("simple_repeat", "Simple repeats:"),
    ("low_complexity", "Low complexity:"),
)
def _open_text(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    if path.endswith(".bz2"):
        return bz2.open(path, "rt")
    return open(path, "rt")


def _sha256(path: str) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _as_bool(value: str) -> bool:
    return str(value or "").strip().lower() in {"1", "true", "yes", "y", "complete"}


def _flip(orientation: str) -> str:
    return "-" if orientation == "+" else "+"


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


def _overlap_bp(left: Iterable[Tuple[int, int]], right: Iterable[Tuple[int, int]]) -> int:
    a, b = _merge_intervals(left), _merge_intervals(right)
    i = j = total = 0
    while i < len(a) and j < len(b):
        total += max(0, min(a[i][1], b[j][1]) - max(a[i][0], b[j][0]))
        if a[i][1] <= b[j][1]:
            i += 1
        else:
            j += 1
    return total


def _simple_overlap(cigar: Optional[str]) -> Optional[int]:
    """Return sequence overlap for match-only CIGARs, otherwise None."""
    if not cigar or cigar == "*":
        return None
    matches = list(_CIGAR_RE.finditer(cigar))
    if not matches or "".join(match.group(0) for match in matches) != cigar:
        return None
    if any(match.group(2) not in {"M", "=", "X"} for match in matches):
        return None
    return sum(int(match.group(1)) for match in matches)


def _gapless_cigar(cigar: str, expected: int) -> bool:
    if not cigar:
        return False
    matches = list(_CIGAR_RE.finditer(cigar))
    return bool(matches) and "".join(x.group(0) for x in matches) == cigar and all(
        x.group(2) in {"M", "=", "X"} for x in matches
    ) and sum(int(x.group(1)) for x in matches) == expected


def _attributes(text: str) -> Dict[str, str]:
    result = {}
    for token in text.strip().split(";"):
        if not token:
            continue
        if "=" in token:
            key, value = token.split("=", 1)
        elif " " in token:
            key, value = token.split(" ", 1)
        else:
            continue
        result[unquote(key.strip())] = unquote(value.strip().strip('"'))
    return result


def _clean(value: object) -> str:
    return str(value if value is not None else "").replace("\t", " ").replace("\n", " ")


def _repeat_summary_buckets(category: str, feature_type: str, family: str,
                            te_class: str) -> Tuple[str, ...]:
    """Map one annotation record to RepeatMasker-style aggregate/subtype rows."""
    # An explicit class is authoritative. Family names are intentionally not
    # used to override it: RepeatMasker libraries can contain LINE/LTR-like
    # strings in names that are formally reported as Unknown/SubclassI.
    text = " ".join((feature_type, te_class or family)).lower()
    normalized = re.sub(r"[^a-z0-9]+", " ", text)
    words = set(normalized.split())
    if category != "TE":
        if any(token in words for token in ("rrna", "trna", "snrna", "scrna", "srprna")):
            return ("small_rna", "masked")
        if "satellite" in words:
            return ("satellite", "masked")
        if "simple" in words and "repeat" in words:
            return ("simple_repeat", "masked")
        if "low" in words and "complexity" in words:
            return ("low_complexity", "masked")
        return ()

    buckets = {"masked"}
    raw_class = (te_class or family).strip().lower()
    line_class = raw_class == "line" or raw_class.startswith("line/")
    ltr_class = raw_class == "ltr" or raw_class.startswith("ltr/")
    dna_class = raw_class == "dna" or raw_class.startswith("dna/")
    if ("helitron" in words or "rolling" in words or
            re.search(r"(?:^|[/_ -])rc(?:$|[/_ -])", text)):
        buckets.add("rolling")
    elif raw_class == "sine" or raw_class.startswith("sine/"):
        buckets.update(("retroelements", "sines"))
    elif "penelope" in words:
        buckets.update(("retroelements", "penelope"))
    elif line_class:
        buckets.update(("retroelements", "lines"))
        if words.intersection({"cre", "slacs"}):
            buckets.add("line_cre")
        elif words.intersection({"l2", "cr1", "rex"}):
            buckets.add("line_l2")
        elif words.intersection({"r1", "loa", "jockey"}):
            buckets.add("line_r1")
        elif words.intersection({"r2", "r4", "nesl"}):
            buckets.add("line_r2")
        elif words.intersection({"rte", "bov", "bovb"}):
            buckets.add("line_rte")
        elif words.intersection({"l1", "cin4"}):
            buckets.add("line_l1")
    elif ltr_class:
        buckets.update(("retroelements", "ltr"))
        if words.intersection({"bel", "pao"}):
            buckets.add("ltr_bel")
        elif "copia" in words:
            buckets.add("ltr_copia")
        elif words.intersection({"gypsy", "dirs", "dirs1"}):
            buckets.add("ltr_gypsy")
        elif "retroviral" in words or any(word.startswith("erv") for word in words):
            buckets.add("retroviral")
    elif raw_class in {"bel-pao", "bel/pao"}:
        buckets.add("ltr_bel")
    elif raw_class in {"dirs", "dirs1"}:
        buckets.add("ltr_gypsy")
    elif raw_class in {"pif-harbinger", "pif/harbinger"}:
        buckets.add("dna_tourist")
    elif dna_class:
        buckets.add("dna")
        if words.intersection({"hobo", "activator", "hat"}):
            buckets.add("dna_hobo")
        elif words.intersection({"tc1", "tcmar", "is630", "pogo", "mariner"}):
            buckets.add("dna_tc1")
        elif raw_class.startswith(("dna/ensp", "dna/en-spm")):
            buckets.add("dna_ensp")
        elif words.intersection({"mule", "mudr", "is905"}):
            buckets.add("dna_mule")
        elif "piggybac" in words:
            buckets.add("dna_piggybac")
        elif words.intersection({"tourist", "harbinger", "pif"}):
            buckets.add("dna_tourist")
        elif words.intersection({"mirage", "transib"}) or "p-element" in raw_class:
            buckets.add("dna_other")
    else:
        buckets.add("unclassified")
    if buckets.intersection({"retroelements", "dna", "unclassified"}):
        buckets.add("interspersed")
    ordered = tuple(key for key, _label in _RM_ROWS if key in buckets)
    return ordered + (("masked",) if "masked" in buckets else ())


class _StreamingUnion:
    """Accumulate union lengths for independently sorted interval streams."""

    def __init__(self):
        self.last: Dict[Tuple, Tuple[int, int]] = {}
        self.total: Counter = Counter()

    def add(self, key: Tuple, start: int, end: int) -> None:
        if end <= start:
            return
        previous = self.last.get(key)
        if previous is None or start > previous[1]:
            self.total[key] += end - start
            self.last[key] = (start, end)
        elif end > previous[1]:
            self.total[key] += end - previous[1]
            self.last[key] = (previous[0], end)


def _fasta_sequence_stats(path: str) -> Tuple[int, int, int, Optional[float]]:
    """Return sequence count, total bp, N/X bp, and GC percentage."""
    sequences = total = nx = gc = 0
    with _open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                sequences += 1
                continue
            sequence = "".join(line.split()).upper()
            total += len(sequence)
            nx += sequence.count("N") + sequence.count("X")
            gc += sequence.count("G") + sequence.count("C")
    denominator = total - nx
    return sequences, total, nx, (100.0 * gc / denominator if denominator else None)


@dataclass(frozen=True, slots=True)
class Node:
    node_id: str
    sequence: Optional[str]
    length: int
    tags: Dict[str, str]
    # Position in GFA segment order; column stores reference nodes by it.
    index: int = -1


@dataclass(frozen=True)
class ManifestEntry:
    sample: str
    haplotype: str
    fasta: Optional[str]
    gff3: Optional[str]
    te_annotation: Optional[str]
    seqid_map: Optional[str]
    path_prefix: str
    gaf: Optional[str]
    annotation_complete: bool

    @property
    def key(self) -> Tuple[str, str]:
        return self.sample, self.haplotype


@dataclass(frozen=True, slots=True)
class Feature:
    entry_key: Tuple[str, str]
    seqid: str
    start: int
    end: int
    strand: str
    phase: str
    feature_type: str
    category: str
    feature_id: str
    parent: str
    family: str
    te_class: str
    source: str
    attributes: str


@dataclass(frozen=True, slots=True)
class Occurrence:
    occurrence_id: str
    entry_key: Tuple[str, str]
    sample: str
    haplotype: str
    seqid: str
    node_id: str
    occurrence_index: int
    genome_start: int
    genome_end: int
    node_start: int
    node_end: int
    orientation: str
    path_id: str
    evidence: str
    coordinate_status: str


@dataclass(frozen=True, slots=True)
class ProjectedFeature:
    occurrence_id: str
    node_id: str
    node_start: int
    node_end: int
    sample: str
    haplotype: str
    seqid: str
    genome_start: int
    genome_end: int
    orientation: str
    feature_type: str
    category: str
    feature_id: str
    parent: str
    family: str
    te_class: str
    source: str
    strand: str
    phase: str
    attributes: str


class _Interner:
    """Map repeated values (sources, seqids, path ids) to small integer codes."""

    __slots__ = ("values", "codes")

    def __init__(self):
        self.values: List = []
        self.codes: Dict = {}

    def code(self, value) -> int:
        code = self.codes.get(value)
        if code is None:
            code = self.codes[value] = len(self.values)
            self.values.append(value)
        return code


class _OccurrenceStore:
    """Column store for occurrences in insertion order.

    A pangenome graph yields one occurrence per path step (hundreds of millions
    for PGGB graphs), so rows are kept as packed columns and materialized as
    ``Occurrence`` only when read.  Row ``i`` has ``occurrence_id`` ``occ_{i+1}``.
    """

    def __init__(self, annotator: "GfaAnnotator"):
        self._annotator = annotator
        self.sources = _Interner()      # (entry_key, seqid)
        self.paths = _Interner()
        self.labels = _Interner()       # evidence and coordinate_status strings
        self.source = array("I")
        self.node = array("I")
        self.index = array("q")
        self.genome_start = array("q")
        self.genome_end = array("q")
        self.node_start = array("q")
        self.node_end = array("q")
        self.reverse = array("B")
        self.path = array("I")
        self.evidence = array("H")
        self.status = array("H")

    def append(self, entry_key: Tuple[str, str], seqid: str, node: int, index: int,
               genome_start: int, genome_end: int, node_start: int, node_end: int,
               orientation: str, path_id: str, evidence: str,
               coordinate_status: str) -> None:
        self.source.append(self.sources.code((entry_key, seqid)))
        self.node.append(node)
        self.index.append(index)
        self.genome_start.append(genome_start)
        self.genome_end.append(genome_end)
        self.node_start.append(node_start)
        self.node_end.append(node_end)
        self.reverse.append(orientation == "-")
        self.path.append(self.paths.code(path_id))
        self.evidence.append(self.labels.code(evidence))
        self.status.append(self.labels.code(coordinate_status))

    def key(self, i: int) -> Tuple:
        entry_key, seqid = self.sources.values[self.source[i]]
        return (entry_key, seqid, self._annotator.node_ids[self.node[i]],
                self.genome_start[i], self.genome_end[i], self.node_start[i],
                self.node_end[i], "-" if self.reverse[i] else "+")

    def __len__(self) -> int:
        return len(self.source)

    def __getitem__(self, i: int) -> Occurrence:
        if i < 0:
            i += len(self)
        entry_key, seqid = self.sources.values[self.source[i]]
        entry = self._annotator.entries[entry_key]
        return Occurrence(
            f"occ_{i + 1:012d}", entry_key, entry.sample, entry.haplotype, seqid,
            self._annotator.node_ids[self.node[i]], self.index[i],
            self.genome_start[i], self.genome_end[i], self.node_start[i],
            self.node_end[i], "-" if self.reverse[i] else "+",
            self.paths.values[self.path[i]], self.labels.values[self.evidence[i]],
            self.labels.values[self.status[i]],
        )

    def __iter__(self) -> Iterator[Occurrence]:
        return (self[i] for i in range(len(self)))

    def __eq__(self, other) -> bool:
        if isinstance(other, (list, tuple)):
            return list(self) == list(other)
        return NotImplemented

    __hash__ = None


class _ProjectionStore:
    """Column store for feature projections, kept in occurrence order."""

    def __init__(self, annotator: "GfaAnnotator"):
        self._annotator = annotator
        self.occurrence = array("I")
        self.feature = array("I")
        self.node_start = array("q")
        self.node_end = array("q")
        self.genome_start = array("q")
        self.genome_end = array("q")

    def __len__(self) -> int:
        return len(self.occurrence)

    def __getitem__(self, i: int) -> ProjectedFeature:
        if i < 0:
            i += len(self)
        annotator = self._annotator
        occurrences = annotator.occurrences
        o = self.occurrence[i]
        feature = annotator.features[self.feature[i]]
        entry_key, seqid = occurrences.sources.values[occurrences.source[o]]
        entry = annotator.entries[entry_key]
        if occurrences.reverse[o]:
            orientation = "-"
            strand = _flip(feature.strand) if feature.strand in "+-" else feature.strand
        else:
            orientation = "+"
            strand = feature.strand
        return ProjectedFeature(
            f"occ_{o + 1:012d}", annotator.node_ids[occurrences.node[o]],
            self.node_start[i], self.node_end[i], entry.sample, entry.haplotype, seqid,
            self.genome_start[i], self.genome_end[i], orientation,
            feature.feature_type, feature.category, feature.feature_id,
            feature.parent, feature.family, feature.te_class, feature.source,
            strand, feature.phase, feature.attributes,
        )

    def __iter__(self) -> Iterator[ProjectedFeature]:
        return (self[i] for i in range(len(self)))


def _column(values: array, dtype) -> np.ndarray:
    """Zero-copy numpy view of an ``array.array`` column (empty-safe)."""
    if not len(values):
        return np.zeros(0, dtype=dtype)
    return np.frombuffer(values, dtype=dtype)


def _group_offsets(sorted_keys: np.ndarray, groups: int) -> np.ndarray:
    return np.searchsorted(sorted_keys, np.arange(groups + 1), side="left")


class GfaAnnotator:
    """Project all-genome annotations and summarize GFA node identities."""

    def __init__(
        self,
        gfa_file: str,
        output_dir: str,
        annotation_manifest: str,
        gfa_gaf: Optional[str] = None,
        gfa_gaf_min_mapq: int = 20,
        gfa_gaf_min_identity: float = 0.90,
        te_min_fraction: float = 0.80,
        partial_min_fraction: float = 0.10,
    ):
        self.gfa_file = os.path.abspath(gfa_file)
        self.output_dir = os.path.abspath(output_dir)
        self.annotation_manifest = os.path.abspath(annotation_manifest)
        self.default_gaf = os.path.abspath(gfa_gaf) if gfa_gaf else None
        self.min_mapq = int(gfa_gaf_min_mapq)
        self.min_identity = float(gfa_gaf_min_identity)
        self.te_min_fraction = float(te_min_fraction)
        self.partial_min_fraction = float(partial_min_fraction)
        if not 0 <= self.partial_min_fraction <= self.te_min_fraction <= 1:
            raise ValueError("annotation fractions must satisfy 0 <= partial <= TE <= 1")

        self.entries: Dict[Tuple[str, str], ManifestEntry] = {}
        self.seqid_maps: Dict[Tuple[str, str], Dict[str, str]] = defaultdict(dict)
        self.features: List[Feature] = []
        self.features_by_source: Dict[Tuple[Tuple[str, str], str], List[int]] = defaultdict(list)
        self.feature_bins: Dict[Tuple[Tuple[str, str], str, int], List[int]] = defaultdict(list)
        self.annotation_seqids: Dict[Tuple[str, str], Set[str]] = defaultdict(set)
        self.nodes: Dict[str, Node] = {}
        # Segment ids by Node.index, followed by ids that occur only in L lines
        # or alignments (indexed through extra_node_index).
        self.node_ids: List[str] = []
        self.extra_node_index: Dict[str, int] = {}
        # GFA links as sorted packed keys (see _link_key) with interned CIGARs.
        self.link_keys = np.zeros(0, dtype=np.int64)
        self.link_cigars = np.zeros(0, dtype=np.int32)
        self.link_cigar_values: List[str] = []
        self.occurrences = _OccurrenceStore(self)
        # Built only once a duplicate occurrence becomes possible (see _add_occurrence).
        self.occurrence_keys: Optional[Set[Tuple]] = None
        self.embedded_sources: Set[Tuple[Tuple[str, str], str]] = set()
        # Genome spans of embedded traversals per source, for duplicate-free appends.
        self._source_spans: Dict[Tuple[Tuple[str, str], str], List[Tuple[int, int]]] = {}
        self.projected = _ProjectionStore(self)
        self.unresolved: List[Tuple[str, str, str]] = []
        # Per Node.index: (identity, is_te, te_fraction, family, class, gene_structure, confidence).
        self.node_summaries: List[Optional[Tuple[str, str, str, str, str, str, str]]] = []
        self.gaf_stats: Counter = Counter()
        self.graph_mode = "unknown"
        os.makedirs(self.output_dir, exist_ok=True)

    @property
    def interval_table_path(self) -> str:
        return os.path.join(self.output_dir, "gfa_node_annotation_intervals.tsv")

    @property
    def occurrence_table_path(self) -> str:
        return os.path.join(self.output_dir, "gfa_node_occurrence_annotations.tsv")

    @property
    def node_table_path(self) -> str:
        return os.path.join(self.output_dir, "gfa_node_annotations.tsv")

    @property
    def matrix_path(self) -> str:
        return os.path.join(self.output_dir, "gfa_node_sample_matrix.tsv")

    @property
    def repeat_table_path(self) -> str:
        return os.path.join(self.output_dir, "gfa_annotation.tbl")

    @property
    def repeat_sample_table_path(self) -> str:
        return os.path.join(self.output_dir, "gfa_annotation_by_sample.tbl")

    @property
    def repeat_summary_path(self) -> str:
        return os.path.join(self.output_dir, "gfa_annotation_summary.tsv")

    @property
    def gff3_path(self) -> str:
        return os.path.join(self.output_dir, "gfa_node_annotation.gff3")

    @property
    def conflict_path(self) -> str:
        return os.path.join(self.output_dir, "gfa_node_annotation_conflicts.tsv")

    @property
    def unresolved_path(self) -> str:
        return os.path.join(self.output_dir, "gfa_annotation_unresolved.tsv")

    @property
    def annotated_gfa_path(self) -> str:
        return os.path.join(self.output_dir, "pangenome.annotated.gfa")

    @property
    def provenance_path(self) -> str:
        return os.path.join(self.output_dir, "gfa_annotation_provenance.json")

    def annotate(self) -> str:
        logger.info("Loading annotation manifest and source annotations")
        self.load_manifest_and_annotations()
        logger.info("Loaded %d manifest entries and %d annotation features",
                    len(self.entries), len(self.features))
        logger.info("Parsing graph nodes and sample occurrences")
        self.parse_gfa_occurrences()
        logger.info("Parsed %d nodes and %d occurrences (%s)",
                    len(self.nodes), len(self.occurrences), self.graph_mode)
        logger.info("Projecting genome-coordinate features to node coordinates")
        self.project_annotations()
        logger.info("Writing node, occurrence, RepeatMasker-style, GFF3, GFA and provenance outputs")
        self.write_outputs()
        logger.info(
            "Projected %d annotations through %d occurrences onto %d GFA nodes",
            len(self.projected), len(self.occurrences), len(self.nodes),
        )
        return self.node_table_path

    def load_manifest_and_annotations(self) -> None:
        if not os.path.isfile(self.annotation_manifest):
            raise FileNotFoundError(f"Annotation manifest not found: {self.annotation_manifest}")
        base = os.path.dirname(self.annotation_manifest)
        with _open_text(self.annotation_manifest) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if not reader.fieldnames or "sample" not in reader.fieldnames:
                raise ValueError("Annotation manifest requires a tab-delimited sample column")
            for number, row in enumerate(reader, 2):
                sample = (row.get("sample") or "").strip()
                haplotype = (row.get("haplotype") or "0").strip()
                if not sample:
                    raise ValueError(f"Manifest line {number} has an empty sample")
                def resolve(name: str) -> Optional[str]:
                    value = (row.get(name) or "").strip()
                    if not value:
                        return None
                    path = value if os.path.isabs(value) else os.path.join(base, value)
                    path = os.path.abspath(path)
                    if not os.path.isfile(path):
                        raise FileNotFoundError(f"Manifest line {number} {name} not found: {path}")
                    return path
                entry = ManifestEntry(
                    sample=sample, haplotype=haplotype, fasta=resolve("fasta"),
                    gff3=resolve("gff3"), te_annotation=resolve("te_annotation"),
                    seqid_map=resolve("seqid_map"),
                    path_prefix=(row.get("path_prefix") or "").strip(),
                    gaf=resolve("gaf"),
                    annotation_complete=_as_bool(row.get("annotation_complete") or ""),
                )
                if entry.key in self.entries:
                    raise ValueError(f"Duplicate manifest sample/haplotype: {entry.key}")
                if not entry.gff3 and not entry.te_annotation:
                    raise ValueError(f"Manifest entry {entry.key} has no GFF3 or TE annotation")
                self.entries[entry.key] = entry

        for entry in self.entries.values():
            if entry.seqid_map:
                self._read_seqid_map(entry)
            before = len(self.features)
            if entry.gff3:
                self._read_gff3(entry, entry.gff3, force_te=False)
            if entry.te_annotation:
                lower = entry.te_annotation.lower()
                for suffix in (".gz", ".bz2"):
                    if lower.endswith(suffix):
                        lower = lower[:-len(suffix)]
                if lower.endswith(".out"):
                    self._read_repeatmasker(entry, entry.te_annotation)
                elif lower.endswith(".bed"):
                    self._read_bed(entry, entry.te_annotation)
                else:
                    self._read_gff3(entry, entry.te_annotation, force_te=True)
            self._derive_introns(entry, before)

        seen = set()
        unique = []
        for feature in self.features:
            key = (feature.entry_key, feature.seqid, feature.start, feature.end,
                   feature.strand, feature.feature_type, feature.category,
                   feature.feature_id, feature.parent, feature.family, feature.te_class)
            if key not in seen:
                seen.add(key)
                unique.append(feature)
        self.features = unique
        for index, feature in enumerate(self.features):
            source_key = (feature.entry_key, feature.seqid)
            self.features_by_source[source_key].append(index)
            self.annotation_seqids[feature.entry_key].add(feature.seqid)
            for bin_id in range(feature.start // _BIN_SIZE, (feature.end - 1) // _BIN_SIZE + 1):
                self.feature_bins[(feature.entry_key, feature.seqid, bin_id)].append(index)

    def _read_gff3(self, entry: ManifestEntry, path: str, force_te: bool) -> None:
        with _open_text(path) as handle:
            for line_number, line in enumerate(handle, 1):
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) != 9:
                    raise ValueError(f"Malformed GFF3 {path}:{line_number}")
                try:
                    start, end = int(fields[3]) - 1, int(fields[4])
                except ValueError as exc:
                    raise ValueError(f"Invalid GFF3 coordinates {path}:{line_number}") from exc
                if start < 0 or end <= start:
                    raise ValueError(f"Invalid GFF3 interval {path}:{line_number}")
                attrs = _attributes(fields[8])
                lower = {key.lower(): value for key, value in attrs.items()}
                feature_type = fields[2]
                evidence = " ".join((fields[1], feature_type, fields[8]))
                is_te = force_te or bool(_TE_WORD_RE.search(evidence))
                category = "TE" if is_te else ("gene_structure" if feature_type.lower() in _GENE_TYPES else "other")
                feature_id = lower.get("id", lower.get("name", f"{feature_type}:{fields[0]}:{start}-{end}"))
                parent = lower.get("parent", "")
                family = ""
                te_class = ""
                if is_te:
                    family = next((lower[k] for k in ("family", "repeat_family", "name", "id", "target") if k in lower), feature_id)
                    te_class = next((lower[k] for k in ("classification", "repeat_class", "class", "superfamily") if k in lower), feature_type)
                self.features.append(Feature(
                    entry.key, self._normalize_seqid(entry, fields[0]), start, end, fields[6], fields[7],
                    feature_type, category, feature_id, parent, family, te_class,
                    os.path.basename(path), fields[8],
                ))

    def _read_repeatmasker(self, entry: ManifestEntry, path: str) -> None:
        with _open_text(path) as handle:
            for line in handle:
                fields = line.split()
                if len(fields) < 14:
                    continue
                try:
                    start, end = int(fields[5]) - 1, int(fields[6])
                except ValueError:
                    continue
                if start < 0 or end <= start:
                    continue
                strand = "-" if fields[8] == "C" else "+"
                is_te = not bool(_NON_TE_RM_CLASS_RE.match(fields[10]))
                self.features.append(Feature(
                    entry.key, self._normalize_seqid(entry, fields[4]), start, end, strand, ".",
                    "transposable_element" if is_te else fields[10],
                    "TE" if is_te else "other",
                    f"RM:{fields[4]}:{start}-{end}:{fields[9]}",
                    f"RMID:{fields[14]}" if len(fields) > 14 else "",
                    fields[9] if is_te else "", fields[10] if is_te else "",
                    os.path.basename(path), ";".join((
                        f"SW_score={fields[0]}", f"perc_div={fields[1]}",
                        f"perc_del={fields[2]}", f"perc_ins={fields[3]}",
                        f"RepeatID={fields[14] if len(fields) > 14 else ''}",
                    )),
                ))

    def _read_bed(self, entry: ManifestEntry, path: str) -> None:
        with _open_text(path) as handle:
            for line_number, line in enumerate(handle, 1):
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 3:
                    raise ValueError(f"Malformed BED {path}:{line_number}")
                start, end = int(fields[1]), int(fields[2])
                name = fields[3] if len(fields) > 3 else f"TE:{fields[0]}:{start}-{end}"
                family, _, te_class = name.partition("#")
                self.features.append(Feature(
                    entry.key, self._normalize_seqid(entry, fields[0]), start, end,
                    fields[5] if len(fields) > 5 and fields[5] in "+-" else ".", ".",
                    "transposable_element", "TE", name, "", family,
                    te_class or "Unknown", os.path.basename(path), "BED",
                ))

    def _read_seqid_map(self, entry: ManifestEntry) -> None:
        with _open_text(entry.seqid_map) as handle:
            for line_number, line in enumerate(handle, 1):
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 2:
                    raise ValueError(
                        f"Malformed seqid_map {entry.seqid_map}:{line_number}")
                if fields[0].lower() in {"annotation_seqid", "source_seqid"}:
                    continue
                previous = self.seqid_maps[entry.key].get(fields[0])
                if previous is not None and previous != fields[1]:
                    raise ValueError(
                        f"Conflicting seqid_map value for {entry.key} {fields[0]}")
                self.seqid_maps[entry.key][fields[0]] = fields[1]

    def _normalize_seqid(self, entry: ManifestEntry, seqid: str) -> str:
        """Normalize PanSN annotation seqids to the contig part used by P/W paths."""
        if seqid in self.seqid_maps.get(entry.key, {}):
            return self.seqid_maps[entry.key][seqid]
        prefix = f"{entry.sample}#{entry.haplotype}#"
        if seqid.startswith(prefix):
            return seqid[len(prefix):]
        if entry.path_prefix and seqid.startswith(entry.path_prefix):
            stripped = seqid[len(entry.path_prefix):].lstrip("#:/")
            return stripped or seqid
        return seqid

    def _derive_introns(self, entry: ManifestEntry, start_index: int) -> None:
        exons: Dict[Tuple[str, str], List[Feature]] = defaultdict(list)
        for feature in self.features[start_index:]:
            if feature.entry_key != entry.key or feature.feature_type.lower() != "exon":
                continue
            for parent in filter(None, feature.parent.split(",")):
                exons[(feature.seqid, parent)].append(feature)
        for (seqid, parent), values in exons.items():
            values.sort(key=lambda item: (item.start, item.end))
            for left, right in zip(values, values[1:]):
                if right.start <= left.end:
                    continue
                self.features.append(Feature(
                    entry.key, seqid, left.end, right.start, left.strand, ".", "intron",
                    "gene_structure", f"intron:{parent}:{left.end}-{right.start}", parent,
                    "", "", "derived_from_exons", f"Parent={parent}",
                ))

    def parse_gfa_occurrences(self) -> None:
        # Three streaming passes keep memory proportional to the graph, not to
        # its paths: segments, then links, then P/W traversals one at a time.
        # Results (occurrence order, unresolved order) match a single pass that
        # buffers every path and resolves P lines before W lines.
        last_path_line = self._read_gfa_segments()
        if not self.nodes:
            raise ValueError("GFA contains no segments")
        self._read_gfa_links()

        malformed: List[Tuple[str, str, str]] = []
        path_unresolved: List[Tuple[str, str, str]] = []
        walk_unresolved: List[Tuple[str, str, str]] = []
        deferred_walks: List[Tuple] = []
        before = self.unresolved
        try:
            with open(self.gfa_file, "rt") as handle:
                for line_number, line in enumerate(handle, 1):
                    head = line[:1]
                    if head != "P" and head != "W":
                        continue
                    fields = line.rstrip("\n").split("\t")
                    if fields[0] == "P" and len(fields) >= 3:
                        tokens = [token for token in re.split(r"[,;]", fields[2]) if token]
                        if not tokens or any(len(token) < 2 or token[-1] not in "+-"
                                             for token in tokens):
                            malformed.append((
                                f"line:{line_number}", "malformed_path_steps", fields[2]))
                            continue
                        steps = [(token[:-1], token[-1]) for token in tokens]
                        del tokens
                        overlaps = None if len(fields) < 4 or fields[3] == "*" else fields[3].split(",")
                        self.unresolved = path_unresolved
                        self._add_path(fields[1], steps, overlaps)
                    elif fields[0] == "W" and len(fields) >= 7:
                        parsed = _WALK_STEP_RE.findall(fields[6])
                        if not parsed or "".join(a + b for a, b in parsed) != fields[6]:
                            malformed.append((f"line:{line_number}", "malformed_walk", fields[6]))
                            continue
                        walk = (fields[1], fields[2], fields[3], fields[4],
                                [(node, "+" if sign == ">" else "-") for sign, node in parsed])
                        del parsed
                        if line_number < last_path_line:
                            deferred_walks.append(walk)
                            continue
                        self.unresolved = walk_unresolved
                        for item in deferred_walks:
                            self._add_walk(*item)
                        deferred_walks.clear()
                        self._add_walk(*walk)
            self.unresolved = walk_unresolved
            for item in deferred_walks:
                self._add_walk(*item)
        finally:
            self.unresolved = before
        self.unresolved.extend(malformed)
        self.unresolved.extend(path_unresolved)
        self.unresolved.extend(walk_unresolved)

        self.graph_mode = "PGGB_paths" if self.occurrences else "rGFA"
        self._add_gaf_occurrences()
        self._add_rgfa_origin_occurrences()
        self.link_keys = np.zeros(0, dtype=np.int64)
        self.link_cigars = np.zeros(0, dtype=np.int32)

    def _read_gfa_segments(self) -> int:
        """Load S lines; return the line number of the last P line."""
        last_path_line = 0
        with open(self.gfa_file, "rt") as handle:
            for line_number, line in enumerate(handle, 1):
                if line.startswith("P\t"):
                    last_path_line = line_number
                    continue
                if line[:1] != "S":
                    continue
                fields = line.rstrip("\n").split("\t")
                if fields[0] != "S":
                    continue
                if len(fields) < 3:
                    raise ValueError(f"Malformed S line {line_number}")
                if fields[1] in self.nodes:
                    raise ValueError(f"Duplicate GFA segment {fields[1]} at line {line_number}")
                tags = {tag[:2]: tag[5:] for tag in fields[3:] if len(tag) >= 6 and tag[2:5] in {":i:", ":Z:"}}
                sequence = None if fields[2] == "*" else fields[2]
                length = len(sequence) if sequence is not None else int(tags.get("LN", 0))
                if length <= 0:
                    raise ValueError(
                        f"GFA segment {fields[1]} has no sequence or positive LN")
                node = Node(fields[1], sequence, length, tags or _NO_TAGS, len(self.node_ids))
                self.nodes[node.node_id] = node
                self.node_ids.append(node.node_id)
        return last_path_line

    def _read_gfa_links(self) -> None:
        keys = array("q")
        cigars = array("i")
        cigar_codes: Dict[str, int] = {}
        with open(self.gfa_file, "rt") as handle:
            for line in handle:
                if line[:1] != "L":
                    continue
                fields = line.rstrip("\n").split("\t")
                if fields[0] != "L" or len(fields) < 6:
                    continue
                # Traversals only ever look links up with +/- orientations.
                if fields[2] not in ("+", "-") or fields[4] not in ("+", "-"):
                    continue
                keys.append(self._link_key(
                    self._node_ref(fields[1], create=True), fields[2],
                    self._node_ref(fields[3], create=True), fields[4]))
                code = cigar_codes.get(fields[5])
                if code is None:
                    code = cigar_codes[fields[5]] = len(self.link_cigar_values)
                    self.link_cigar_values.append(fields[5])
                cigars.append(code)
        key_array = _column(keys, np.int64)
        order = np.argsort(key_array, kind="stable")
        sorted_keys = key_array[order]
        del key_array, keys
        # A repeated link keeps its last CIGAR, as a dict assignment would.
        last = np.ones(len(sorted_keys), dtype=bool)
        last[:-1] = sorted_keys[:-1] != sorted_keys[1:]
        self.link_keys = sorted_keys[last]
        self.link_cigars = _column(cigars, np.int32)[order][last]

    def _node_ref(self, node_id: str, create: bool = False) -> int:
        node = self.nodes.get(node_id)
        if node is not None:
            return node.index
        index = self.extra_node_index.get(node_id)
        if index is None:
            if not create:
                return -1
            index = self.extra_node_index[node_id] = len(self.node_ids)
            self.node_ids.append(node_id)
        return index

    @staticmethod
    def _link_key(left: int, left_orientation: str, right: int,
                  right_orientation: str) -> int:
        if left >= 1 << 30 or right >= 1 << 30:
            raise ValueError("GFA has too many segments for packed link keys")
        return (((left << 1) | (left_orientation == "-")) << 32) | (
            (right << 1) | (right_orientation == "-"))

    def _link_get(self, left_id: str, left_orientation: str, right_id: str,
                  right_orientation: str) -> Optional[str]:
        left, right = self._node_ref(left_id), self._node_ref(right_id)
        if left < 0 or right < 0 or not len(self.link_keys):
            return None
        key = self._link_key(left, left_orientation, right, right_orientation)
        position = int(np.searchsorted(self.link_keys, key))
        if position < len(self.link_keys) and int(self.link_keys[position]) == key:
            return self.link_cigar_values[int(self.link_cigars[position])]
        return None

    def _path_link_cigars(self, path_name: str,
                          steps: Sequence[Tuple[str, str]]) -> Optional[List[str]]:
        """Return the link CIGAR between consecutive steps, or None if one is missing."""
        if len(steps) < 2:
            return []
        refs = np.fromiter((self._node_ref(node_id) for node_id, _o in steps),
                           dtype=np.int64, count=len(steps))
        reverse = np.fromiter((orientation == "-" for _n, orientation in steps),
                              dtype=np.int64, count=len(steps))
        known = (refs[:-1] >= 0) & (refs[1:] >= 0)
        code = (refs << 1) | reverse
        flipped = (refs << 1) | (1 - reverse)
        found = np.zeros(len(steps) - 1, dtype=bool)
        chosen = np.zeros(len(steps) - 1, dtype=np.int64)
        # GFA links are bidirected: A+ -> B+ is equivalent to B- -> A-.  PGGB
        # paths legitimately use either view; the forward view wins.
        for query in ((code[:-1] << 32) | code[1:], (flipped[1:] << 32) | flipped[:-1]):
            if not len(self.link_keys):
                break
            position = np.minimum(np.searchsorted(self.link_keys, query),
                                  len(self.link_keys) - 1)
            hit = known & ~found & (self.link_keys[position] == query)
            chosen[hit] = self.link_cigars[position[hit]]
            found |= hit
        if not found.all():
            index = int(np.argmin(found))
            left, right = steps[index], steps[index + 1]
            self.unresolved.append((
                f"P:{path_name}", "missing_path_link",
                f"index={index};{left[0]}{left[1]}->{right[0]}{right[1]}",
            ))
            return None
        values = self.link_cigar_values
        return [values[code] for code in chosen.tolist()]

    def _add_path(self, path_name: str, steps: List[Tuple[str, str]],
                  overlaps: Optional[List[str]]) -> None:
        resolved = self._resolve_embedded(path_name)
        if not resolved:
            self.unresolved.append((f"P:{path_name}", "manifest_path_unmatched", path_name))
            return
        entry, seqid = resolved
        if overlaps is None:
            overlaps = self._path_link_cigars(path_name, steps)
            if overlaps is None:
                return
        self._add_embedded_traversal(entry, seqid, f"P:{path_name}", 0, steps, overlaps)

    def _add_walk(self, sample: str, haplotype: str, seqid: str, start_text: str,
                  steps: List[Tuple[str, str]]) -> None:
        entry = self.entries.get((sample, haplotype))
        if not entry:
            candidates = [x for x in self.entries.values() if x.sample == sample]
            entry = candidates[0] if len(candidates) == 1 else None
        if not entry or start_text == "*":
            self.unresolved.append((f"W:{sample}#{haplotype}#{seqid}", "manifest_or_start_unresolved", start_text))
            return
        self._add_embedded_traversal(
            entry, seqid, f"W:{sample}#{haplotype}#{seqid}:{start_text}", int(start_text),
            steps, ["0M"] * max(0, len(steps) - 1),
        )

    def _resolve_embedded(self, path_name: str) -> Optional[Tuple[ManifestEntry, str]]:
        parts = path_name.split("#")
        if len(parts) >= 3 and (parts[0], parts[1]) in self.entries:
            return self.entries[(parts[0], parts[1])], "#".join(parts[2:])
        matches = []
        for entry in self.entries.values():
            if entry.path_prefix and path_name.startswith(entry.path_prefix):
                seqid = path_name[len(entry.path_prefix):].lstrip("#:/")
                matches.append((entry, seqid))
        if len(matches) == 1:
            return matches[0]
        if len(self.entries) == 1:
            return next(iter(self.entries.values())), path_name
        return None

    def _resolve_source(self, source: str, allowed: Optional[Set[Tuple[str, str]]] = None) -> Optional[Tuple[ManifestEntry, str]]:
        parts = source.split("#")
        if len(parts) >= 3 and (parts[0], parts[1]) in self.entries:
            entry = self.entries[(parts[0], parts[1])]
            if allowed is None or entry.key in allowed:
                return entry, "#".join(parts[2:])
        candidates = []
        for entry in self.entries.values():
            if allowed is not None and entry.key not in allowed:
                continue
            if source in self.annotation_seqids.get(entry.key, set()):
                candidates.append((entry, source))
            elif entry.path_prefix and source.startswith(entry.path_prefix):
                candidates.append((entry, source[len(entry.path_prefix):].lstrip("#:/")))
        return candidates[0] if len(candidates) == 1 else None

    def _add_embedded_traversal(self, entry: ManifestEntry, seqid: str, path_id: str,
                                source_start: int, steps: Sequence[Tuple[str, str]],
                                overlaps: Sequence[Optional[str]]) -> None:
        if len(overlaps) != max(0, len(steps) - 1):
            self.unresolved.append((path_id, "overlap_count_mismatch", str(len(overlaps))))
            return
        overlap_lengths = []
        parsed: Dict[Optional[str], Optional[int]] = {}
        for cigar in overlaps:
            if cigar in parsed:
                overlap = parsed[cigar]
            else:
                overlap = parsed[cigar] = _simple_overlap(cigar)
            if overlap is None:
                self.unresolved.append((path_id, "complex_or_missing_overlap", str(cigar)))
                return
            overlap_lengths.append(overlap)
        layout = []
        missing = None
        cursor = source_start
        for index, (node_id, orientation) in enumerate(steps):
            node = self.nodes.get(node_id)
            if node is None:
                missing = node_id
                break
            overlap = overlap_lengths[index - 1] if index else 0
            start = cursor - overlap
            end = start + node.length
            layout.append((node, orientation, start, end))
            cursor = end
        del overlap_lengths

        # Every key carries (source, genome_start), so this traversal can only
        # duplicate an earlier occurrence of the same source whose genome span
        # intersects its own, or one of its own steps when overlaps let starts
        # repeat.  Otherwise rows are appended without the global key set.
        source = (entry.key, seqid)
        spans = self._source_spans.setdefault(source, [])
        if layout:
            low = min(start for _n, _o, start, _e in layout)
            high = max(end for _n, _o, _s, end in layout)
            trusted = self.occurrence_keys is None and not any(
                low < span_high and span_low < high for span_low, span_high in spans)
            spans.append((low, high))
        else:
            trusted = True
        if not trusted:
            for index, (node, orientation, start, end) in enumerate(layout):
                self._add_occurrence(entry, seqid, node.node_id, index, start, end, 0,
                                     node.length, orientation, path_id,
                                     "embedded_path", "exact")
        else:
            local = set() if any(parsed.values()) else None
            append = self.occurrences.append
            for index, (node, orientation, start, end) in enumerate(layout):
                if local is not None:
                    key = (node.node_id, start, end, orientation)
                    if key in local:
                        continue
                    local.add(key)
                append(entry.key, seqid, node.index, index, start, end, 0, node.length,
                       orientation, path_id, "embedded_path", "exact")
        if missing is not None:
            self.unresolved.append((path_id, "missing_node", missing))
            return
        self.embedded_sources.add((entry.key, seqid))

    def _add_occurrence(self, entry: ManifestEntry, seqid: str, node_id: str, index: int,
                        genome_start: int, genome_end: int, node_start: int, node_end: int,
                        orientation: str, path_id: str, evidence: str,
                        coordinate_status: str) -> None:
        if genome_end <= genome_start or node_end <= node_start:
            return
        if self.occurrence_keys is None:
            occurrences = self.occurrences
            self.occurrence_keys = {occurrences.key(i) for i in range(len(occurrences))}
        key = (entry.key, seqid, node_id, genome_start, genome_end, node_start, node_end, orientation)
        if key in self.occurrence_keys:
            return
        self.occurrence_keys.add(key)
        self.occurrences.append(
            entry.key, seqid, self._node_ref(node_id, create=True), index,
            genome_start, genome_end, node_start, node_end, orientation, path_id,
            evidence, coordinate_status,
        )

    def _add_gaf_occurrences(self) -> None:
        by_path: Dict[str, Set[Tuple[str, str]]] = defaultdict(set)
        for entry in self.entries.values():
            gaf = entry.gaf or self.default_gaf
            if gaf:
                by_path[gaf].add(entry.key)
        strict_rgfa = bool(by_path) and all(
            node.sequence is not None and {"SN", "SO", "SR"}.issubset(node.tags)
            for node in self.nodes.values()
        )
        if strict_rgfa:
            self._add_minigraph_gaf_occurrences(by_path)
            return
        for gaf, allowed in by_path.items():
            with _open_text(gaf) as handle:
                for line_number, line in enumerate(handle, 1):
                    if not line.strip() or line.startswith("#"):
                        continue
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) < 12:
                        self.unresolved.append((f"{gaf}:{line_number}", "malformed_gaf", "<12 columns"))
                        continue
                    try:
                        qlen, qstart, qend = int(fields[1]), int(fields[2]), int(fields[3])
                        plen = int(fields[6])
                        pstart, pend = int(fields[7]), int(fields[8])
                        matches, block, mapq = int(fields[9]), int(fields[10]), int(fields[11])
                    except ValueError:
                        self.unresolved.append((f"{gaf}:{line_number}", "invalid_gaf_number", ""))
                        continue
                    if not (0 <= qstart < qend <= qlen and
                            0 <= pstart < pend <= plen and block > 0 and mapq >= 0):
                        self.unresolved.append((
                            f"{gaf}:{line_number}", "invalid_gaf_coordinates", ""))
                        continue
                    identity = matches / block if block else 0.0
                    if mapq < self.min_mapq or identity < self.min_identity:
                        continue
                    resolved = self._resolve_source(fields[0], allowed)
                    if not resolved:
                        self.unresolved.append((f"{gaf}:{line_number}", "gaf_query_unmatched", fields[0]))
                        continue
                    entry, seqid = resolved
                    parsed = _WALK_STEP_RE.findall(fields[5])
                    if not parsed or "".join(a + b for a, b in parsed) != fields[5]:
                        self.unresolved.append((f"{gaf}:{line_number}", "named_or_malformed_gaf_path", fields[5]))
                        continue
                    qspan, pspan = qend - qstart, pend - pstart
                    tags = {token[:2]: token[5:] for token in fields[12:] if len(token) >= 6 and token[2:5] in {":Z:", ":i:"}}
                    inferred = qspan == pspan == block
                    gapless = qspan == pspan and (
                        _gapless_cigar(tags.get("cg", ""), qspan) or inferred)
                    if not gapless:
                        self.unresolved.append((f"{gaf}:{line_number}", "gapped_gaf_requires_base_transform", tags.get("cg", "")))
                        continue
                    missing_node = next((node_id for _sign, node_id in parsed
                                         if node_id not in self.nodes), None)
                    if missing_node:
                        self.unresolved.append(
                            (f"{gaf}:{line_number}", "missing_node", missing_node))
                        continue
                    raw_length = sum(self.nodes[node_id].length for _sign, node_id in parsed)
                    layout = []
                    cursor = 0
                    layout_failed = False
                    for index, (sign, node_id) in enumerate(parsed):
                        orientation = "+" if sign == ">" else "-"
                        overlap = 0
                        if index:
                            previous_sign, previous_id = parsed[index - 1]
                            previous_orientation = "+" if previous_sign == ">" else "-"
                            cigar = self._link_get(
                                previous_id, previous_orientation, node_id, orientation)
                            if cigar is None:
                                cigar = self._link_get(
                                    node_id, _flip(orientation), previous_id,
                                    _flip(previous_orientation),
                                )
                            if cigar is None and plen != raw_length:
                                self.unresolved.append((
                                    f"{gaf}:{line_number}", "missing_gaf_path_link",
                                    f"index={index - 1};{previous_id}{previous_orientation}"
                                    f"->{node_id}{orientation}",
                                ))
                                layout_failed = True
                                break
                            if cigar is not None:
                                parsed_overlap = _simple_overlap(cigar)
                                if parsed_overlap is None:
                                    self.unresolved.append((
                                        f"{gaf}:{line_number}", "complex_gaf_path_overlap", cigar))
                                    layout_failed = True
                                    break
                                overlap = parsed_overlap
                        node = self.nodes[node_id]
                        left, right = cursor - overlap, cursor - overlap + node.length
                        cursor = right
                        layout.append((index, sign, node_id, left, right))
                    if layout_failed:
                        continue
                    if cursor != plen:
                        self.unresolved.append((
                            f"{gaf}:{line_number}", "gaf_path_length_mismatch",
                            f"computed={cursor};declared={plen}",
                        ))
                        continue
                    for index, sign, node_id, left, right in layout:
                        node = self.nodes[node_id]
                        overlap_start, overlap_end = max(left, pstart), min(right, pend)
                        if overlap_end <= overlap_start:
                            continue
                        if fields[4] == "+":
                            genome_start = qstart + overlap_start - pstart
                            genome_end = qstart + overlap_end - pstart
                            orientation = "+" if sign == ">" else "-"
                        else:
                            genome_start = qend - (overlap_end - pstart)
                            genome_end = qend - (overlap_start - pstart)
                            orientation = "-" if sign == ">" else "+"
                        if sign == ">":
                            node_start, node_end = overlap_start - left, overlap_end - left
                        else:
                            node_start, node_end = right - overlap_end, right - overlap_start
                        self._add_occurrence(
                            entry, seqid, node_id, index, genome_start, genome_end,
                            node_start, node_end, orientation, f"GAF:{os.path.basename(gaf)}:{line_number}",
                            "GAF", "gapless_cigar" if tags.get("cg") else "gapless_length_inferred",
                        )
                    if qstart == 0 and qend == qlen:
                        # Still not used to infer absence: GAF may omit unmapped contigs.
                        pass
            if self.occurrences and self.graph_mode == "rGFA":
                self.graph_mode = "rGFA_GAF"

    def _add_minigraph_gaf_occurrences(
            self, by_path: Dict[str, Set[Tuple[str, str]]]) -> None:
        """Project strict minigraph rGFA/GAF alignments through cg blocks.

        Real minigraph GAF paths may be named source paths or walks whose tokens
        are ``sample#hap#contig:start-end`` intervals rather than GFA segment
        IDs.  The shared rGFA parser resolves both forms and validates every
        path edge.  Each match CIGAR operation is then intersected with the
        traversed nodes; query-only and path-only gaps advance only their own
        coordinate system and therefore cannot shift downstream annotations.
        """
        try:
            segments, rgfa_links, coordinates, _references, source_intervals = \
                load_rgfa_graph(self.gfa_file)
        except RgfaGafError as exc:
            raise ValueError(f"Cannot index minigraph rGFA: {exc}") from exc

        for gaf, allowed in by_path.items():
            try:
                records = iter_gaf_records(
                    gaf, segments, rgfa_links, coordinates, source_intervals,
                    require_primary=True, require_cigar=True,
                    stats=self.gaf_stats,
                )
                for record, steps in records:
                    if record.mapq < self.min_mapq or record.identity < self.min_identity:
                        self.gaf_stats["records_skipped_low_quality"] += 1
                        continue
                    if record.path.startswith((">", "<")):
                        if re.search(r"[><][^><]+:\d+-\d+", record.path):
                            self.gaf_stats["records_coordinate_walk"] += 1
                        else:
                            self.gaf_stats["records_node_walk"] += 1
                    else:
                        self.gaf_stats["records_named_path"] += 1
                    resolved = self._resolve_source(record.query, allowed)
                    if not resolved:
                        self.gaf_stats["records_skipped_query_unmatched"] += 1
                        self.unresolved.append((
                            f"{gaf}:{record.record_id}", "gaf_query_unmatched", record.query))
                        continue
                    entry, seqid = resolved
                    before = len(self.occurrences)
                    self._add_cigar_block_occurrences(
                        entry, seqid, gaf, record, steps, segments)
                    self.gaf_stats["records_projected"] += 1
                    self.gaf_stats["occurrences_added"] += len(self.occurrences) - before
            except RgfaGafError as exc:
                raise ValueError(f"Invalid minigraph GAF {gaf}: {exc}") from exc
        if self.occurrences and self.graph_mode == "rGFA":
            self.graph_mode = "rGFA_GAF"

    def _add_cigar_block_occurrences(
            self, entry: ManifestEntry, seqid: str, gaf: str, record,
            steps: Sequence, segments: Dict) -> None:
        cigar = record.tags.get("cg", "")
        parsed = list(_CIGAR_RE.finditer(cigar))
        if not parsed or "".join(item.group(0) for item in parsed) != cigar:
            self.unresolved.append((
                f"{gaf}:{record.record_id}", "malformed_cg", cigar))
            return
        operations: List[Tuple[int, str]] = []
        for item in parsed:
            length, operation = int(item.group(1)), item.group(2)
            if length <= 0:
                self.unresolved.append((
                    f"{gaf}:{record.record_id}", "zero_length_cg_operation",
                    item.group(0)))
                return
            # M, = and X all advance query and graph coordinates one-for-one.
            # Adjacent runs can therefore be merged without changing a single
            # mapped coordinate; only I/S or D/N break linearity.
            if (operations and operation in {"M", "=", "X"} and
                    operations[-1][1] in {"M", "=", "X"}):
                operations[-1] = (operations[-1][0] + length, "M")
            else:
                operations.append((length, operation))

        layout = []
        path_cursor = 0
        for index, step in enumerate(steps):
            node = segments[step.node_id]
            layout.append((index, step.node_id, step.orientation,
                           path_cursor, path_cursor + node.length))
            path_cursor += node.length
        if path_cursor != record.path_length:
            # iter_gaf_records already enforces this; retain a local guard so
            # this coordinate-transform routine remains safe if reused.
            self.unresolved.append((
                f"{gaf}:{record.record_id}", "gaf_path_length_mismatch",
                f"computed={path_cursor};declared={record.path_length}"))
            return

        query_consumed = 0
        target_cursor = record.path_start
        layout_index = 0
        path_id = f"GAF:{os.path.basename(gaf)}:{record.record_id}"
        for length, operation in operations:
            if operation in {"M", "=", "X"}:
                self.gaf_stats["cigar_match_runs"] += 1
                block_start, block_end = target_cursor, target_cursor + length
                while layout_index < len(layout) and layout[layout_index][4] <= block_start:
                    layout_index += 1
                node_index = layout_index
                while node_index < len(layout) and layout[node_index][3] < block_end:
                    index, node_id, step_orientation, node_left, node_right = layout[node_index]
                    overlap_start = max(block_start, node_left)
                    overlap_end = min(block_end, node_right)
                    if overlap_end > overlap_start:
                        delta = overlap_start - block_start
                        span = overlap_end - overlap_start
                        if record.strand == "+":
                            genome_start = record.query_start + query_consumed + delta
                            genome_end = genome_start + span
                            orientation = step_orientation
                        else:
                            genome_end = record.query_end - query_consumed - delta
                            genome_start = genome_end - span
                            orientation = _flip(step_orientation)
                        if step_orientation == "+":
                            node_start = overlap_start - node_left
                            node_end = overlap_end - node_left
                        else:
                            node_start = node_right - overlap_end
                            node_end = node_right - overlap_start
                        self._add_occurrence(
                            entry, seqid, node_id, index, genome_start, genome_end,
                            node_start, node_end, orientation, path_id, "GAF",
                            "cigar_block_exact",
                        )
                    node_index += 1
            if operation in {"M", "=", "X", "I", "S"}:
                query_consumed += length
            if operation in {"M", "=", "X", "D", "N"}:
                target_cursor += length

        if (query_consumed != record.query_end - record.query_start or
                target_cursor != record.path_end):
            # This should be unreachable because the shared parser validates
            # both spans before yielding a record.
            raise ValueError(
                f"Internal GAF CIGAR transform mismatch at {gaf}:{record.record_id}")

    def _add_rgfa_origin_occurrences(self) -> None:
        for node in self.nodes.values():
            source, offset = node.tags.get("SN"), node.tags.get("SO")
            if not source or offset is None:
                continue
            resolved = self._resolve_source(source)
            if not resolved:
                continue
            try:
                start = int(offset)
            except ValueError:
                self.unresolved.append((node.node_id, "invalid_rgfa_SO", offset))
                continue
            entry, seqid = resolved
            self._add_occurrence(entry, seqid, node.node_id, 0, start, start + node.length,
                                 0, node.length, "+", f"rGFA:{source}:{start}",
                                 "rGFA_origin", "origin_only")

    def project_annotations(self) -> None:
        occurrences, projected = self.occurrences, self.projected
        sources = occurrences.sources.values
        features, feature_bins = self.features, self.feature_bins
        total = len(occurrences)
        for o in range(total):
            entry_key, seqid = sources[occurrences.source[o]]
            genome_start, genome_end = occurrences.genome_start[o], occurrences.genome_end[o]
            candidates: Set[int] = set()
            for bin_id in range(genome_start // _BIN_SIZE,
                                (genome_end - 1) // _BIN_SIZE + 1):
                candidates.update(feature_bins.get((entry_key, seqid, bin_id), ()))
            if candidates:
                node_start, node_end = occurrences.node_start[o], occurrences.node_end[o]
                reverse = occurrences.reverse[o]
                for index in sorted(candidates):
                    feature = features[index]
                    overlap_start = max(feature.start, genome_start)
                    overlap_end = min(feature.end, genome_end)
                    if overlap_end <= overlap_start:
                        continue
                    left = overlap_start - genome_start
                    right = overlap_end - genome_start
                    projected.occurrence.append(o)
                    projected.feature.append(index)
                    if reverse:
                        projected.node_start.append(node_end - right)
                        projected.node_end.append(node_end - left)
                    else:
                        projected.node_start.append(node_start + left)
                        projected.node_end.append(node_start + right)
                    projected.genome_start.append(overlap_start)
                    projected.genome_end.append(overlap_end)
            if (o + 1) % 1_000_000 == 0:
                logger.info("Projection progress: %d/%d occurrences; %d intervals",
                            o + 1, total, len(projected))

    def _entry_complete(self, occurrence: Occurrence) -> bool:
        return self.entries[occurrence.entry_key].annotation_complete

    def _node_sort_rank(self) -> np.ndarray:
        """Rank of every node id (by node index) in Python string order."""
        ids = self.node_ids
        rank = np.empty(len(ids), dtype=np.int64)
        rank[np.fromiter(sorted(range(len(ids)), key=ids.__getitem__),
                         dtype=np.int64, count=len(ids))] = np.arange(len(ids))
        return rank

    def _projection_columns(self) -> Dict[str, np.ndarray]:
        projected, occurrences = self.projected, self.occurrences
        occurrence = _column(projected.occurrence, np.uint32)
        return {
            "occurrence": occurrence,
            "feature": _column(projected.feature, np.uint32),
            "node": _column(occurrences.node, np.uint32)[occurrence],
            "node_start": _column(projected.node_start, np.int64),
            "node_end": _column(projected.node_end, np.int64),
        }

    @staticmethod
    def _chunks(total: int) -> Iterator[Tuple[int, int]]:
        for start in range(0, total, _CHUNK):
            yield start, min(total, start + _CHUNK)

    def write_outputs(self) -> None:
        occurrences, projected = self.occurrences, self.projected
        columns = self._projection_columns()
        node_rank = self._node_sort_rank()
        rank_rows = node_rank[columns["node"]]
        feature_types = sorted({feature.feature_type for feature in self.features})
        type_rank = {value: rank for rank, value in enumerate(feature_types)}
        feature_type_rank = np.fromiter(
            (type_rank[feature.feature_type] for feature in self.features),
            dtype=np.int32, count=len(self.features))
        # Same order as a stable sort on (node_id, node_start, node_end,
        # occurrence_id, feature_type); occurrence ids are zero-padded row numbers.
        sorted_projected = np.lexsort((
            np.arange(len(projected)), feature_type_rank[columns["feature"]],
            columns["occurrence"], columns["node_end"], columns["node_start"], rank_rows,
        ))

        interval_fields = [
            "node_id", "node_start", "node_end", "overlap_bp", "occurrence_id",
            "sample", "haplotype", "seqid", "genome_start", "genome_end",
            "orientation", "category", "feature_type", "feature_id", "parent",
            "te_family", "te_class", "source", "strand", "phase", "attributes",
        ]
        with open(self.interval_table_path, "wt", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=interval_fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for start, end in self._chunks(len(sorted_projected)):
                for index in sorted_projected[start:end].tolist():
                    row = projected[index]
                    writer.writerow({
                        "node_id": row.node_id, "node_start": row.node_start, "node_end": row.node_end,
                        "overlap_bp": row.node_end - row.node_start, "occurrence_id": row.occurrence_id,
                        "sample": row.sample, "haplotype": row.haplotype, "seqid": row.seqid,
                        "genome_start": row.genome_start, "genome_end": row.genome_end,
                        "orientation": row.orientation, "category": row.category,
                        "feature_type": row.feature_type, "feature_id": row.feature_id,
                        "parent": row.parent, "te_family": row.family, "te_class": row.te_class,
                        "source": row.source, "strand": row.strand, "phase": row.phase,
                        "attributes": _clean(row.attributes),
                    })

        # Projections were appended in occurrence order: each occurrence owns a
        # contiguous run, in the order the original per-occurrence lists had.
        occurrence_offsets = _group_offsets(columns["occurrence"], len(occurrences))
        occurrence_fields = [
            "occurrence_id", "node_id", "sample", "haplotype", "seqid",
            "genome_start", "genome_end", "node_start", "node_end", "orientation",
            "path_id", "evidence", "coordinate_status", "annotation_callable",
            "identity", "te_covered_bp", "te_fraction", "feature_types",
        ]
        conflicts = []
        conflict_nodes = set()
        with open(self.occurrence_table_path, "wt", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=occurrence_fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for start, end in self._chunks(len(occurrences)):
                offsets = occurrence_offsets[start:end + 1].tolist()
                for position, o in enumerate(range(start, end)):
                    occurrence = occurrences[o]
                    rows = [projected[j] for j in range(offsets[position], offsets[position + 1])]
                    te_intervals = [(x.node_start, x.node_end) for x in rows if x.category == "TE"]
                    covered = _covered_bp(te_intervals)
                    span = occurrence.node_end - occurrence.node_start
                    fraction = covered / span if span else 0.0
                    callable_ = self._entry_complete(occurrence) or bool(te_intervals)
                    identity = "TE" if fraction >= self.te_min_fraction else (
                        "PartialTE" if fraction >= self.partial_min_fraction else (
                            "NonTE" if callable_ else "Unknown"))
                    values = {
                        "occurrence_id": occurrence.occurrence_id, "node_id": occurrence.node_id,
                        "sample": occurrence.sample, "haplotype": occurrence.haplotype,
                        "seqid": occurrence.seqid, "genome_start": occurrence.genome_start,
                        "genome_end": occurrence.genome_end, "node_start": occurrence.node_start,
                        "node_end": occurrence.node_end, "orientation": occurrence.orientation,
                        "path_id": occurrence.path_id, "evidence": occurrence.evidence,
                        "coordinate_status": occurrence.coordinate_status,
                        "annotation_callable": int(callable_), "identity": identity,
                        "te_covered_bp": covered, "te_fraction": f"{fraction:.6f}",
                        "feature_types": ",".join(sorted({x.feature_type for x in rows})) or "None",
                    }
                    writer.writerow(values)
                    cds = [(x.node_start, x.node_end) for x in rows if x.feature_type.lower() == "cds"]
                    overlap = _overlap_bp(te_intervals, cds)
                    if overlap:
                        conflict_nodes.add(occurrence.node_id)
                        conflicts.append((occurrence.node_id, occurrence.occurrence_id,
                                          occurrence.sample, "TE_CDS_overlap", overlap))
        del occurrence_offsets

        groups = len(self.node_ids)
        occurrence_node = _column(occurrences.node, np.uint32)
        occurrence_by_node = np.argsort(occurrence_node, kind="stable")
        occurrence_node_offsets = _group_offsets(occurrence_node[occurrence_by_node], groups)
        projected_by_node = np.argsort(columns["node"], kind="stable")
        projected_node_offsets = _group_offsets(columns["node"][projected_by_node], groups)

        samples = sorted({entry.sample for entry in self.entries.values()})
        summary_fields = [
            "node_id", "length", "sequence_status", "identity", "is_te",
            "te_covered_bp", "te_fraction", "te_occurrence_support",
            "te_genome_support", "annotated_occurrences", "total_occurrences",
            "callable_genomes", "primary_te_family", "primary_te_class",
            "feature_types", "gene_structure", "conflict", "confidence",
        ]
        sample_fields = [f"sample:{sample}" for sample in samples]
        self.node_summaries = [None] * len(self.nodes)
        with open(self.node_table_path, "wt", newline="") as handle, \
                open(self.matrix_path, "wt") as matrix_handle:
            writer = csv.DictWriter(handle, fieldnames=summary_fields + sample_fields,
                                    delimiter="\t", lineterminator="\n")
            writer.writeheader()
            matrix_handle.write("node_id\tsample\tstate\n")
            for node_id, node in self.nodes.items():
                group = node.index
                rows = [projected[j] for j in projected_by_node[
                    projected_node_offsets[group]:projected_node_offsets[group + 1]].tolist()]
                occurrences_here = [occurrences[j] for j in occurrence_by_node[
                    occurrence_node_offsets[group]:occurrence_node_offsets[group + 1]].tolist()]
                te_rows = [x for x in rows if x.category == "TE"]
                te_intervals = [(x.node_start, x.node_end) for x in te_rows]
                te_covered = min(node.length, _covered_bp(te_intervals))
                te_fraction = te_covered / node.length if node.length else 0.0
                # A gapped GAF traversal can map one biological node context
                # through several linear occurrence blocks.  Support must be
                # counted per traversal context, otherwise harmless =/X/M
                # fragmentation or an indel changes TE/Mixed classification.
                contexts: Dict[Tuple, List[Occurrence]] = defaultdict(list)
                occurrence_context: Dict[str, Tuple] = {}
                for occurrence in occurrences_here:
                    context = (
                        occurrence.entry_key, occurrence.seqid, occurrence.node_id,
                        occurrence.occurrence_index, occurrence.path_id,
                    )
                    contexts[context].append(occurrence)
                    occurrence_context[occurrence.occurrence_id] = context
                te_contexts = {
                    occurrence_context[row.occurrence_id] for row in te_rows
                    if row.occurrence_id in occurrence_context
                }
                callable_contexts = {
                    context for context, blocks in contexts.items()
                    if context in te_contexts or any(self._entry_complete(x) for x in blocks)
                }
                te_occ_support = (
                    len(te_contexts) / len(callable_contexts) if callable_contexts else 0.0)
                callable_samples = {
                    contexts[context][0].sample for context in callable_contexts}
                te_samples = {x.sample for x in te_rows}
                te_genome_support = len(te_samples) / len(callable_samples) if callable_samples else 0.0
                feature_types = sorted({x.feature_type for x in rows})
                gene_types = sorted({x.feature_type for x in rows if x.category == "gene_structure"})
                disagreement = bool(te_contexts) and any(
                    context not in te_contexts and
                    any(self._entry_complete(x) for x in blocks)
                    for context, blocks in contexts.items())
                if node_id in conflict_nodes:
                    identity = "Conflict"
                elif te_rows and (gene_types or disagreement):
                    identity = "Mixed"
                elif te_fraction >= self.te_min_fraction and te_occ_support >= 0.5:
                    identity = "TE"
                elif te_fraction >= self.partial_min_fraction or te_rows:
                    identity = "PartialTE"
                elif gene_types:
                    identity = "GeneStructure"
                elif callable_contexts:
                    identity = "NonTE"
                else:
                    identity = "Unknown"
                family_coverage: Dict[Tuple[str, str], List[Tuple[int, int]]] = defaultdict(list)
                for row in te_rows:
                    family_coverage[(row.family or "Unknown", row.te_class or "Unknown")].append((row.node_start, row.node_end))
                labels = sorted(family_coverage, key=lambda label: (-_covered_bp(family_coverage[label]), label))
                primary_family, primary_class = labels[0] if labels else ("None", "None")
                confidence = "low"
                if identity in {"TE", "NonTE", "GeneStructure"} and len(callable_samples) >= 2:
                    confidence = "high"
                elif identity not in {"Unknown", "Conflict"}:
                    confidence = "medium"
                sample_states = {}
                for sample in samples:
                    sample_occ = [x for x in occurrences_here if x.sample == sample]
                    if not sample_occ:
                        entry_keys = {entry.key for entry in self.entries.values() if entry.sample == sample}
                        node_seqids = {item.seqid for item in occurrences_here}
                        complete_paths = bool(entry_keys and node_seqids) and all(
                            (entry_key, seqid) in self.embedded_sources
                            for entry_key in entry_keys for seqid in node_seqids)
                        sample_states[sample] = "Absent" if complete_paths else "NA"
                        continue
                    sample_te = [(x.node_start, x.node_end) for x in te_rows if x.sample == sample]
                    fraction = _covered_bp(sample_te) / node.length if node.length else 0.0
                    if fraction >= self.te_min_fraction:
                        sample_states[sample] = "TE"
                    elif sample_te:
                        sample_states[sample] = "PartialTE"
                    elif any(self._entry_complete(x) for x in sample_occ):
                        sample_states[sample] = "NonTE"
                    else:
                        sample_states[sample] = "Unknown"
                values = {
                    "node_id": node_id, "length": node.length,
                    "sequence_status": "available" if node.sequence is not None else "unavailable",
                    "identity": identity,
                    "is_te": "1" if te_rows else ("0" if callable_contexts else "NA"),
                    "te_covered_bp": te_covered, "te_fraction": f"{te_fraction:.6f}",
                    "te_occurrence_support": f"{te_occ_support:.6f}",
                    "te_genome_support": f"{te_genome_support:.6f}",
                    "annotated_occurrences": len(callable_contexts),
                    "total_occurrences": len(contexts),
                    "callable_genomes": len(callable_samples), "primary_te_family": primary_family,
                    "primary_te_class": primary_class,
                    "feature_types": ",".join(feature_types) or "None",
                    "gene_structure": ",".join(gene_types) or "None",
                    "conflict": int(node_id in conflict_nodes or disagreement),
                    "confidence": confidence,
                }
                values.update({f"sample:{sample}": sample_states[sample] for sample in samples})
                self.node_summaries[group] = (
                    identity, values["is_te"], values["te_fraction"],
                    primary_family, primary_class, values["gene_structure"], confidence,
                )
                writer.writerow(values)
                for sample in samples:
                    matrix_handle.write(
                        f"{_clean(node_id)}\t{_clean(sample)}\t{sample_states[sample]}\n")
        del occurrence_by_node, occurrence_node_offsets, projected_by_node, projected_node_offsets
        with open(self.conflict_path, "wt") as handle:
            handle.write("node_id\toccurrence_id\tsample\treason\toverlap_bp\n")
            for row in conflicts:
                handle.write("\t".join(map(str, row)) + "\n")
        with open(self.unresolved_path, "wt") as handle:
            handle.write("record\treason\tdetail\n")
            for row in self.unresolved:
                handle.write("\t".join(_clean(x) for x in row) + "\n")
        self._write_repeat_tables(sorted_projected, columns)
        del sorted_projected
        self._write_gff3(np.lexsort((
            np.arange(len(projected)), columns["node_end"], columns["node_start"], rank_rows)))
        del columns, rank_rows
        self._write_annotated_gfa()
        self._write_provenance(conflicts)

    def _write_repeat_tables(
            self, sorted_projected: np.ndarray, columns: Dict[str, np.ndarray]) -> None:
        """Write RepeatMasker-style graph-space and source-space summaries."""
        occurrences = self.occurrences
        source = _column(occurrences.source, np.uint32)
        genome_start = _column(occurrences.genome_start, np.int64)
        genome_end = _column(occurrences.genome_end, np.int64)
        by_source = np.argsort(source, kind="stable")
        source_offsets = _group_offsets(source[by_source], len(occurrences.sources.values)).tolist()
        callable_intervals: Dict[
            Tuple[Tuple[str, str], str], List[Tuple[int, int]]] = {}
        for code, key in enumerate(occurrences.sources.values):
            members = by_source[source_offsets[code]:source_offsets[code + 1]]
            if len(members):
                callable_intervals[key] = _merge_intervals(zip(
                    genome_start[members].tolist(), genome_end[members].tolist()))
        del by_source
        callable_ends = {
            key: [end for _start, end in intervals]
            for key, intervals in callable_intervals.items()
        }

        source_union = _StreamingUnion()
        source_element_counts: Counter = Counter()
        graph_element_counts: Counter = Counter()
        counted_repeatmasker_ids: Dict[Tuple[str, str], Set[str]] = defaultdict(set)
        for feature in sorted(
                self.features,
                key=lambda x: (x.entry_key, x.seqid, x.start, x.end,
                               x.feature_id, x.feature_type)):
            buckets = _repeat_summary_buckets(
                feature.category, feature.feature_type,
                feature.family, feature.te_class)
            source_key = (feature.entry_key, feature.seqid)
            intervals = callable_intervals.get(source_key)
            if not buckets or not intervals:
                continue
            index = bisect_right(callable_ends[source_key], feature.start)
            overlaps: List[Tuple[int, int]] = []
            while index < len(intervals) and intervals[index][0] < feature.end:
                start = max(feature.start, intervals[index][0])
                end = min(feature.end, intervals[index][1])
                if end > start:
                    overlaps.append((start, end))
                index += 1
            if not overlaps:
                continue
            count_element = True
            if feature.parent.startswith("RMID:"):
                seen_ids = counted_repeatmasker_ids[feature.entry_key]
                fragment_key = f"{feature.te_class}\0{feature.parent}"
                count_element = fragment_key not in seen_ids
                seen_ids.add(fragment_key)
            for bucket in buckets:
                if count_element:
                    source_element_counts[(feature.entry_key, bucket)] += 1
                    graph_element_counts[bucket] += 1
                for start, end in overlaps:
                    source_union.add(
                        (feature.entry_key, bucket, feature.seqid), start, end)

        source_lengths: Counter = Counter()
        for (entry_key, bucket, _seqid), length in source_union.total.items():
            source_lengths[(entry_key, bucket)] += length

        # Rows arrive grouped by node and sorted by node_start, so the per-node
        # interval union state can be dropped whenever the node changes.
        graph_lengths: Counter = Counter()
        bucket_cache: Dict[int, Tuple[str, ...]] = {}
        node_union = _StreamingUnion()
        current_node = None
        for start, end in self._chunks(len(sorted_projected)):
            order = sorted_projected[start:end]
            for node, node_start, node_end, feature_index in zip(
                    columns["node"][order].tolist(), columns["node_start"][order].tolist(),
                    columns["node_end"][order].tolist(), columns["feature"][order].tolist()):
                if node != current_node:
                    graph_lengths.update(node_union.total)
                    node_union = _StreamingUnion()
                    current_node = node
                buckets = bucket_cache.get(feature_index)
                if buckets is None:
                    feature = self.features[feature_index]
                    buckets = bucket_cache[feature_index] = _repeat_summary_buckets(
                        feature.category, feature.feature_type, feature.family,
                        feature.te_class)
                for bucket in buckets:
                    node_union.add(bucket, node_start, node_end)
        graph_lengths.update(node_union.total)

        graph_total = sum(node.length for node in self.nodes.values())
        complete_node_sequence = all(
            node.sequence is not None for node in self.nodes.values())
        if complete_node_sequence:
            graph_nx = sum(
                node.sequence.upper().count("N") + node.sequence.upper().count("X")
                for node in self.nodes.values())
            graph_gc_bp = sum(
                node.sequence.upper().count("G") + node.sequence.upper().count("C")
                for node in self.nodes.values())
            graph_non_nx = graph_total - graph_nx
            graph_gc = (100.0 * graph_gc_bp / graph_non_nx
                        if graph_non_nx else None)
        else:
            graph_non_nx = None
            graph_gc = None
        graph_stats = {
            "scope": "graph", "sample": ".", "haplotype": ".",
            "coordinate_space": "unique GFA node sequence space",
            "denominator_kind": (
                "complete_node_sequences" if complete_node_sequence
                else "node_lengths_with_incomplete_sequences"),
            "sequence_count": len(self.nodes), "total_length": graph_total,
            "non_nx_length": graph_non_nx, "gc_percent": graph_gc,
        }
        with open(self.repeat_table_path, "wt") as handle:
            handle.write(self._format_repeat_table_block(
                os.path.basename(self.gfa_file), graph_stats,
                graph_element_counts, graph_lengths, graph_space=True))

        sample_stats = []
        with open(self.repeat_sample_table_path, "wt") as handle:
            for entry in self.entries.values():
                source_keys = [key for key in callable_intervals
                               if key[0] == entry.key]
                if entry.fasta:
                    sequence_count, total, nx, gc = _fasta_sequence_stats(entry.fasta)
                    non_nx = total - nx
                    denominator_kind = "source_fasta"
                else:
                    sequence_count = len(source_keys)
                    total = sum(
                        end - start for key in source_keys
                        for start, end in callable_intervals[key])
                    non_nx = None
                    gc = None
                    denominator_kind = "graph_callable_source_intervals"
                stats = {
                    "scope": "sample", "sample": entry.sample,
                    "haplotype": entry.haplotype,
                    "coordinate_space": "source genome coordinates",
                    "denominator_kind": denominator_kind,
                    "sequence_count": sequence_count, "total_length": total,
                    "non_nx_length": non_nx, "gc_percent": gc,
                }
                sample_stats.append((entry, stats))
                counts, lengths = self._source_repeat_counters(
                    entry.key, source_element_counts, source_lengths)
                handle.write(self._format_repeat_table_block(
                    os.path.basename(entry.fasta) if entry.fasta else
                    f"{entry.sample}#{entry.haplotype}",
                    stats, counts, lengths, graph_space=False))
                handle.write("\n")

        with open(self.repeat_summary_path, "wt", newline="") as handle:
            fields = [
                "scope", "sample", "haplotype", "coordinate_space",
                "denominator_kind", "sequence_count", "total_length",
                "non_nx_length", "gc_percent", "category_key",
                "category_label", "element_count", "occupied_bp", "percentage",
            ]
            writer = csv.DictWriter(
                handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            scopes = [(graph_stats, graph_element_counts, graph_lengths)]
            for entry, stats in sample_stats:
                counts, lengths = self._source_repeat_counters(
                    entry.key, source_element_counts, source_lengths)
                scopes.append((stats, counts, lengths))
            for stats, counts, lengths in scopes:
                denominator = stats["total_length"]
                for category_key, category_label in _RM_ROWS + (("masked", "Bases masked"),):
                    occupied = lengths[category_key]
                    writer.writerow({
                        **stats,
                        "non_nx_length": (
                            stats["non_nx_length"]
                            if stats["non_nx_length"] is not None else "NA"),
                        "gc_percent": (
                            f"{stats['gc_percent']:.6f}"
                            if stats["gc_percent"] is not None else "NA"),
                        "category_key": category_key,
                        "category_label": category_label.strip().rstrip(":"),
                        "element_count": (
                            "NA" if category_key in {"interspersed", "masked"}
                            else counts[category_key]),
                        "occupied_bp": occupied,
                        "percentage": (
                            f"{100.0 * occupied / denominator:.6f}"
                            if denominator else "0.000000"),
                    })

    @staticmethod
    def _source_repeat_counters(
            entry_key: Tuple[str, str], element_counts: Counter,
            occupied_lengths: Counter) -> Tuple[Counter, Counter]:
        keys = [key for key, _label in _RM_ROWS] + ["masked"]
        return (
            Counter({key: element_counts[(entry_key, key)] for key in keys}),
            Counter({key: occupied_lengths[(entry_key, key)] for key in keys}),
        )

    @staticmethod
    def _format_repeat_table_block(
            file_name: str, stats: Dict[str, object], counts: Counter,
            lengths: Counter, graph_space: bool) -> str:
        total = int(stats["total_length"])
        non_nx = stats["non_nx_length"]
        masked = lengths["masked"]
        masked_percent = 100.0 * masked / total if total else 0.0
        lines = [
            "=" * 50, f"file name: {file_name}",
            f"scope: {stats['coordinate_space']}",
            f"denominator: {stats['denominator_kind']}",
            f"sequences: {int(stats['sequence_count']):13d}",
        ]
        if non_nx is None:
            lines.extend((
                f"total length: {total:13d} bp (N/X-excluded length unavailable)",
                "GC level: NA",
            ))
        else:
            lines.append(
                f"total length: {total:13d} bp ({int(non_nx)} bp excl N/X-runs)")
            gc = stats["gc_percent"]
            lines.append(f"GC level: {gc:.2f} %" if gc is not None else "GC level: NA")
        lines.extend((
            f"bases masked: {masked:13d} bp ({masked_percent:6.2f} %)",
            "=" * 50,
            "                       number of      length       percentage",
            "                       elements*    occupied     of sequence",
            "------------------------------------------------------------",
        ))
        for key, label in _RM_ROWS:
            occupied = lengths[key]
            percentage = 100.0 * occupied / total if total else 0.0
            count = "" if key == "interspersed" else str(counts[key])
            lines.append(
                f"{label:<40}{count:>10} {occupied:>12} bp {percentage:>7.2f} %")
            if key in {"retroelements", "ltr", "dna", "interspersed"}:
                lines.append("")
        lines.extend((
            "* RepeatMasker fragments sharing class and ID count as one element;",
            "  other inputs count distinct source records overlapping callable sequence.",
            "  Occupied lengths are coordinate-space interval unions.",
        ))
        if graph_space:
            lines.extend((
                "  Shared nodes are counted once. Class rows are independent unions and",
                "  can overlap when assemblies assign different identities to one node.",
            ))
        else:
            lines.append(
                "  Source annotations are intersected with graph-callable genome intervals.")
        lines.extend((
            "  Total interspersed repeats unions retroelement, DNA and unclassified",
            "  aggregate intervals; Rolling-circles remains a separate RepeatMasker row.",
            "Generated by Pan_TE GFA annotation; this is a RepeatMasker-compatible",
            "summary, not a RepeatMasker run.",
        ))
        return "\n".join(lines) + "\n"

    def _write_gff3(self, order: np.ndarray) -> None:
        """Write projections in a stable (node_id, node_start, node_end) order."""
        with open(self.gff3_path, "wt") as handle:
            handle.write("##gff-version 3\n")
            rows = (self.projected[i] for start, end in self._chunks(len(order))
                    for i in order[start:end].tolist())
            for index, row in enumerate(rows, 1):
                def esc(value: object) -> str:
                    return quote(_clean(value), safe="._:-")
                attrs = ";".join((
                    f"ID=proj{index}", f"Occurrence={esc(row.occurrence_id)}",
                    f"Sample={esc(row.sample)}", f"Haplotype={esc(row.haplotype)}",
                    f"GenomeSeqid={esc(row.seqid)}", f"GenomeStart={row.genome_start}",
                    f"GenomeEnd={row.genome_end}", f"Orientation={row.orientation}",
                    f"OriginalID={esc(row.feature_id)}", f"OriginalParent={esc(row.parent)}",
                    f"Category={esc(row.category)}", f"TEFamily={esc(row.family)}",
                    f"TEClass={esc(row.te_class)}",
                    f"OriginalAttributes={esc(row.attributes)}",
                ))
                handle.write("\t".join(map(str, (
                    row.node_id, "Pan_TE_projection", row.feature_type,
                    row.node_start + 1, row.node_end, ".", row.strand, row.phase, attrs,
                ))) + "\n")

    def _write_annotated_gfa(self) -> None:
        reserved = {"ai", "it", "tf", "tc", "fm", "gs", "cf"}
        with open(self.gfa_file, "rt") as source, open(self.annotated_gfa_path, "wt") as target:
            for line in source:
                if not line.startswith("S\t"):
                    target.write(line)
                    continue
                fields = line.rstrip("\n").split("\t")
                node = self.nodes.get(fields[1])
                summary = (self.node_summaries[node.index]
                           if node is not None and node.index < len(self.node_summaries)
                           else None)
                if not summary:
                    target.write(line)
                    continue
                fields = fields[:3] + [x for x in fields[3:] if x[:2] not in reserved or len(x) < 3 or x[2] != ":"]
                (identity, is_te, te_fraction, te_family, te_class,
                 gene_structure, confidence) = summary
                fields.extend([
                    f"ai:Z:{identity}", f"it:Z:{is_te}",
                    f"tf:f:{te_fraction}",
                    f"tc:Z:{_clean(te_class)}",
                    f"fm:Z:{_clean(te_family)}",
                    f"gs:Z:{_clean(gene_structure)}", f"cf:Z:{confidence}",
                ])
                target.write("\t".join(fields) + "\n")

    def _write_provenance(self, conflicts: Sequence[Tuple]) -> None:
        payload = {
            "contract": "pan-te-gfa-all-genome-annotation-v1",
            "coordinate_system": "internal_and_tsv_0_based_half_open;gff3_1_based_inclusive",
            "graph_mode": self.graph_mode,
            "gfa": {"path": self.gfa_file, "sha256": _sha256(self.gfa_file)},
            "manifest": {"path": self.annotation_manifest, "sha256": _sha256(self.annotation_manifest)},
            "global_gaf": ({"path": self.default_gaf, "sha256": _sha256(self.default_gaf)}
                           if self.default_gaf else None),
            "thresholds": {"te_min_fraction": self.te_min_fraction,
                           "partial_min_fraction": self.partial_min_fraction,
                           "gaf_min_mapq": self.min_mapq,
                           "gaf_min_identity": self.min_identity},
            "counts": {"manifest_entries": len(self.entries), "features": len(self.features),
                       "nodes": len(self.nodes), "occurrences": len(self.occurrences),
                       "projected_intervals": len(self.projected), "conflicts": len(conflicts),
                       "unresolved": len(self.unresolved)},
            "gaf_statistics": dict(sorted(self.gaf_stats.items())),
            "outputs": {
                "graph_repeatmasker_style_table": self.repeat_table_path,
                "sample_repeatmasker_style_table": self.repeat_sample_table_path,
                "machine_readable_repeat_summary": self.repeat_summary_path,
            },
            "inputs": [
                {"sample": entry.sample, "haplotype": entry.haplotype,
                 "annotation_complete": entry.annotation_complete,
                 "fasta": entry.fasta, "gff3": entry.gff3,
                 "te_annotation": entry.te_annotation, "seqid_map": entry.seqid_map,
                 "gaf": entry.gaf,
                 "checksums": {name: _sha256(path) for name, path in (
                     ("fasta", entry.fasta), ("gff3", entry.gff3),
                     ("te_annotation", entry.te_annotation),
                     ("seqid_map", entry.seqid_map), ("gaf", entry.gaf)) if path}}
                for entry in self.entries.values()
            ],
            "semantics": {
                "occurrence_table_is_authoritative": True,
                "missing_annotation_is_not_negative_without_annotation_complete": True,
                "rgfa_without_embedded_paths_never_infers_absence": True,
                "graph_tbl_uses_unique_node_coordinate_unions": True,
                "sample_tbl_uses_source_coordinate_unions_intersected_with_graph_callable_intervals": True,
                "tbl_repeatmasker_elements_are_distinct_class_and_fragment_id": True,
                "tbl_other_elements_are_distinct_source_annotation_records": True,
            },
        }
        with open(self.provenance_path, "wt") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")


__all__ = [
    "GfaAnnotator", "Feature", "ManifestEntry", "Node", "Occurrence",
    "ProjectedFeature", "_covered_bp", "_merge_intervals", "_simple_overlap",
]


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Project per-genome TE/GFF3 annotations onto GFA nodes")
    parser.add_argument("--gfa", required=True, help="PGGB GFA or minigraph rGFA")
    parser.add_argument("--manifest", required=True,
                        help="Tab-delimited all-genome annotation manifest")
    parser.add_argument("--out-dir", required=True, help="Annotation output directory")
    parser.add_argument("--gaf", help="Optional GAF for rGFA sample occurrences")
    parser.add_argument("--gaf-min-mapq", type=int, default=20)
    parser.add_argument("--gaf-min-identity", type=float, default=0.90)
    parser.add_argument("--te-min-fraction", type=float, default=0.80)
    parser.add_argument("--partial-min-fraction", type=float, default=0.10)
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="[%(asctime)s] %(levelname)s: %(message)s")
    annotator = GfaAnnotator(
        gfa_file=args.gfa, output_dir=args.out_dir,
        annotation_manifest=args.manifest, gfa_gaf=args.gaf,
        gfa_gaf_min_mapq=args.gaf_min_mapq,
        gfa_gaf_min_identity=args.gaf_min_identity,
        te_min_fraction=args.te_min_fraction,
        partial_min_fraction=args.partial_min_fraction,
    )
    annotator.annotate()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
