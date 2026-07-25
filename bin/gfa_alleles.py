#!/usr/bin/env python3
"""Locus-aware, streaming PGGB allele extraction for Pan_TE.

Ordered path occurrences, rather than global node membership, define alternate
loci.  The implementation preserves copy support, tiles long alleles without
loss, and can scan alternate P lines one at a time for cohort-scale graphs.
"""

from array import array
from bisect import bisect_left
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass, field
import csv
import hashlib
import json
import os
import sqlite3
import sys
from typing import (Callable, Dict, Iterable, Iterator, List, Optional,
                    Sequence, Tuple)


Step = Tuple[str, str]


def path_group(path_name: str) -> str:
    parts = path_name.split('#')
    return '#'.join(parts[:2]) if len(parts) >= 3 else path_name


def path_contig(path_name: str) -> str:
    parts = path_name.split('#')
    return parts[-1] if len(parts) >= 3 else path_name


def reverse_complement(seq: str) -> str:
    complement = str.maketrans(
        'ACGTacgtRYSWKMBDHVNryswkmbdhvn',
        'TGCAtgcaYRSWMKVHDBNyrswmkvhdbn')
    return seq.translate(complement)[::-1]


def build_sequence(steps: Sequence[Step], sequences: Dict[str, str]) -> str:
    parts = []
    for node_id, orientation in steps:
        if node_id not in sequences:
            raise ValueError(f'path references missing segment {node_id!r}')
        sequence = sequences[node_id]
        if not sequence or sequence == '*':
            raise ValueError(f'segment {node_id!r} has no embedded sequence')
        if orientation == '-':
            sequence = reverse_complement(sequence)
        parts.append(sequence)
    return ''.join(parts)


def parse_p_steps(fields: Sequence[str]) -> List[Step]:
    """Parse a P record and reject unsupported non-zero path overlaps."""
    if len(fields) < 3 or not fields[2]:
        raise ValueError('malformed P record: missing path steps')
    steps = []
    for token in fields[2].split(','):
        if len(token) < 2 or token[-1] not in '+-':
            raise ValueError(f'malformed PGGB path step {token!r}')
        steps.append((sys.intern(token[:-1]), token[-1]))
    if len(fields) >= 4 and fields[3] not in ('', '*'):
        overlaps = fields[3].split(',')
        if len(overlaps) != max(0, len(steps) - 1):
            raise ValueError('P overlap count does not match path step count')
        unsupported = [overlap for overlap in overlaps if overlap != '0M']
        if unsupported:
            raise ValueError(
                'locus-aware PGGB mode requires 0M or * path overlaps; '
                f'found {unsupported[0]!r}')
    return steps


def window_starts(length: int, window: int, overlap: int,
                  min_size: int = 1) -> List[int]:
    """Cover a sequence without emitting a terminal window below min_size."""
    if length <= window:
        return [0]
    step = window - overlap
    if step <= 0:
        raise ValueError('window overlap must be smaller than the window size')
    starts = []
    start = 0
    while start + window < length:
        starts.append(start)
        start += step
    if length - start < min_size:
        start = max(0, length - window)
    if not starts or starts[-1] != start:
        starts.append(start)
    return starts


@dataclass
class AlleleCandidate:
    candidate_id: str
    parent_allele_id: str
    sequence: str = field(repr=False)
    sample: str = ''
    path: str = ''
    contig: str = ''
    ref_path: str = ''
    path_step_start: int = 0
    path_step_end: int = 0
    left_anchor: str = ''
    right_anchor: str = ''
    left_ref_rank: Optional[int] = None
    right_ref_rank: Optional[int] = None
    ref_span_bp: Optional[int] = None
    alt_span_bp: int = 0
    net_change_bp: Optional[int] = None
    allele_type: str = ''
    chunk_index: int = 1
    chunk_start: int = 0
    chunk_end: int = 0
    locus_key: str = ''
    sequence_sha256: str = ''
    representative_id: str = ''
    occurrence_count: int = 1
    support_sample_count: int = 1
    support_samples: str = ''
    qc_status: str = 'pending'
    # rGFA+GAF provenance. PGGB records retain these defaults, while the GAF
    # backend records the query alignment from which an occurrence was called.
    sequence_source: str = 'graph_path'
    source_query: str = ''
    query_start: Optional[int] = None
    query_end: Optional[int] = None
    gaf_record_id: Optional[int] = None
    mapq: Optional[int] = None
    alignment_identity: Optional[float] = None

    def fasta_header(self) -> str:
        ref_bp = self.ref_span_bp if self.ref_span_bp is not None else 'NA'
        return (
            f'>{self.candidate_id} type={self.allele_type} '
            f'support={self.support_sample_count} '
            f'occurrences={self.occurrence_count} '
            f'alt_bp={self.alt_span_bp} ref_bp={ref_bp}')

    def manifest_row(self) -> dict:
        row = asdict(self)
        row.pop('sequence')
        return row


MANIFEST_FIELDS = [
    'candidate_id', 'parent_allele_id', 'sample', 'path', 'contig',
    'ref_path', 'path_step_start', 'path_step_end', 'left_anchor',
    'right_anchor', 'left_ref_rank', 'right_ref_rank', 'ref_span_bp',
    'alt_span_bp', 'net_change_bp', 'allele_type', 'chunk_index',
    'chunk_start', 'chunk_end', 'locus_key', 'sequence_sha256',
    'representative_id', 'occurrence_count', 'support_sample_count',
    'support_samples', 'qc_status', 'sequence_source', 'source_query',
    'query_start', 'query_end', 'gaf_record_id', 'mapq',
    'alignment_identity',
]


class SQLiteAlleleStore:
    """Disk-backed exact locus aggregation for cohort-scale path streaming.

    Only one sequence per exact ``(locus, sha256)`` allele is retained.  Every
    occurrence remains in the temporary database as compact JSON metadata, so
    the final manifest preserves copy and sample provenance without keeping
    cohort-wide occurrence sequences in RAM.
    """

    def __init__(self, path: str):
        self.path = path
        if os.path.exists(path):
            os.remove(path)
        self.connection = sqlite3.connect(path)
        self.connection.execute('PRAGMA journal_mode=OFF')
        self.connection.execute('PRAGMA synchronous=OFF')
        self.connection.executescript("""
            CREATE TABLE representatives (
                representative_id TEXT PRIMARY KEY,
                locus_key TEXT NOT NULL,
                sequence_sha256 TEXT NOT NULL,
                sequence TEXT NOT NULL,
                payload TEXT NOT NULL,
                occurrence_count INTEGER NOT NULL DEFAULT 1,
                support_sample_count INTEGER NOT NULL DEFAULT 1,
                support_samples TEXT NOT NULL DEFAULT '',
                qc_status TEXT NOT NULL DEFAULT 'pending',
                UNIQUE(locus_key, sequence_sha256)
            );
            CREATE TABLE occurrences (
                ordinal INTEGER PRIMARY KEY AUTOINCREMENT,
                payload TEXT NOT NULL,
                representative_id TEXT NOT NULL
            );
            CREATE TABLE supports (
                representative_id TEXT NOT NULL,
                sample TEXT NOT NULL,
                PRIMARY KEY(representative_id, sample)
            );
        """)
        self.occurrence_count = 0
        self.representative_count = 0

    def add(self, occurrence: AlleleCandidate) -> None:
        key = (occurrence.locus_key, occurrence.sequence_sha256)
        row = self.connection.execute(
            'SELECT representative_id FROM representatives '
            'WHERE locus_key=? AND sequence_sha256=?', key).fetchone()
        payload = json.dumps(occurrence.manifest_row(), separators=(',', ':'))
        if row is None:
            representative_id = occurrence.candidate_id
            self.connection.execute(
                'INSERT INTO representatives '
                '(representative_id,locus_key,sequence_sha256,sequence,payload) '
                'VALUES (?,?,?,?,?)',
                (representative_id, occurrence.locus_key,
                 occurrence.sequence_sha256, occurrence.sequence, payload))
            self.representative_count += 1
        else:
            representative_id = row[0]
            self.connection.execute(
                'UPDATE representatives SET occurrence_count=occurrence_count+1 '
                'WHERE representative_id=?', (representative_id,))
        self.connection.execute(
            'INSERT OR IGNORE INTO supports (representative_id,sample) '
            'VALUES (?,?)', (representative_id, occurrence.sample))
        self.connection.execute(
            'INSERT INTO occurrences (payload,representative_id) VALUES (?,?)',
            (payload, representative_id))
        self.occurrence_count += 1

    def finalize(self) -> None:
        cursor = self.connection.execute(
            "SELECT representative_id,group_concat(sample, ','),count(*) "
            'FROM (SELECT representative_id,sample FROM supports '
            'ORDER BY representative_id,sample) GROUP BY representative_id')
        while True:
            rows = cursor.fetchmany(10000)
            if not rows:
                break
            self.connection.executemany(
                'UPDATE representatives SET support_samples=?, '
                'support_sample_count=? WHERE representative_id=?',
                ((samples, count, representative_id)
                 for representative_id, samples, count in rows))
        self.connection.commit()

    def representatives(self) -> List[AlleleCandidate]:
        records = []
        query = (
            'SELECT payload,sequence,representative_id,occurrence_count,'
            'support_sample_count,support_samples,qc_status '
            'FROM representatives ORDER BY rowid')
        for (payload, sequence, representative_id, occurrence_count,
             support_count, support_samples, qc_status) in self.connection.execute(query):
            values = json.loads(payload)
            values.update(
                representative_id=representative_id,
                occurrence_count=occurrence_count,
                support_sample_count=support_count,
                support_samples=support_samples,
                qc_status=qc_status)
            records.append(AlleleCandidate(sequence=sequence, **values))
        return records

    def set_qc(self, passed_ids: Iterable[str]) -> None:
        self.connection.execute(
            "UPDATE representatives SET qc_status='filtered'")
        self.connection.executemany(
            "UPDATE representatives SET qc_status='passed' "
            'WHERE representative_id=?', ((value,) for value in passed_ids))
        self.connection.commit()

    def write_manifest(self, path: str) -> None:
        query = (
            'SELECT o.payload,o.representative_id,r.occurrence_count,'
            'r.support_sample_count,r.support_samples,r.qc_status '
            'FROM occurrences o JOIN representatives r '
            'ON o.representative_id=r.representative_id ORDER BY o.ordinal')
        with open(path, 'w', newline='') as handle:
            writer = csv.DictWriter(
                handle, fieldnames=MANIFEST_FIELDS, delimiter='\t',
                extrasaction='ignore')
            writer.writeheader()
            for (payload, representative_id, occurrence_count,
                 support_count, support_samples, qc_status) in self.connection.execute(query):
                row = json.loads(payload)
                row.update(
                    representative_id=representative_id,
                    occurrence_count=occurrence_count,
                    support_sample_count=support_count,
                    support_samples=support_samples,
                    qc_status=qc_status)
                writer.writerow(row)

    def close(self, remove: bool = False) -> None:
        if self.connection is not None:
            self.connection.close()
            self.connection = None
        if remove and os.path.exists(self.path):
            os.remove(self.path)


@dataclass
class ReferenceModel:
    name: str
    contig: str
    steps: Sequence[Step] = field(repr=False)
    occurrences: dict = field(repr=False)
    prefix: array = field(repr=False)
    unique_keys: set = field(repr=False)
    highly_repetitive_nodes: int = 0


def _classify_allele(ref_span: Optional[int], alt_span: int,
                     min_change: int) -> Tuple[str, Optional[int]]:
    if ref_span is None:
        return 'unresolved', None
    delta = alt_span - ref_span
    if ref_span == 0:
        return 'insertion', delta
    if delta >= min_change:
        return 'insertion_like', delta
    if delta <= -min_change:
        return 'deletion_associated', delta
    return 'replacement', delta


def _reference_model(name: str, steps: Sequence[Step],
                     sequences: Dict[str, str], max_anchor_occurrences: int
                     ) -> ReferenceModel:
    occurrences = defaultdict(lambda: array('I'))
    prefix = array('Q', [0])
    total = 0
    for rank, (node_id, orientation) in enumerate(steps):
        if node_id not in sequences:
            raise ValueError(
                f'reference path {name!r} references missing segment {node_id!r}')
        sequence = sequences[node_id]
        if not sequence or sequence == '*':
            raise ValueError(f'segment {node_id!r} has no embedded sequence')
        total += len(sequence)
        prefix.append(total)
        occurrences[(node_id, orientation)].append(rank)
    highly_repetitive = 0
    for key in list(occurrences):
        if len(occurrences[key]) > max_anchor_occurrences:
            del occurrences[key]
            highly_repetitive += 1
    unique_keys = {
        key for key, ranks in occurrences.items() if len(ranks) == 1}
    return ReferenceModel(
        name, path_contig(name), steps, dict(occurrences), prefix,
        unique_keys, highly_repetitive)


def _maximum_collinear_anchors(
        steps: Sequence[Step], reference: ReferenceModel
        ) -> List[Tuple[int, int, str, str]]:
    """Return a maximum collinear chain using occurrence-aware LCS chaining.

    Reference nodes with a bounded number of occurrences contribute every
    possible occurrence.  Hunt-Szymanski chaining resolves repeated nodes and
    avoids the first-fit failure mode in rearranged paths.
    """
    tail_ranks = array('I')
    tail_trace = array('q')
    trace_path_rank = array('I')
    trace_ref_rank = array('I')
    trace_previous = array('q')

    for path_rank, step in enumerate(steps):
        ranks = reference.occurrences.get(step)
        if not ranks:
            continue
        # Reverse order prevents one alternate step from extending itself more
        # than once when it has multiple matching reference occurrences.
        for ref_rank in reversed(ranks):
            chain_rank = bisect_left(tail_ranks, ref_rank)
            previous = tail_trace[chain_rank - 1] if chain_rank else -1
            trace_index = len(trace_path_rank)
            trace_path_rank.append(path_rank)
            trace_ref_rank.append(ref_rank)
            trace_previous.append(previous)
            if chain_rank == len(tail_ranks):
                tail_ranks.append(ref_rank)
                tail_trace.append(trace_index)
            else:
                tail_ranks[chain_rank] = ref_rank
                tail_trace[chain_rank] = trace_index

    if not tail_trace:
        return []
    result = []
    trace_index = tail_trace[-1]
    while trace_index >= 0:
        path_rank = trace_path_rank[trace_index]
        ref_rank = trace_ref_rank[trace_index]
        node_id, orientation = steps[path_rank]
        result.append((path_rank, ref_rank, node_id, orientation))
        trace_index = trace_previous[trace_index]
    result.reverse()
    # Unique anchors provide stable locus boundaries and avoid splitting one
    # alternate allele around a low-copy repeat.  Repeated occurrences remain
    # available as a fallback in regions lacking two unique anchors.
    unique_result = [
        anchor for anchor in result
        if (anchor[2], anchor[3]) in reference.unique_keys]
    return unique_result if len(unique_result) >= 2 else result


def _resolve_reference_models(
        ref_paths: Dict[str, Sequence[Step]], sequences: Dict[str, str],
        max_anchor_occurrences: int, counters: Counter):
    exact = {}
    folded = defaultdict(list)
    models = {}
    for name, steps in ref_paths.items():
        contig = path_contig(name)
        if contig in exact:
            raise ValueError(
                f'multiple selected reference paths have contig {contig!r}')
        model = _reference_model(
            name, steps, sequences, max_anchor_occurrences)
        exact[contig] = model
        folded[contig.casefold()].append(model)
        models[name] = model
        counters['reference_anchor_keys'] += len(model.occurrences)
        counters['highly_repetitive_anchor_keys_excluded'] += (
            model.highly_repetitive_nodes)
    return exact, folded, models


def _match_reference(contig: str, exact: dict, folded: dict
                     ) -> Optional[ReferenceModel]:
    if contig in exact:
        return exact[contig]
    matches = folded.get(contig.casefold(), [])
    return matches[0] if len(matches) == 1 else None


def _extract_from_records(
        alternate_records: Iterable[Tuple[str, Sequence[Step]]],
        ref_paths: Dict[str, Sequence[Step]], sequences: Dict[str, str],
        min_size: int, window_size: int, window_overlap: int,
        max_anchor_occurrences: int, stats: Optional[dict],
        candidate_sink: Optional[Callable[[AlleleCandidate], None]] = None
        ) -> List[AlleleCandidate]:
    if window_size < min_size:
        raise ValueError('window_size must be at least min_size')
    if window_overlap < 0 or window_overlap >= window_size:
        raise ValueError('window_overlap must satisfy 0 <= overlap < window_size')

    counters = Counter()
    exact, folded, _ = _resolve_reference_models(
        ref_paths, sequences, max_anchor_occurrences, counters)
    candidates = []
    candidate_number = 0
    allele_number = 0
    unmatched = []

    def emit_allele(
            alt_steps: Sequence[Step], ref_steps: Sequence[Step], *,
            sample: str, path_name: str, contig: str,
            reference: ReferenceModel, path_start: int, path_end: int,
            locus_key: str, left_anchor: str = '', right_anchor: str = '',
            left_ref_rank: Optional[int] = None,
            right_ref_rank: Optional[int] = None,
            ref_span: Optional[int] = None,
            allele_type_override: Optional[str] = None) -> None:
        nonlocal candidate_number, allele_number
        if alt_steps == ref_steps:
            counters['unchanged_alleles_skipped'] += 1
            return
        alt_sequence = build_sequence(alt_steps, sequences)
        ref_sequence = build_sequence(ref_steps, sequences)
        if alt_sequence.upper() == ref_sequence.upper():
            counters['unchanged_sequence_alleles_skipped'] += 1
            return
        alt_span = len(alt_sequence)
        counters['allele_occurrences_total'] += 1
        counters['allele_bp_total'] += alt_span
        if alt_span == 0:
            counters['alleles_empty_alt'] += 1
            return
        if alt_span < min_size:
            counters['alleles_below_min'] += 1
            counters['alleles_below_min_bp'] += alt_span
            return

        allele_type, delta = _classify_allele(
            ref_span, alt_span, min_size)
        if allele_type_override:
            allele_type = allele_type_override
        counters[f'allele_type_{allele_type}'] += 1
        allele_number += 1
        parent_id = f'gfa_allele_{allele_number:08d}'
        starts = window_starts(
            alt_span, window_size, window_overlap, min_size=min_size)
        if len(starts) > 1:
            counters['long_alleles_windowed'] += 1
            counters['long_allele_bp'] += alt_span

        for chunk_index, start in enumerate(starts, 1):
            end = min(start + window_size, alt_span)
            sequence = alt_sequence[start:end]
            candidate_number += 1
            candidate_id = f'gfa_candidate_{candidate_number:08d}'
            digest = hashlib.sha256(sequence.upper().encode()).hexdigest()
            candidate = AlleleCandidate(
                candidate_id=candidate_id,
                parent_allele_id=parent_id,
                sequence=sequence,
                sample=sample,
                path=path_name,
                contig=contig,
                ref_path=reference.name,
                path_step_start=path_start,
                path_step_end=path_end,
                left_anchor=left_anchor,
                right_anchor=right_anchor,
                left_ref_rank=left_ref_rank,
                right_ref_rank=right_ref_rank,
                ref_span_bp=ref_span,
                alt_span_bp=alt_span,
                net_change_bp=delta,
                allele_type=allele_type,
                chunk_index=chunk_index,
                chunk_start=start,
                chunk_end=end,
                locus_key=f'{locus_key}:chunk:{chunk_index}',
                sequence_sha256=digest,
                support_samples=sample,
            )
            if candidate_sink is None:
                candidates.append(candidate)
            else:
                candidate_sink(candidate)
            counters['candidate_occurrences'] += 1
            counters['candidate_bp'] += len(sequence)

    for path_name, steps in alternate_records:
        contig = path_contig(path_name)
        reference = _match_reference(contig, exact, folded)
        if reference is None:
            unmatched.append((path_name, contig))
            continue
        if contig != reference.contig:
            counters['casefold_contig_matches'] += 1
        sample = path_group(path_name)
        anchors = _maximum_collinear_anchors(steps, reference)
        counters['alternate_paths_processed'] += 1
        counters['collinear_anchors'] += len(anchors)

        if not anchors:
            emit_allele(
                steps, reference.steps, sample=sample, path_name=path_name,
                contig=contig, reference=reference,
                path_start=0, path_end=len(steps),
                locus_key=f'{reference.name}:UNRESOLVED_FULL_PATH',
                ref_span=int(reference.prefix[-1]),
                allele_type_override='unresolved_no_anchor')
            counters['paths_without_anchors'] += 1
            continue

        first_path_rank, first_ref_rank, first_node, _ = anchors[0]
        if first_path_rank > 0 or first_ref_rank > 0:
            emit_allele(
                steps[:first_path_rank], reference.steps[:first_ref_rank],
                sample=sample, path_name=path_name, contig=contig,
                reference=reference, path_start=0,
                path_end=first_path_rank,
                locus_key=f'{reference.name}:START:{first_ref_rank}',
                right_anchor=first_node,
                right_ref_rank=first_ref_rank,
                ref_span=int(reference.prefix[first_ref_rank]))

        for left, right in zip(anchors, anchors[1:]):
            left_path_rank, left_ref_rank, left_node, _ = left
            right_path_rank, right_ref_rank, right_node, _ = right
            alt_internal = steps[left_path_rank + 1:right_path_rank]
            ref_internal = reference.steps[left_ref_rank + 1:right_ref_rank]
            ref_span = int(
                reference.prefix[right_ref_rank]
                - reference.prefix[left_ref_rank + 1])
            emit_allele(
                alt_internal, ref_internal, sample=sample,
                path_name=path_name, contig=contig, reference=reference,
                path_start=left_path_rank + 1,
                path_end=right_path_rank,
                locus_key=(
                    f'{reference.name}:{left_ref_rank}:{right_ref_rank}'),
                left_anchor=left_node, right_anchor=right_node,
                left_ref_rank=left_ref_rank,
                right_ref_rank=right_ref_rank,
                ref_span=ref_span)

        last_path_rank, last_ref_rank, last_node, _ = anchors[-1]
        if (last_path_rank + 1 < len(steps)
                or last_ref_rank + 1 < len(reference.steps)):
            emit_allele(
                steps[last_path_rank + 1:],
                reference.steps[last_ref_rank + 1:],
                sample=sample, path_name=path_name, contig=contig,
                reference=reference, path_start=last_path_rank + 1,
                path_end=len(steps),
                locus_key=f'{reference.name}:{last_ref_rank}:END',
                left_anchor=last_node, left_ref_rank=last_ref_rank,
                ref_span=int(
                    reference.prefix[-1]
                    - reference.prefix[last_ref_rank + 1]))

    if unmatched:
        counters['unmatched_paths'] = len(unmatched)
        examples = ', '.join(
            f'{name} (contig={contig})' for name, contig in unmatched[:5])
        raise ValueError(
            f'{len(unmatched)} alternate path(s) have no matching reference '
            f'contig; examples: {examples}')
    if stats is not None:
        stats.update(counters)
        stats['candidate_occurrences'] = counters['candidate_occurrences']
    return candidates


def extract_pggb_alleles(
        paths: Dict[str, List[Step]], ref_path_names: Sequence[str],
        sequences: Dict[str, str], min_size: int = 100,
        window_size: int = 30000, window_overlap: int = 1000,
        max_anchor_occurrences: int = 64,
        stats: Optional[dict] = None) -> List[AlleleCandidate]:
    ref_set = set(ref_path_names)
    ref_paths = {name: paths[name] for name in ref_path_names}
    alternates = (
        (name, steps) for name, steps in paths.items() if name not in ref_set)
    return _extract_from_records(
        alternates, ref_paths, sequences, min_size, window_size,
        window_overlap, max_anchor_occurrences, stats)


def extract_pggb_alleles_stream(
        gfa_file: str, ref_paths: Dict[str, Sequence[Step]],
        sequences: Dict[str, str], min_size: int = 100,
        window_size: int = 30000, window_overlap: int = 1000,
        max_anchor_occurrences: int = 64,
        stats: Optional[dict] = None,
        candidate_sink: Optional[Callable[[AlleleCandidate], None]] = None
        ) -> List[AlleleCandidate]:
    """Scan alternate P records and optionally stream candidates to a sink."""
    ref_set = set(ref_paths)

    def records() -> Iterator[Tuple[str, Sequence[Step]]]:
        with open(gfa_file) as handle:
            for line_number, line in enumerate(handle, 1):
                if not line.startswith('P\t'):
                    continue
                fields = line.rstrip('\n').split('\t')
                name = fields[1]
                steps = parse_p_steps(fields)
                if name not in ref_set:
                    yield name, steps

    return _extract_from_records(
        records(), ref_paths, sequences, min_size, window_size,
        window_overlap, max_anchor_occurrences, stats, candidate_sink)


def deduplicate_locus_alleles(
        occurrences: Sequence[AlleleCandidate], stats: Optional[dict] = None
        ) -> List[AlleleCandidate]:
    """Collapse exact sequence duplicates only within the same graph locus."""
    representatives = []
    by_key = {}
    samples = {}
    counts = Counter()
    for occurrence in occurrences:
        key = (occurrence.locus_key, occurrence.sequence_sha256)
        representative = by_key.get(key)
        if representative is None:
            representative = occurrence
            by_key[key] = representative
            representatives.append(representative)
            samples[representative.candidate_id] = {occurrence.sample}
        else:
            samples[representative.candidate_id].add(occurrence.sample)
            occurrence.sequence = ''
        occurrence.representative_id = representative.candidate_id
        counts[representative.candidate_id] += 1

    representative_by_id = {record.candidate_id: record
                            for record in representatives}
    for representative in representatives:
        support = sorted(samples[representative.candidate_id])
        representative.occurrence_count = counts[representative.candidate_id]
        representative.support_sample_count = len(support)
        representative.support_samples = ','.join(support)
    for occurrence in occurrences:
        representative = representative_by_id[occurrence.representative_id]
        occurrence.occurrence_count = representative.occurrence_count
        occurrence.support_sample_count = representative.support_sample_count
        occurrence.support_samples = representative.support_samples

    if stats is not None:
        stats['candidate_occurrences'] = len(occurrences)
        stats['unique_locus_alleles'] = len(representatives)
        stats['exact_locus_duplicates_collapsed'] = (
            len(occurrences) - len(representatives))
    return representatives


def write_manifest(path: str, occurrences: Iterable[AlleleCandidate]) -> None:
    with open(path, 'w', newline='') as handle:
        writer = csv.DictWriter(
            handle, fieldnames=MANIFEST_FIELDS, delimiter='\t',
            extrasaction='ignore')
        writer.writeheader()
        for occurrence in occurrences:
            writer.writerow(occurrence.manifest_row())
