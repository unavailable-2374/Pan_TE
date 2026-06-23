#!/usr/bin/env python3
"""
PAF TE Enrichment Preprocessor for Pan_TE

Extracts TE insertion sequences from all-vs-all genome alignment PAF files.
Asymmetric gaps between syntenic alignment chains indicate TE insertions.
Extracted sequences are concatenated with N-spacers into synthetic scaffolds
and appended to the reference genome, enriching the input for the three
detection paths (LTR, mdl-repeat, RepGraph) so the final TE library covers
all haplotypes in the population.

Supports PAF from FastGA, wfmash, and minimap2 asm20.
"""

import os
import sys
import re
import json
import logging
import subprocess
import shutil
from collections import defaultdict
from dataclasses import dataclass
from typing import List, Dict, Tuple, Optional, Set

logger = logging.getLogger(__name__)

# Import complexity utilities from Refiner
BIN_DIR = os.path.dirname(os.path.abspath(__file__))
_refiner_dir = os.path.join(BIN_DIR, 'Refiner', 'utils')
if _refiner_dir not in sys.path:
    sys.path.insert(0, _refiner_dir)
from complexity_utils import calculate_shannon_entropy, calculate_dust_score

_CIGAR_OP = re.compile(r'(\d+)([MIDNSHP=X])')


# --- Data classes -----------------------------------------------------------

@dataclass
class PAFRecord:
    """Single PAF alignment record."""
    qname: str
    qlen: int
    qstart: int
    qend: int
    strand: str
    tname: str
    tlen: int
    tstart: int
    tend: int
    matches: int
    block_len: int
    mapq: int
    cigar: str


@dataclass
class InsertionCandidate:
    """A candidate TE insertion from a chain gap or CIGAR."""
    seq_name: str
    start: int       # 0-based
    end: int         # 0-based exclusive
    source: str      # 'gap', 'cigar', or combined like 'cigar+gap' after merge


# --- Main processor ----------------------------------------------------------

class PAFProcessor:
    """Extract TE insertion sequences from all-vs-all PAF alignments and
    build synthetic scaffolds for genome enrichment.

    The processor identifies TE insertions by detecting asymmetric gaps in
    syntenic alignment chains: when one genome has a large gap relative to
    another at a syntenic breakpoint, this suggests an insertion (often a TE).

    For minimap2, CIGAR-level I/D operations provide additional fine-grained
    insertion detection within individual alignments.

    The final output is a synthetic scaffold FASTA: all extracted TE insertions
    concatenated with N-spacers. This can be used either as:
    - A standalone genome input for the three detection engines (reference-free
      mode: --paf + --paf-fasta only), or
    - An enrichment appended to an existing genome (--genome + --paf mode).

    In both modes, the detection engines (LTR, mdl-repeat, RepGraph) discover
    families and build consensus from the population-level sequence content.
    """

    def __init__(self, paf_file: str, paf_fasta: str, output_dir: str,
                 reference_genome: str = None, threads: int = 4,
                 mapping_file: str = None,
                 assembly_fastas: List[str] = None):
        self.paf_file = paf_file
        self.paf_fasta = paf_fasta
        self.output_dir = output_dir
        self.reference_genome = reference_genome  # Optional: skip ref-seq candidates
        self.threads = threads
        self.mapping_file = mapping_file
        self.assembly_fastas = assembly_fastas or []  # Per-assembly FASTAs for base selection

        # --- Core size parameters ---
        self.min_insertion_size = 100        # MITEs ~100bp, SINEs ~100-500bp
        self.max_insertion_size = 30000      # LTR 5-15kb, LINE ~6kb

        # --- Chaining parameters ---
        self.max_chain_gap = 500000          # Max gap between chain members (bp)
        self.min_chain_alignments = 3        # Min alignments per chain
        self.overlap_tolerance = 200         # Overlap tolerance between chain members
        self.min_alignment_length = 5000     # Substantial syntenic anchor
        # Ensure chain has sufficient flanking synteny relative to max insertion
        self.min_chain_aligned = max(50000, 2 * self.max_insertion_size)

        # --- Filtering parameters ---
        self.gap_ratio_threshold = 0.3       # Stricter asymmetry: ~3:1 gap ratio
        self.min_mapq = 20                   # Only for minimap2

        # --- QC parameters ---
        self.cdhit_identity = 0.8
        self.entropy_threshold_long = 1.0    # Shannon k=1, seq >= 100bp
        self.entropy_threshold_short = 0.8   # Shannon k=1, seq < 100bp
        self.dust_threshold = 0.7            # Consistent with Refiner_mdl
        self.max_n_percent = 0.20

        # --- Scaffold construction ---
        self.n_spacer_length = 5000          # N's between insertions in scaffold
        self.max_scaffold_length = 50000000  # 50MB per scaffold

        # --- Runtime state ---
        self.aligner_type = 'unknown'
        self.seq_lengths: Dict[str, int] = {}
        self.ref_seqnames: Set[str] = set()
        self.seq_mapping: Dict[str, Tuple[str, str]] = {}
        self._stats: Dict[str, int] = defaultdict(int)
        self._base_assembly: Optional[str] = None  # Selected base assembly path

        os.makedirs(self.output_dir, exist_ok=True)

    # ------------------------------------------------------------------
    # Public entry point
    # ------------------------------------------------------------------
    def process(self) -> str:
        """Run full PAF TE extraction and build enriched genome.

        Returns:
            In reference-free mode (no reference_genome): path to the
            comprehensive genome FASTA (base assembly + PAF scaffolds).
            In enrichment mode: path to the scaffold-only FASTA for
            appending to an existing genome.
        """
        self._log_config()
        self._validate_paf_format()
        self._validate_and_index()

        grouped, self.aligner_type = self._parse_and_filter_paf()
        logger.info(f"Aligner detected: {self.aligner_type}")

        chains = self._build_all_chains(grouped)
        candidates = self._extract_gap_insertions(chains)

        if self.aligner_type == 'minimap2':
            cigar_cands = self._extract_cigar_insertions(grouped)
            candidates.extend(cigar_cands)

        candidates = self._merge_overlapping(candidates)
        candidates = self._validate_coordinates(candidates)
        logger.info(f"Final candidates after merge+validation: {len(candidates)}")

        scaffold_path = os.path.join(self.output_dir, 'paf_scaffolds.fa')

        if not candidates:
            logger.warning("No TE insertion candidates found in PAF")
            with open(scaffold_path, 'w'):
                pass
            self._write_stats()
            # Reference-free: return base assembly even if no PAF candidates
            if not self.reference_genome and self._base_assembly:
                return self._base_assembly
            return scaffold_path

        raw_fasta = self._extract_sequences(candidates)
        filtered_fasta = self._qc_filter(raw_fasta)
        deduped_fasta = self._cdhit_dedup(filtered_fasta)

        n_seqs = count_fasta_sequences(deduped_fasta)
        self._stats['final_sequences'] = n_seqs

        if n_seqs == 0:
            with open(scaffold_path, 'w'):
                pass
            self._write_stats()
            if not self.reference_genome and self._base_assembly:
                return self._base_assembly
            return scaffold_path

        # Build synthetic scaffolds: concatenate insertions with N-spacers
        self._build_synthetic_scaffolds(deduped_fasta, scaffold_path)

        logger.info(f"PAF enrichment complete: {n_seqs} TE insertions "
                    f"assembled into scaffolds")

        # Cleanup intermediate files
        for tmp_name in ['paf_raw_insertions.fa', 'paf_filtered_insertions.fa',
                         'paf_te_insertions.fa']:
            tmp_path = os.path.join(self.output_dir, tmp_name)
            if os.path.exists(tmp_path):
                os.remove(tmp_path)

        # Reference-free mode: combine base assembly + scaffolds into one genome
        if not self.reference_genome and self._base_assembly:
            combined = os.path.join(self.output_dir, 'pan_genome.fa')
            with open(combined, 'wb') as out:
                with open(self._base_assembly, 'rb') as f:
                    shutil.copyfileobj(f, out)
                # Ensure newline before scaffold
                out.write(b'\n')
                with open(scaffold_path, 'rb') as f:
                    shutil.copyfileobj(f, out)
            n_base = count_fasta_sequences(self._base_assembly)
            n_scaffold = count_fasta_sequences(scaffold_path)
            logger.info(f"Pan-genome constructed: {n_base} base sequences + "
                        f"{n_scaffold} scaffolds")
            self._stats['base_assembly'] = os.path.basename(self._base_assembly)
            self._stats['base_sequences'] = n_base
            self._write_stats()
            return combined

        self._write_stats()
        return scaffold_path

    # ------------------------------------------------------------------
    # Step 0: Validation + indexing
    # ------------------------------------------------------------------
    def _validate_paf_format(self):
        """Validate first few PAF records for format sanity."""
        checked = 0
        with open(self.paf_file) as fh:
            for line in fh:
                if line.startswith('#') or not line.strip():
                    continue
                fields = line.rstrip('\n').split('\t')
                if len(fields) < 12:
                    raise ValueError(
                        f"PAF line has {len(fields)} fields, expected >=12: "
                        f"{line[:100]}")
                try:
                    qstart, qend = int(fields[2]), int(fields[3])
                    qlen = int(fields[1])
                    tstart, tend = int(fields[7]), int(fields[8])
                    tlen = int(fields[6])
                except ValueError as e:
                    raise ValueError(f"Non-integer coordinate in PAF: {e}")
                if qstart > qend or tstart > tend:
                    raise ValueError(
                        f"PAF coordinates inverted: q={qstart}-{qend}, "
                        f"t={tstart}-{tend}")
                if qend > qlen or tend > tlen:
                    raise ValueError(
                        f"PAF coordinates exceed sequence length: "
                        f"q={qend}>{qlen} or t={tend}>{tlen}")
                checked += 1
                if checked >= 10:
                    break
        logger.info(f"PAF format validation passed ({checked} records checked)")

    def _validate_and_index(self):
        """Index FASTA files and build metadata."""
        # Index multi-assembly FASTA
        fai = self.paf_fasta + '.fai'
        if not os.path.exists(fai):
            subprocess.run(['samtools', 'faidx', self.paf_fasta],
                           check=True, capture_output=True)
        self.seq_lengths = self._parse_fai(fai)
        logger.info(f"Multi-assembly FASTA: {len(self.seq_lengths)} sequences")

        # Identify reference/base seqnames (used to skip redundant candidates)
        if self.reference_genome:
            # Enrichment mode: skip candidates already in the provided genome
            ref_fai = self.reference_genome + '.fai'
            if not os.path.exists(ref_fai):
                subprocess.run(['samtools', 'faidx', self.reference_genome],
                               check=True, capture_output=True)
            self.ref_seqnames = set(self._parse_fai(ref_fai).keys())
            logger.info(f"Enrichment mode: {len(self.ref_seqnames)} reference sequences")
        elif self.assembly_fastas:
            # Reference-free mode: auto-select largest assembly as base
            self._base_assembly = self._select_base_assembly()
            base_fai = self._base_assembly + '.fai'
            if not os.path.exists(base_fai):
                subprocess.run(['samtools', 'faidx', self._base_assembly],
                               check=True, capture_output=True)
            self.ref_seqnames = set(self._parse_fai(base_fai).keys())
            logger.info(f"Reference-free mode: base assembly "
                        f"{os.path.basename(self._base_assembly)} "
                        f"({len(self.ref_seqnames)} sequences)")
        else:
            logger.info("Reference-free mode: no base assembly, "
                        "extracting from all assemblies")

        # Load chromosome mapping if provided
        if self.mapping_file:
            self._load_mapping()

    def _select_base_assembly(self) -> str:
        """Select the largest assembly as the base genome.

        In reference-free mode, one assembly serves as the base genome for
        the three detection engines. PAF-derived TE insertions from all
        OTHER assemblies are appended as scaffolds, enriching the base
        with population-level TE content.

        The largest assembly is chosen because it likely has the most
        complete TE representation (fewer collapsed repeats).
        """
        best_path = None
        best_size = 0
        for fasta in self.assembly_fastas:
            if not os.path.exists(fasta):
                logger.warning(f"Assembly FASTA not found: {fasta}")
                continue
            size = os.path.getsize(fasta)
            if size > best_size:
                best_size = size
                best_path = fasta
        if not best_path:
            raise ValueError("No valid assembly FASTA found in assembly_fastas")
        logger.info(f"Selected base assembly: {os.path.basename(best_path)} "
                    f"({best_size / 1e6:.1f} Mb)")
        return best_path

    def _parse_fai(self, fai_path: str) -> Dict[str, int]:
        """Parse .fai index: returns {seqname: length}."""
        result = {}
        with open(fai_path) as f:
            for line in f:
                parts = line.split('\t')
                if len(parts) >= 2:
                    result[parts[0]] = int(parts[1])
        return result

    def _load_mapping(self):
        """Load chromosome mapping TSV: seqname\\tassembly_id\\tchr_id."""
        with open(self.mapping_file) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split('\t')
                if len(parts) >= 3:
                    self.seq_mapping[parts[0]] = (parts[1], parts[2])
        logger.info(f"Loaded chromosome mapping: {len(self.seq_mapping)} entries")

    # ------------------------------------------------------------------
    # Step 1: Parse + filter PAF
    # ------------------------------------------------------------------
    def _parse_and_filter_paf(self) -> Tuple[Dict, str]:
        """Stream-parse PAF, filter, detect aligner, group by (qname, tname).

        Both directions of each pair are kept to preserve haplotype-specific
        TE insertions. Deduplication occurs at the candidate level via
        _merge_overlapping() and at the sequence level via CD-HIT-EST.
        """
        grouped = defaultdict(list)
        has_tp = False
        total = 0
        kept = 0
        skipped_short_fields = 0

        with open(self.paf_file) as fh:
            for line in fh:
                if line.startswith('#'):
                    continue
                fields = line.rstrip('\n').split('\t')
                if len(fields) < 12:
                    skipped_short_fields += 1
                    continue
                total += 1

                qname, qlen = fields[0], int(fields[1])
                qstart, qend = int(fields[2]), int(fields[3])
                strand = fields[4]
                tname, tlen = fields[5], int(fields[6])
                tstart, tend = int(fields[7]), int(fields[8])
                matches, block_len = int(fields[9]), int(fields[10])
                mapq = int(fields[11])

                # Filter: self-alignment
                if qname == tname:
                    continue

                # Filter: alignment too short
                aln_len = qend - qstart
                if aln_len < self.min_alignment_length:
                    continue

                # Parse tags for aligner detection and CIGAR
                cigar = ''
                record_has_tp = False
                for f in fields[12:]:
                    if f.startswith('tp:'):
                        record_has_tp = True
                        has_tp = True  # global flag for aligner detection
                    elif f.startswith('cg:Z:'):
                        cigar = f[5:]

                # minimap2: per-record mapq filter (only if THIS record has tp: tag)
                if record_has_tp and mapq < self.min_mapq:
                    continue

                # Mapping-based filters (when mapping file is provided)
                if self.seq_mapping:
                    q_map = self.seq_mapping.get(qname)
                    t_map = self.seq_mapping.get(tname)
                    if q_map and t_map:
                        # Same assembly → paralog alignment, skip
                        if q_map[0] == t_map[0]:
                            continue
                        # Different chromosomes → not homologous, skip
                        if q_map[1] != t_map[1]:
                            continue

                rec = PAFRecord(qname, qlen, qstart, qend, strand,
                                tname, tlen, tstart, tend,
                                matches, block_len, mapq, cigar)

                # Keep both directions — dedup at candidate level
                grouped[(qname, tname)].append(rec)
                kept += 1

        aligner = 'minimap2' if has_tp else 'fastga'
        if skipped_short_fields > 0:
            logger.warning(f"Skipped {skipped_short_fields} malformed PAF lines (<12 fields)")
        logger.info(f"PAF: {total} total, {kept} after filter, "
                    f"{len(grouped)} sequence pairs")
        self._stats['paf_total_records'] = total
        self._stats['paf_kept_records'] = kept
        self._stats['paf_sequence_pairs'] = len(grouped)
        return dict(grouped), aligner

    # ------------------------------------------------------------------
    # Step 2: Synteny chaining
    # ------------------------------------------------------------------
    def _build_all_chains(self, grouped: Dict) -> List[List[PAFRecord]]:
        """Build collinear chains for all (query, target) pairs.

        Separates alignments by strand, then applies greedy collinear chaining
        sorted by target coordinate. Only chains meeting minimum alignment
        count and total aligned base thresholds are retained.
        """
        all_chains = []
        for (qname, tname), records in grouped.items():
            plus = [r for r in records if r.strand == '+']
            minus = [r for r in records if r.strand == '-']
            for strand_recs in [plus, minus]:
                if len(strand_recs) < self.min_chain_alignments:
                    continue
                if len(strand_recs) > 10000:
                    logger.warning(
                        f"Large alignment set for {qname}:{tname} "
                        f"({len(strand_recs)} records), chaining may be slow")
                chains = self._greedy_chain(strand_recs)
                all_chains.extend(chains)
        logger.info(f"Built {len(all_chains)} synteny chains")
        self._stats['chains_built'] = len(all_chains)
        return all_chains

    def _greedy_chain(self, records: List[PAFRecord]) -> List[List[PAFRecord]]:
        """Greedy collinear chaining sorted by target start.

        O(n^2) in worst case, but acceptable for typical inter-assembly pairs
        (hundreds to low thousands of alignments per pair). Records are sorted
        by target start; early termination on large target gaps limits the
        effective inner loop.
        """
        records = sorted(records, key=lambda r: r.tstart)
        chains = []
        used = [False] * len(records)

        for i in range(len(records)):
            if used[i]:
                continue
            chain = [records[i]]
            used[i] = True

            for j in range(i + 1, len(records)):
                if used[j]:
                    continue
                prev = chain[-1]
                curr = records[j]

                # Target must be monotonically increasing
                t_gap = curr.tstart - prev.tend
                if t_gap > self.max_chain_gap:
                    break  # Sorted by tstart: all subsequent are further
                if t_gap < -self.overlap_tolerance:
                    continue

                # Query must also be monotonically ordered (strand-aware)
                if curr.strand == '+':
                    q_gap = curr.qstart - prev.qend
                else:
                    # Minus strand: query coords decrease along the chain.
                    # Records sorted by tstart → query positions decrease.
                    # prev covers higher query coords than curr.
                    q_gap = prev.qstart - curr.qend

                if q_gap < -self.overlap_tolerance or q_gap > self.max_chain_gap:
                    continue

                chain.append(curr)
                used[j] = True

            if len(chain) >= self.min_chain_alignments:
                total_aligned = sum(r.tend - r.tstart for r in chain)
                if total_aligned >= self.min_chain_aligned:
                    chains.append(chain)

        return chains

    # ------------------------------------------------------------------
    # Step 3: Gap-based insertion extraction
    # ------------------------------------------------------------------
    def _extract_gap_insertions(self, chains: List[List[PAFRecord]]) -> List[InsertionCandidate]:
        """Extract TE candidates from asymmetric gaps between chain members.

        For each consecutive pair of alignments in a chain, computes the gap
        on both query and target sides. An asymmetric gap (ratio < threshold)
        indicates an insertion on the side with the larger gap.
        """
        candidates = []

        for chain in chains:
            for i in range(len(chain) - 1):
                curr = chain[i]
                nxt = chain[i + 1]

                t_gap = nxt.tstart - curr.tend
                if curr.strand == '+':
                    q_gap = nxt.qstart - curr.qend
                else:
                    q_gap = curr.qstart - nxt.qend

                max_gap = max(abs(t_gap), abs(q_gap))
                min_gap = min(abs(t_gap), abs(q_gap))

                if max_gap < self.min_insertion_size or max_gap > self.max_insertion_size:
                    continue

                # Asymmetry check
                ratio = min_gap / max_gap if max_gap > 0 else 1.0
                if ratio > self.gap_ratio_threshold:
                    continue

                # Determine which side has the insertion
                if abs(q_gap) > abs(t_gap):
                    # Insertion on query side — require positive gap
                    if curr.strand == '+':
                        ins_start, ins_end = curr.qend, nxt.qstart
                    else:
                        ins_start, ins_end = nxt.qend, curr.qstart
                    if ins_end <= ins_start:
                        continue
                    seq_name = curr.qname
                else:
                    # Insertion on target side — require positive gap
                    ins_start, ins_end = curr.tend, nxt.tstart
                    if ins_end <= ins_start:
                        continue
                    seq_name = curr.tname

                ins_len = ins_end - ins_start
                if ins_len < self.min_insertion_size or ins_len > self.max_insertion_size:
                    continue

                # Skip candidates on reference sequences (already in detection paths)
                if seq_name in self.ref_seqnames:
                    continue

                candidates.append(InsertionCandidate(
                    seq_name=seq_name, start=ins_start, end=ins_end, source='gap'))

        logger.info(f"Gap insertions: {len(candidates)} candidates")
        self._stats['gap_candidates'] = len(candidates)
        return candidates

    # ------------------------------------------------------------------
    # Step 3b: CIGAR-based insertion extraction (minimap2 only)
    # ------------------------------------------------------------------
    def _extract_cigar_insertions(self, grouped: Dict) -> List[InsertionCandidate]:
        """Walk CIGAR strings for large I/D operations (minimap2 only).

        PAF CIGAR convention: qstart/qend are always on the forward strand
        regardless of alignment strand. CIGAR walks from qstart to qend,
        consuming exactly (qend - qstart) query bases. No reverse-complement
        coordinate transformation is needed.
        """
        candidates = []
        for records in grouped.values():
            for rec in records:
                if not rec.cigar:
                    continue

                # PAF CIGAR should not contain S/H ops (those are SAM concepts).
                # If present, skip the record as the input may be malformed.
                if 'S' in rec.cigar or 'H' in rec.cigar:
                    logger.debug(
                        f"Skipping record with S/H ops in CIGAR: "
                        f"{rec.qname}:{rec.qstart}-{rec.qend}")
                    continue

                q_pos = rec.qstart
                t_pos = rec.tstart

                for m in _CIGAR_OP.finditer(rec.cigar):
                    length = int(m.group(1))
                    op = m.group(2)

                    if op in ('M', '=', 'X'):
                        q_pos += length
                        t_pos += length
                    elif op == 'I':
                        if self.min_insertion_size <= length <= self.max_insertion_size:
                            seq_name = rec.qname
                            if seq_name not in self.ref_seqnames:
                                candidates.append(InsertionCandidate(
                                    seq_name=seq_name, start=q_pos,
                                    end=q_pos + length, source='cigar'))
                        q_pos += length
                    elif op == 'D':
                        if self.min_insertion_size <= length <= self.max_insertion_size:
                            seq_name = rec.tname
                            if seq_name not in self.ref_seqnames:
                                candidates.append(InsertionCandidate(
                                    seq_name=seq_name, start=t_pos,
                                    end=t_pos + length, source='cigar'))
                        t_pos += length
                    elif op == 'N':
                        t_pos += length
                    # P op: padding, no coordinate advance

                # Consistency check
                if q_pos != rec.qend:
                    logger.warning(
                        f"CIGAR walk mismatch for {rec.qname}: "
                        f"expected q_pos={rec.qend}, got {q_pos}")

        logger.info(f"CIGAR insertions: {len(candidates)} candidates")
        self._stats['cigar_candidates'] = len(candidates)
        return candidates

    # ------------------------------------------------------------------
    # Step 4: Merge overlapping candidates
    # ------------------------------------------------------------------
    def _merge_overlapping(self, candidates: List[InsertionCandidate]) -> List[InsertionCandidate]:
        """Merge overlapping intervals on the same sequence.

        Tracks source provenance through merging (e.g. 'cigar+gap').
        Filters merged intervals that exceed max_insertion_size.
        """
        if not candidates:
            return []

        by_seq: Dict[str, List[InsertionCandidate]] = defaultdict(list)
        for c in candidates:
            by_seq[c.seq_name].append(c)

        merged = []
        for seq_name, cands in by_seq.items():
            if not cands:
                continue
            cands.sort(key=lambda c: c.start)
            cur_start = cands[0].start
            cur_end = cands[0].end
            cur_sources: Set[str] = {cands[0].source}

            for c in cands[1:]:
                if c.start <= cur_end:
                    cur_end = max(cur_end, c.end)
                    cur_sources.add(c.source)
                else:
                    size = cur_end - cur_start
                    if self.min_insertion_size <= size <= self.max_insertion_size:
                        merged.append(InsertionCandidate(
                            seq_name=seq_name, start=cur_start,
                            end=cur_end,
                            source='+'.join(sorted(cur_sources))))
                    cur_start, cur_end = c.start, c.end
                    cur_sources = {c.source}

            size = cur_end - cur_start
            if self.min_insertion_size <= size <= self.max_insertion_size:
                merged.append(InsertionCandidate(
                    seq_name=seq_name, start=cur_start,
                    end=cur_end,
                    source='+'.join(sorted(cur_sources))))

        self._stats['merged_candidates'] = len(merged)
        return merged

    # ------------------------------------------------------------------
    # Step 4b: Coordinate validation
    # ------------------------------------------------------------------
    def _validate_coordinates(self, candidates: List[InsertionCandidate]) -> List[InsertionCandidate]:
        """Clamp coordinates to sequence boundaries and filter invalid candidates."""
        valid = []
        for c in candidates:
            seq_len = self.seq_lengths.get(c.seq_name)
            if seq_len is None:
                logger.debug(f"Skipping candidate on unknown sequence: {c.seq_name}")
                continue
            start = max(0, c.start)
            end = min(seq_len, c.end)
            if end - start < self.min_insertion_size:
                continue
            valid.append(InsertionCandidate(
                seq_name=c.seq_name, start=start, end=end, source=c.source))
        if len(valid) < len(candidates):
            logger.info(f"Coordinate validation: {len(candidates)} -> {len(valid)} candidates")
        self._stats['validated_candidates'] = len(valid)
        return valid

    # ------------------------------------------------------------------
    # Step 5: Sequence extraction
    # ------------------------------------------------------------------
    def _extract_sequences(self, candidates: List[InsertionCandidate]) -> str:
        """Batch-extract sequences with samtools faidx.

        Returns path to the raw extraction FASTA. Uses 1-based inclusive
        coordinates for samtools; candidate coordinates are 0-based half-open.
        """
        out_fasta = os.path.join(self.output_dir, 'paf_raw_insertions.fa')
        BATCH = 500
        global_idx = 0
        failed_regions = 0

        with open(out_fasta, 'w') as out:
            for i in range(0, len(candidates), BATCH):
                batch = candidates[i:i + BATCH]
                # samtools faidx uses 1-based inclusive coords
                regions = [f"{c.seq_name}:{c.start + 1}-{c.end}" for c in batch]
                # Build region->candidate lookup for header-based matching
                region_to_cand = {}
                for c, region in zip(batch, regions):
                    region_to_cand[region] = c
                cmd = ['samtools', 'faidx', self.paf_fasta] + regions
                try:
                    result = subprocess.run(cmd, capture_output=True, text=True,
                                           timeout=300)
                    # samtools may return non-zero if SOME regions fail but
                    # still output valid sequences for others; process stdout
                    # regardless of returncode using header-based matching
                    if not result.stdout.strip():
                        if result.returncode != 0:
                            logger.warning(
                                f"samtools faidx batch failed: "
                                f"{result.stderr[:200]}")
                        continue
                    # Match output headers to candidates by parsing samtools
                    # header format ">seqname:start-end". Buffer each entry
                    # to skip headers with empty sequences (samtools 1.4.1
                    # emits a header for missing regions but no sequence).
                    pending_header = None
                    pending_cand = None
                    seq_lines = []
                    batch_extracted = 0

                    for line in result.stdout.split('\n'):
                        if line.startswith('>'):
                            # Flush previous entry if it had sequence data
                            if pending_cand and seq_lines:
                                global_idx += 1
                                out.write(
                                    f">paf_ins_{global_idx} "
                                    f"{pending_cand.seq_name}:"
                                    f"{pending_cand.start}-"
                                    f"{pending_cand.end}"
                                    f"(0-based) "
                                    f"len={pending_cand.end - pending_cand.start}"
                                    f" source={pending_cand.source}\n")
                                out.write('\n'.join(seq_lines) + '\n')
                                batch_extracted += 1
                            # Start new entry
                            sam_region = line[1:].split()[0]
                            pending_cand = region_to_cand.get(sam_region)
                            seq_lines = []
                            if not pending_cand:
                                logger.debug(
                                    f"Unmatched samtools region: {sam_region}")
                        elif line.strip() and pending_cand is not None:
                            seq_lines.append(line.strip().upper())

                    # Flush last entry
                    if pending_cand and seq_lines:
                        global_idx += 1
                        out.write(
                            f">paf_ins_{global_idx} "
                            f"{pending_cand.seq_name}:"
                            f"{pending_cand.start}-{pending_cand.end}"
                            f"(0-based) "
                            f"len={pending_cand.end - pending_cand.start}"
                            f" source={pending_cand.source}\n")
                        out.write('\n'.join(seq_lines) + '\n')
                        batch_extracted += 1

                    failed_regions += len(batch) - batch_extracted
                except subprocess.TimeoutExpired:
                    logger.warning("samtools faidx batch timed out")

        if failed_regions > 0:
            logger.warning(f"samtools: {failed_regions} regions could not be extracted")

        n = count_fasta_sequences(out_fasta)
        logger.info(f"Extracted {n} sequences")
        self._stats['extracted_sequences'] = n
        return out_fasta

    # ------------------------------------------------------------------
    # Step 6: QC filter
    # ------------------------------------------------------------------
    def _qc_filter(self, fasta_path: str) -> str:
        """Filter by entropy, DUST, N-content, length."""
        out_path = os.path.join(self.output_dir, 'paf_filtered_insertions.fa')
        stats = defaultdict(int)

        with open(fasta_path) as fin, open(out_path, 'w') as fout:
            header = ''
            seq_parts = []
            for line in fin:
                if line.startswith('>'):
                    if header:
                        self._qc_one(header, ''.join(seq_parts), fout, stats)
                    header = line.rstrip('\n')
                    seq_parts = []
                else:
                    seq_parts.append(line.strip())
            if header:
                self._qc_one(header, ''.join(seq_parts), fout, stats)

        logger.info(f"QC: {stats['total']} total, {stats['passed']} passed "
                    f"(short={stats['short']}, long={stats['long']}, "
                    f"low_ent={stats['low_entropy']}, dust={stats['high_dust']}, "
                    f"high_N={stats['high_n']})")
        self._stats['qc_passed'] = stats['passed']
        self._stats['qc_total'] = stats['total']
        return out_path

    def _qc_one(self, header: str, seq: str, fout, stats: dict):
        """Apply QC filters to a single sequence."""
        stats['total'] += 1
        slen = len(seq)
        if slen < self.min_insertion_size:
            stats['short'] += 1
            return
        if slen > self.max_insertion_size:
            stats['long'] += 1
            return
        n_frac = seq.count('N') / slen if slen > 0 else 1.0
        if n_frac > self.max_n_percent:
            stats['high_n'] += 1
            return
        # Length-dependent entropy threshold (aligned with Refiner_mdl)
        entropy = calculate_shannon_entropy(seq, k=1)
        threshold = (self.entropy_threshold_long if slen >= 100
                     else self.entropy_threshold_short)
        if entropy < threshold:
            stats['low_entropy'] += 1
            return
        dust = calculate_dust_score(seq)
        if dust > self.dust_threshold:
            stats['high_dust'] += 1
            return
        stats['passed'] += 1
        fout.write(header + '\n' + seq + '\n')

    # ------------------------------------------------------------------
    # Step 7: CD-HIT-EST dedup
    # ------------------------------------------------------------------
    def _cdhit_dedup(self, fasta_path: str) -> str:
        """Deduplicate with cd-hit-est.

        Cross-haplotype deduplication is handled here (80% identity).
        Final cross-source deduplication occurs in combine_results().
        """
        n_input = count_fasta_sequences(fasta_path)
        if n_input <= 1:
            final = os.path.join(self.output_dir, 'paf_te_insertions.fa')
            shutil.copy2(fasta_path, final)
            return final

        out_path = os.path.join(self.output_dir, 'paf_te_insertions.fa')
        cmd = [
            'cd-hit-est',
            '-i', fasta_path,
            '-o', out_path,
            '-c', str(self.cdhit_identity),
            '-aS', '0.8',
            '-aL', '0.8',
            '-G', '0',
            '-n', '8',
            '-r', '1',
            '-mask', 'NX',
            '-M', '0',
            '-T', str(self.threads),
            '-d', '0'
        ]
        try:
            subprocess.run(cmd, check=True, capture_output=True,
                           text=True, timeout=3600)
        except subprocess.CalledProcessError as e:
            logger.error(f"cd-hit-est failed: {e.stderr[:300]}")
            shutil.copy2(fasta_path, out_path)
            return out_path

        n_output = count_fasta_sequences(out_path)
        logger.info(f"CD-HIT-EST: {n_input} -> {n_output}")

        # Cleanup cluster file
        clstr = out_path + '.clstr'
        if os.path.exists(clstr):
            os.remove(clstr)

        return out_path

    # ------------------------------------------------------------------
    # Step 8: Build synthetic scaffolds
    # ------------------------------------------------------------------
    def _build_synthetic_scaffolds(self, fasta_path: str, output_path: str):
        """Concatenate TE insertions with N-spacers into synthetic scaffolds.

        Each scaffold is named paf_scaffold_N and contains multiple TE
        insertions separated by long N-gaps (default 5000 N's). This prevents
        detection tools from connecting adjacent insertions while keeping them
        in a single contiguous sequence that the tools can process normally.

        Scaffolds are split at max_scaffold_length to stay within tool limits.
        """
        spacer = 'N' * self.n_spacer_length
        line_width = 80

        # Read all sequences
        sequences = []
        with open(fasta_path) as f:
            header = ''
            seq_parts = []
            for line in f:
                if line.startswith('>'):
                    if header and seq_parts:
                        sequences.append(''.join(seq_parts))
                    header = line.rstrip('\n')
                    seq_parts = []
                else:
                    seq_parts.append(line.strip())
            if header and seq_parts:
                sequences.append(''.join(seq_parts))

        if not sequences:
            with open(output_path, 'w'):
                pass
            return

        # Build scaffolds with N-spacers, split at max length
        # Each entry is (joined_sequence, n_parts)
        scaffolds = []
        current_parts = []
        current_length = 0

        for seq in sequences:
            addition = len(seq) + (self.n_spacer_length if current_parts else 0)
            if current_length + addition > self.max_scaffold_length and current_parts:
                scaffolds.append((spacer.join(current_parts), len(current_parts)))
                current_parts = []
                current_length = 0
            if current_parts:
                current_length += self.n_spacer_length
            current_parts.append(seq)
            current_length += len(seq)

        if current_parts:
            scaffolds.append((spacer.join(current_parts), len(current_parts)))

        # Write scaffolds
        with open(output_path, 'w') as out:
            for i, (scaffold_seq, n_parts) in enumerate(scaffolds, 1):
                out.write(f">paf_scaffold_{i} "
                          f"n_insertions={n_parts} "
                          f"len={len(scaffold_seq)}\n")
                for j in range(0, len(scaffold_seq), line_width):
                    out.write(scaffold_seq[j:j + line_width] + '\n')

        total_bp = sum(len(seq) for seq, _ in scaffolds)
        self._stats['scaffold_count'] = len(scaffolds)
        self._stats['scaffold_total_bp'] = total_bp
        logger.info(f"Built {len(scaffolds)} synthetic scaffold(s), "
                    f"total {total_bp / 1e6:.1f} Mb")

    # ------------------------------------------------------------------
    # Utilities
    # ------------------------------------------------------------------
    def _log_config(self):
        """Log current configuration for reproducibility."""
        logger.info(
            f"PAF Processor v1.0 | min_ins={self.min_insertion_size} "
            f"max_ins={self.max_insertion_size} min_aln={self.min_alignment_length} "
            f"min_chain_aln={self.min_chain_aligned} gap_ratio={self.gap_ratio_threshold} "
            f"entropy={self.entropy_threshold_long}/{self.entropy_threshold_short} "
            f"dust={self.dust_threshold}")

    def _write_stats(self):
        """Write processing statistics to JSON."""
        self._stats['aligner_type'] = self.aligner_type
        stats_path = os.path.join(self.output_dir, 'paf_processing_stats.json')
        with open(stats_path, 'w') as f:
            json.dump(dict(self._stats), f, indent=2)


def count_fasta_sequences(fasta: str) -> int:
    """Count sequences in a FASTA file."""
    if not os.path.exists(fasta):
        return 0
    count = 0
    with open(fasta) as f:
        for line in f:
            if line.startswith('>'):
                count += 1
    return count


# ------------------------------------------------------------------
# Standalone CLI
# ------------------------------------------------------------------
def main():
    """Standalone CLI for PAF TE insertion extraction."""
    import argparse
    parser = argparse.ArgumentParser(
        description='Extract TE insertion sequences from all-vs-all PAF alignments')
    parser.add_argument('--paf', required=True, help='PAF alignment file')
    parser.add_argument('--fasta', required=True, help='Multi-assembly FASTA')
    parser.add_argument('--ref', help='Reference genome FASTA (enrichment mode)')
    parser.add_argument('--assemblies', nargs='+',
                        help='Per-assembly FASTA files (reference-free mode: '
                             'largest is auto-selected as base genome)')
    parser.add_argument('--outdir', required=True, help='Output directory')
    parser.add_argument('--mapping', help='Chromosome mapping TSV')
    parser.add_argument('--threads', type=int, default=4)
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO,
                        format='[%(asctime)s] %(levelname)s: %(message)s')

    proc = PAFProcessor(
        paf_file=args.paf,
        paf_fasta=args.fasta,
        output_dir=args.outdir,
        reference_genome=args.ref,
        threads=args.threads,
        mapping_file=args.mapping,
        assembly_fastas=args.assemblies or []
    )
    result = proc.process()
    print(f"Output: {result}")


if __name__ == '__main__':
    main()
