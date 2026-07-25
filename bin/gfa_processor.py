#!/usr/bin/env python3
"""
GFA TE Enrichment Preprocessor for Pan_TE

Extracts TE insertion candidates from GFA (Graphical Fragment Assembly) pangenome
graph files. Non-reference segments (or chains thereof) in the TE size range are
extracted as candidate TE insertions.

Supports two GFA dialects with auto-detection:
  - rGFA (minigraph-cactus): coarse segments with SR/SN/SO tags, L-line topology
  - PGGB (smoothxg): fine segments with P-line paths, no optional S-line tags

Mirrors the PAFProcessor interface: enrichment mode returns scaffold path;
reference-free mode returns a combined pan_genome FASTA (reconstructed reference
+ scaffolds).
"""

import os
import sys
import gc
import json
import logging
import shutil
import subprocess
from collections import defaultdict
from typing import Dict, List, Set, Tuple

logger = logging.getLogger(__name__)

BIN_DIR = os.path.dirname(os.path.abspath(__file__))
_refiner_dir = os.path.join(BIN_DIR, 'Refiner', 'utils')
if _refiner_dir not in sys.path:
    sys.path.insert(0, _refiner_dir)
from complexity_utils import calculate_shannon_entropy, calculate_dust_score
from gfa_alleles import (AlleleCandidate, SQLiteAlleleStore,
                         deduplicate_locus_alleles,
                         extract_pggb_alleles_stream, parse_p_steps,
                         write_manifest)
from rgfa_gaf import (RgfaGafCallStore, extract_rgfa_gaf_alleles,
                      load_rgfa_graph, populate_rgfa_gaf_locus_calls)

_COMPLEMENT = str.maketrans('ACGTacgtRYSWKMBDHVNryswkmbdhvn',
                            'TGCAtgcaYRSWMKVHDBNyrswmkvhdbn')


def _reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def count_fasta_sequences(fasta: str) -> int:
    """Count sequences in a FASTA file. Return 0 if file doesn't exist."""
    if not os.path.exists(fasta):
        return 0
    count = 0
    with open(fasta) as f:
        for line in f:
            if line.startswith('>'):
                count += 1
    return count


class GfaProcessor:
    """Extract TE insertion sequences from pangenome GFA graphs and build
    synthetic scaffolds for genome enrichment.

    For rGFA (minigraph-cactus), non-reference segments are identified by
    SR:i: rank > 0, and maximal linear chains of such segments are walked
    through the link topology. Chains in the TE size range are emitted.

    For PGGB (smoothxg), the default allele mode pairs each alternate
    chromosome with its reference chromosome and extracts locus-aware
    alternate alleles between ordered, unique reference anchors.  Long
    alleles are tiled with overlap and every occurrence remains traceable in
    a manifest.  The legacy global non-reference-node mode remains available.

    Candidates undergo sequence QC.  Allele mode collapses only exact
    sequence duplicates from the same locus and preserves sample/copy support;
    legacy mode retains the historical CD-HIT-EST behavior.
    """

    def __init__(self, gfa_file: str, output_dir: str,
                 reference_genome: str = None, threads: int = 4,
                 gfa_ref: str = None, gfa_mode: str = 'alleles',
                 window_overlap: int = 1000, spacer_length: int = 0,
                 keep_intermediates: bool = True,
                 gfa_gaf: str = None, gfa_gaf_fasta: str = None,
                 gfa_gaf_min_mapq: int = 20,
                 gfa_gaf_min_identity: float = 0.90):
        """
        Args:
            gfa_file: Path to GFA file (rGFA or PGGB).
            output_dir: Working directory for intermediate and output files.
            reference_genome: If provided, enrichment mode; None means
                reference-free mode.
            threads: Thread count for external tools.
            gfa_ref: Optional explicit reference path/sample name for PGGB.
            gfa_mode: ``alleles`` (default) or historical ``legacy`` mode.
            window_overlap: Overlap used when tiling alleles over 30 kb.
            spacer_length: N spacer length; 0 emits one FASTA record per allele.
            keep_intermediates: Preserve raw/unique/filtered FASTA and manifest.
        """
        if gfa_mode not in {'alleles', 'legacy'}:
            raise ValueError("gfa_mode must be 'alleles' or 'legacy'")
        if window_overlap < 0 or window_overlap >= 30000:
            raise ValueError('window_overlap must satisfy 0 <= overlap < 30000')
        if spacer_length < 0:
            raise ValueError('spacer_length must be non-negative')

        self.gfa_file = gfa_file
        self.output_dir = output_dir
        self.reference_genome = reference_genome
        self.threads = threads
        self.gfa_ref = gfa_ref
        self.gfa_mode = gfa_mode
        self.window_overlap = window_overlap
        self.keep_intermediates = keep_intermediates
        self.gfa_gaf = os.path.abspath(gfa_gaf) if gfa_gaf else None
        self.gfa_gaf_fasta = (os.path.abspath(gfa_gaf_fasta)
                              if gfa_gaf_fasta else None)
        self.gfa_gaf_min_mapq = int(gfa_gaf_min_mapq)
        self.gfa_gaf_min_identity = float(gfa_gaf_min_identity)
        if self.gfa_gaf_min_mapq < 0:
            raise ValueError('gfa_gaf_min_mapq must be non-negative')
        if not 0 <= self.gfa_gaf_min_identity <= 1:
            raise ValueError('gfa_gaf_min_identity must be between 0 and 1')

        # Size parameters
        self.min_insertion_size = 100
        self.max_insertion_size = 30000

        # QC parameters (same as PAFProcessor)
        self.entropy_threshold_long = 1.0
        self.entropy_threshold_short = 0.8
        self.dust_threshold = 0.7
        self.max_n_percent = 0.20
        self.cdhit_identity = 0.8

        # Candidate packing parameters.  Zero avoids artificial sequence and
        # emits one independent FASTA record for every representative allele.
        self.n_spacer_length = spacer_length
        self.max_scaffold_length = 50000000

        # Runtime state
        self.gfa_format: str = 'unknown'
        self._stats: Dict = defaultdict(int)
        self._allele_occurrences: List[AlleleCandidate] = []
        self._allele_representatives: List[AlleleCandidate] = []
        self._allele_store = None
        self._allele_mode_active = False
        self._rgfa_call_store = None
        self._rgfa_reference_samples = set()

        os.makedirs(self.output_dir, exist_ok=True)

    # ------------------------------------------------------------------
    # Public entry point
    # ------------------------------------------------------------------
    def process(self) -> str:
        """Run full GFA TE extraction and build enriched genome.

        Returns:
            In enrichment mode: path to gfa_scaffolds.fa
            In reference-free mode: path to gfa_pan_genome.fa
        """
        try:
            return self._process()
        except BaseException:
            self._cleanup_open_allele_stores()
            raise

    def _process(self) -> str:
        """Implementation of :meth:`process`, separated for failure cleanup."""
        self._log_config()
        self.gfa_format = self._detect_format()
        logger.info(f"GFA format detected: {self.gfa_format}")

        if self.gfa_format == 'rGFA':
            if self.gfa_mode == 'legacy':
                candidates, ref_bp = self._process_rgfa()
            else:
                candidates, ref_bp = self._process_rgfa_gaf()
        else:
            candidates, ref_bp = self._process_pggb()

        self._stats['ref_total_bp'] = ref_bp

        scaffold_path = os.path.join(self.output_dir, 'gfa_scaffolds.fa')

        if not candidates:
            logger.warning("No TE insertion candidates found in GFA")
            with open(scaffold_path, 'w'):
                pass
            self._write_allele_manifest()
            self._cleanup_candidate_intermediates()
            self._write_stats()
            if not self.reference_genome:
                ref_fa = os.path.join(self.output_dir, 'gfa_reference.fa')
                if os.path.exists(ref_fa):
                    return ref_fa
            return scaffold_path

        if self._allele_mode_active:
            # Aggregate before FASTA output: occurrence sequences that are
            # identical at one locus never incur cohort-scale raw FASTA I/O.
            if not self._allele_representatives:
                self._allele_representatives = deduplicate_locus_alleles(
                    self._allele_occurrences, self._stats)
            unique_fasta = os.path.join(
                self.output_dir, 'gfa_unique_locus_alleles.fa')
            self._write_allele_candidates(
                self._allele_representatives, unique_fasta)
            filtered_fasta = self._qc_filter(unique_fasta)
            deduped_fasta = filtered_fasta
            self._mark_allele_qc(filtered_fasta)
            self._stats['candidates_merged'] = len(
                self._allele_representatives)
        else:
            raw_fasta = self._write_raw_candidates(candidates)
            filtered_fasta = self._qc_filter(raw_fasta)
            deduped_fasta = self._cdhit_dedup(filtered_fasta)

        n_seqs = count_fasta_sequences(deduped_fasta)
        self._stats['final_sequences'] = n_seqs

        if n_seqs == 0:
            with open(scaffold_path, 'w'):
                pass
            self._write_allele_manifest()
            self._cleanup_candidate_intermediates()
            self._write_stats()
            if not self.reference_genome:
                ref_fa = os.path.join(self.output_dir, 'gfa_reference.fa')
                if os.path.exists(ref_fa):
                    return ref_fa
            return scaffold_path

        self._build_synthetic_scaffolds(deduped_fasta, scaffold_path)
        logger.info(f"GFA enrichment complete: {n_seqs} TE insertions "
                    f"assembled into scaffolds")

        # Keep provenance by default.  Large deployments may opt into cleanup
        # after archiving the manifest and summary statistics.
        self._cleanup_candidate_intermediates()

        # Reference-free mode: combine reconstructed reference + scaffolds
        if not self.reference_genome:
            ref_fa = os.path.join(self.output_dir, 'gfa_reference.fa')
            if os.path.exists(ref_fa):
                combined = os.path.join(self.output_dir, 'gfa_pan_genome.fa')
                with open(combined, 'wb') as out:
                    with open(ref_fa, 'rb') as f:
                        shutil.copyfileobj(f, out)
                    out.write(b'\n')
                    with open(scaffold_path, 'rb') as f:
                        shutil.copyfileobj(f, out)
                n_ref = count_fasta_sequences(ref_fa)
                n_scaffold = count_fasta_sequences(scaffold_path)
                logger.info(f"Pan-genome constructed: {n_ref} reference "
                            f"sequences + {n_scaffold} scaffolds")
                self._write_allele_manifest()
                self._write_stats()
                return combined

        self._write_allele_manifest()
        self._write_stats()
        return scaffold_path

    # ------------------------------------------------------------------
    # Format detection
    # ------------------------------------------------------------------
    def _detect_format(self) -> str:
        """Auto-detect rGFA vs PGGB by peeking at the file.

        rGFA: S lines carry SR:i: tags (checked on first 50 S lines).
        PGGB: S lines lack SR:i: tags. P lines exist but may appear at
              the end of very large files, so we detect by absence of SR
              tags rather than scanning for P lines.
        """
        has_sr = False
        s_checked = 0

        with open(self.gfa_file) as fh:
            for line in fh:
                if line.startswith('S\t'):
                    if 'SR:i:' in line:
                        has_sr = True
                        break
                    s_checked += 1
                    if s_checked >= 50:
                        break
                elif line.startswith('P\t'):
                    return 'PGGB'

        if has_sr:
            return 'rGFA'
        if s_checked > 0:
            # Many S lines without SR:i: tags -> PGGB
            return 'PGGB'

        raise ValueError(
            "Cannot determine GFA format: no S lines or P lines found")

    # ------------------------------------------------------------------
    # rGFA + GAF backend
    # ------------------------------------------------------------------
    def _process_rgfa_gaf(self) -> Tuple[List[AlleleCandidate], int]:
        """Recover rGFA sample paths from a validated primary GAF file.

        rGFA S/L topology alone has no sample/locus semantics.  The default
        alleles mode therefore requires GAF; callers may still request the
        historical topology walker explicitly through ``--gfa-mode legacy``.
        """
        if not self.gfa_gaf:
            raise ValueError(
                'rGFA alleles mode requires --gfa-gaf; use --gfa-mode legacy '
                'only for the historical topology walker')
        if not os.path.isfile(self.gfa_gaf):
            raise FileNotFoundError(f'GAF file not found: {self.gfa_gaf}')
        if self.gfa_gaf_fasta and not os.path.isfile(self.gfa_gaf_fasta):
            raise FileNotFoundError(
                f'GAF query FASTA not found: {self.gfa_gaf_fasta}')

        segments, _links, _coordinates, references, _source_intervals = load_rgfa_graph(self.gfa_file)
        ref_bp = sum(segment.length for segment in segments.values()
                     if segment.rank == 0)
        self._stats['total_segments'] = len(segments)
        self._stats['ref_segments'] = sum(
            segment.rank == 0 for segment in segments.values())
        self._stats['nonref_segments'] = sum(
            segment.rank != 0 for segment in segments.values())
        self._stats['rgfa_gaf'] = self.gfa_gaf
        self._stats['rgfa_gaf_min_mapq'] = self.gfa_gaf_min_mapq
        self._stats['rgfa_gaf_min_identity'] = self.gfa_gaf_min_identity
        self._stats['sequence_source'] = 'graph_complete_nodes'
        self._rgfa_reference_samples = {
            segment.source_name.split('#', 1)[0] for segment in segments.values()
            if segment.rank == 0
        }
        self._stats['rgfa_reference_samples'] = sorted(self._rgfa_reference_samples)

        if not self.reference_genome:
            ref_fa = os.path.join(self.output_dir, 'gfa_reference.fa')
            with open(ref_fa, 'w') as output:
                for contig in sorted(references):
                    reference = references[contig]
                    sequence = ''.join(segments[node_id].sequence
                                       for node_id in reference.node_ids)
                    output.write(f'>{contig}\n')
                    for start in range(0, len(sequence), 80):
                        output.write(sequence[start:start + 80] + '\n')

        store_path = os.path.join(self.output_dir, '.gfa_alleles.sqlite.tmp')
        calls_path = os.path.join(self.output_dir, '.rgfa_gaf_calls.sqlite.tmp')
        self._allele_store = SQLiteAlleleStore(store_path)
        self._rgfa_call_store = RgfaGafCallStore(calls_path)
        self._allele_mode_active = True
        extract_rgfa_gaf_alleles(
            self.gfa_file, self.gfa_gaf,
            min_size=self.min_insertion_size,
            window_size=self.max_insertion_size,
            window_overlap=self.window_overlap,
            min_mapq=self.gfa_gaf_min_mapq,
            min_identity=self.gfa_gaf_min_identity,
            query_fasta=self.gfa_gaf_fasta,
            candidate_sink=self._allele_store.add,
            call_sink=self._rgfa_call_store.add,
            sample_sink=self._rgfa_call_store.register_sample,
            stats=self._stats,
        )
        self._allele_store.finalize()
        self._allele_representatives = self._allele_store.representatives()
        # A candidate's flanks need not be adjacent in another sample.  Scan
        # GAF once more to call every representative locus only when the same
        # complete quality record spans both outer anchors.
        populate_rgfa_gaf_locus_calls(
            self.gfa_file, self.gfa_gaf, self._allele_representatives,
            self._rgfa_call_store.add,
            min_size=self.min_insertion_size,
            min_mapq=self.gfa_gaf_min_mapq,
            min_identity=self.gfa_gaf_min_identity,
            query_fasta=self.gfa_gaf_fasta,
            sample_sink=self._rgfa_call_store.register_sample,
            stats=self._stats,
        )
        self._stats['candidates_raw'] = self._allele_store.occurrence_count
        self._stats['unique_locus_alleles'] = len(self._allele_representatives)
        self._stats['exact_locus_duplicates_collapsed'] = (
            self._allele_store.occurrence_count - len(self._allele_representatives))
        return self._allele_representatives, ref_bp

    # ------------------------------------------------------------------
    # rGFA backend
    # ------------------------------------------------------------------
    def _process_rgfa(self) -> Tuple[List[Tuple[str, str]], int]:
        """Process rGFA: topology-based bubble walking.

        Returns:
            (candidates, ref_total_bp) where candidates is a list of
            (header, sequence) tuples.
        """
        segments = {}       # seg_id -> (seq, sr, sn, so, length)
        fwd_adj = defaultdict(list)   # (seg_id, orient) -> [(seg_id, orient)]
        bwd_adj = defaultdict(list)
        ref_segs = set()
        samples = set()

        logger.info("Parsing rGFA segments and links...")
        with open(self.gfa_file) as fh:
            for line in fh:
                if line.startswith('S\t'):
                    parts = line.rstrip('\n').split('\t')
                    seg_id = parts[1]
                    seq = parts[2]
                    sr = 0
                    sn = ''
                    so = 0
                    ln = len(seq)
                    for tag in parts[3:]:
                        if tag.startswith('SR:i:'):
                            sr = int(tag[5:])
                        elif tag.startswith('SN:Z:'):
                            sn = tag[5:]
                        elif tag.startswith('SO:i:'):
                            so = int(tag[5:])
                        elif tag.startswith('LN:i:'):
                            ln = int(tag[5:])
                    segments[seg_id] = (seq, sr, sn, so, ln)
                    if sr == 0:
                        ref_segs.add(seg_id)
                    if sn:
                        # Extract sample name: "Sample#hap#chr" -> "Sample"
                        sample = sn.split('#')[0] if '#' in sn else sn
                        samples.add(sample)
                elif line.startswith('L\t'):
                    parts = line.rstrip('\n').split('\t')
                    from_id, from_orient = parts[1], parts[2]
                    to_id, to_orient = parts[3], parts[4]
                    fwd_adj[(from_id, from_orient)].append((to_id, to_orient))
                    bwd_adj[(to_id, to_orient)].append((from_id, from_orient))

        total_segs = len(segments)
        n_ref = len(ref_segs)
        n_nonref = total_segs - n_ref
        logger.info(f"rGFA: {total_segs} segments ({n_ref} ref, {n_nonref} non-ref), "
                    f"{sum(len(v) for v in fwd_adj.values())} links, "
                    f"{len(samples)} samples")

        self._stats['total_segments'] = total_segs
        self._stats['ref_segments'] = n_ref
        self._stats['nonref_segments'] = n_nonref
        self._stats['samples'] = sorted(samples)

        # Identify reference sample from SR:i:0 segments
        ref_samples = set()
        for sid in ref_segs:
            sn = segments[sid][2]
            if sn:
                ref_samples.add(sn.split('#')[0] if '#' in sn else sn)
        self._stats['ref_sample'] = sorted(ref_samples)[0] if ref_samples else 'unknown'

        # Reconstruct reference if reference-free mode
        ref_bp = 0
        if not self.reference_genome:
            ref_bp = self._reconstruct_rgfa_reference(segments, ref_segs)
        else:
            ref_bp = sum(seg[4] for sid, seg in segments.items() if sid in ref_segs)

        # Walk non-ref chains
        candidates = self._walk_rgfa_chains(segments, ref_segs, fwd_adj, bwd_adj)

        return candidates, ref_bp

    def _reconstruct_rgfa_reference(self, segments: dict,
                                     ref_segs: set) -> int:
        """Reconstruct reference genome from SR:i:0 segments.

        Groups by SN chromosome, sorts by SO offset, concatenates.
        Returns total reference bp.
        """
        ref_fa = os.path.join(self.output_dir, 'gfa_reference.fa')
        # Group ref segments by chromosome
        by_chr: Dict[str, List[Tuple[int, str, str]]] = defaultdict(list)
        for seg_id in ref_segs:
            seq, sr, sn, so, ln = segments[seg_id]
            if sn:
                by_chr[sn].append((so, seg_id, seq))

        total_bp = 0
        line_width = 80
        with open(ref_fa, 'w') as out:
            for sn in sorted(by_chr.keys()):
                entries = sorted(by_chr[sn], key=lambda x: x[0])
                # Extract short chromosome name from SN (e.g. "Col-0_T2T#1#Chr1" -> "Chr1")
                parts = sn.split('#')
                chr_name = parts[-1] if '#' in sn else sn
                full_seq = ''.join(e[2] for e in entries)
                total_bp += len(full_seq)
                out.write(f'>{chr_name}\n')
                for i in range(0, len(full_seq), line_width):
                    out.write(full_seq[i:i + line_width] + '\n')

        logger.info(f"Reconstructed reference: {len(by_chr)} chromosomes, "
                    f"{total_bp / 1e6:.1f} Mb")
        return total_bp

    def _walk_rgfa_chains(self, segments: dict, ref_segs: set,
                          fwd_adj: dict, bwd_adj: dict
                          ) -> List[Tuple[str, str]]:
        """Walk maximal linear chains of non-reference segments.

        Starting from each unvisited non-ref segment, walk forward and
        backward through non-ref neighbors. At branching points (multiple
        non-ref successors) or ref segments, stop the chain. Sum chain
        length; if in TE range, emit candidate.

        Returns list of (header, sequence) tuples.
        """
        visited: Set[str] = set()
        candidates = []
        candidate_idx = 0

        nonref_ids = [sid for sid in segments if sid not in ref_segs]

        for start_id in nonref_ids:
            if start_id in visited:
                continue

            # Build chain: walk forward and backward from start
            chain_fwd = self._walk_direction(
                start_id, '+', segments, ref_segs, fwd_adj, visited)
            chain_bwd = self._walk_direction(
                start_id, '+', segments, ref_segs, bwd_adj, visited, reverse=True)

            # Combine: bwd (reversed) + start + fwd
            # chain_bwd is in reverse walk order, reverse it for correct order
            chain = list(reversed(chain_bwd)) + [(start_id, '+')] + chain_fwd

            # Mark all as visited
            for seg_id, _ in chain:
                visited.add(seg_id)

            # Calculate total length
            total_len = sum(segments[sid][4] for sid, _ in chain)

            if total_len < self.min_insertion_size or total_len > self.max_insertion_size:
                continue

            # Build sequence respecting orientations
            seq_parts = []
            for seg_id, orient in chain:
                seq = segments[seg_id][0]
                if orient == '-':
                    seq = _reverse_complement(seq)
                seq_parts.append(seq)
            full_seq = ''.join(seq_parts)

            candidate_idx += 1
            header = (f">gfa_ins_{candidate_idx} chain_len={total_len} "
                      f"n_segs={len(chain)} source=rgfa_chain")
            candidates.append((header, full_seq))

        logger.info(f"rGFA chain walking: {candidate_idx} candidates "
                    f"from {len(nonref_ids)} non-ref segments")
        self._stats['candidates_raw'] = candidate_idx
        self._stats['candidates_merged'] = candidate_idx  # no merge step for GFA
        return candidates

    def _walk_direction(self, start_id: str, start_orient: str,
                        segments: dict, ref_segs: set,
                        adj: dict, visited: set,
                        reverse: bool = False) -> List[Tuple[str, str]]:
        """Walk in one direction through non-ref, non-branching neighbors.

        Returns list of (seg_id, orientation) NOT including start_id.
        """
        chain = []
        cur_id, cur_orient = start_id, start_orient

        while True:
            neighbors = adj.get((cur_id, cur_orient), [])
            # Filter to non-ref, unvisited neighbors
            nonref_neighbors = [
                (nid, norient) for nid, norient in neighbors
                if nid not in ref_segs and nid not in visited and nid != start_id
            ]

            if len(nonref_neighbors) != 1:
                break  # Stop at branching or dead end

            next_id, next_orient = nonref_neighbors[0]
            if next_id in visited:
                break

            chain.append((next_id, next_orient))
            visited.add(next_id)
            cur_id, cur_orient = next_id, next_orient

        return chain

    # ------------------------------------------------------------------
    # PGGB backend
    # ------------------------------------------------------------------
    def _process_pggb(self) -> Tuple[List[Tuple[str, str]], int]:
        """Process PGGB paths with a streaming locus-aware backend."""
        logger.info("Parsing PGGB segments...")
        seg_seqs = {}
        with open(self.gfa_file) as handle:
            for line_number, line in enumerate(handle, 1):
                if not line.startswith('S\t'):
                    continue
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 3:
                    raise ValueError(f"Malformed S line {line_number}")
                seg_id = sys.intern(parts[1])
                sequence = parts[2]
                if not sequence or sequence == '*':
                    raise ValueError(
                        f"PGGB segment {seg_id!r} has no embedded sequence")
                seg_seqs[seg_id] = sequence
        total_segs = len(seg_seqs)
        logger.info(f"PGGB: {total_segs} segments loaded")

        # First path pass retains names only.  Hundreds of cohort paths may
        # contain billions of steps, so alternate steps are streamed later.
        logger.info("Indexing PGGB path names...")
        path_names = []
        with open(self.gfa_file) as handle:
            for line_number, line in enumerate(handle, 1):
                if not line.startswith('P\t'):
                    continue
                parts = line.rstrip('\n').split('\t', 2)
                if len(parts) < 3:
                    raise ValueError(f"Malformed P line {line_number}")
                path_names.append(parts[1])
        if not path_names:
            raise ValueError("PGGB GFA has no P lines")
        if len(path_names) != len(set(path_names)):
            raise ValueError("PGGB GFA contains duplicate P path names")

        logger.info(f"PGGB: {len(path_names)} paths")
        self._stats['samples'] = sorted(path_names)
        self._stats['biological_samples'] = sorted({
            self._pggb_path_group(name) for name in path_names})
        path_catalog = {name: None for name in path_names}
        ref_path_names = self._select_pggb_references(path_catalog)
        ref_path_set = set(ref_path_names)
        ref_sample = self._pggb_path_group(ref_path_names[0])
        logger.info(f"PGGB reference sample: {ref_sample} "
                    f"({len(ref_path_names)} paths)")
        for name in ref_path_names:
            logger.info(f"  reference path: {name}")
        self._stats['ref_sample'] = ref_sample
        self._stats['ref_paths'] = ref_path_names

        # Second pass materializes only the selected backbone paths.
        ref_paths = {}
        legacy_paths = {} if self.gfa_mode == 'legacy' else None
        with open(self.gfa_file) as handle:
            for line_number, line in enumerate(handle, 1):
                if not line.startswith('P\t'):
                    continue
                parts = line.rstrip('\n').split('\t')
                name = parts[1]
                if name in ref_path_set or legacy_paths is not None:
                    steps = parse_p_steps(parts)
                    if name in ref_path_set:
                        ref_paths[name] = steps
                    if legacy_paths is not None:
                        legacy_paths[name] = steps
        missing_refs = [name for name in ref_path_names if name not in ref_paths]
        if missing_refs:
            raise ValueError(f"Missing selected reference paths: {missing_refs}")

        ref_seg_set: Set[str] = {
            seg_id for steps in ref_paths.values() for seg_id, _ in steps}
        self._stats['total_segments'] = total_segs
        self._stats['ref_segments'] = len(ref_seg_set)
        self._stats['nonref_segments'] = total_segs - len(ref_seg_set)

        if not self.reference_genome:
            ref_bp = self._reconstruct_pggb_references(
                ref_path_names, ref_paths, seg_seqs)
        else:
            ref_bp = sum(
                len(seg_seqs[seg_id])
                for steps in ref_paths.values() for seg_id, _ in steps)

        if self.gfa_mode == 'alleles':
            self._allele_mode_active = True
            store_path = os.path.join(
                self.output_dir, '.gfa_alleles.sqlite.tmp')
            self._allele_store = SQLiteAlleleStore(store_path)
            try:
                extract_pggb_alleles_stream(
                    gfa_file=self.gfa_file,
                    ref_paths=ref_paths,
                    sequences=seg_seqs,
                    min_size=self.min_insertion_size,
                    window_size=self.max_insertion_size,
                    window_overlap=self.window_overlap,
                    stats=self._stats,
                    candidate_sink=self._allele_store.add,
                )
                self._allele_store.finalize()
                # Representative sequences are now durable in SQLite.  Drop
                # the embedded graph sequence dictionary before materializing
                # them so graph and allele sequence payloads do not overlap.
                seg_seqs.clear()
                gc.collect()
                self._allele_representatives = (
                    self._allele_store.representatives())
            except Exception:
                self._allele_store.close(remove=False)
                raise
            candidates = self._allele_representatives
            occurrence_count = self._allele_store.occurrence_count
            representative_count = self._allele_store.representative_count
            self._stats['candidates_raw'] = occurrence_count
            self._stats['unique_locus_alleles'] = representative_count
            self._stats['exact_locus_duplicates_collapsed'] = (
                occurrence_count - representative_count)
            logger.info(
                f"PGGB locus-aware extraction: {occurrence_count} candidate "
                f"occurrences, {representative_count} unique locus alleles "
                f"from {len(path_names) - len(ref_path_names)} paths")
        else:
            candidates = self._extract_pggb_chains(
                legacy_paths, ref_path_set, ref_seg_set, seg_seqs)
        return candidates, ref_bp

    @staticmethod
    def _pggb_path_group(path_name: str) -> str:
        """Return the PGGB sample/haplotype prefix for a path name.

        PGGB/wfmash path names conventionally use sample#haplotype#contig.
        Grouping on the first two fields keeps all chromosomes from one
        haplotype together while retaining sensible behavior for simple names.
        """
        parts = path_name.split('#')
        if len(parts) >= 3:
            return '#'.join(parts[:2])
        return path_name

    def _select_pggb_references(self, paths: dict) -> List[str]:
        """Select every path from one PGGB reference sample/haplotype.

        Priority for the anchor path is explicit ``--gfa-ref``, then a T2T
        name, then the first path.  All paths sharing the anchor's
        sample#haplotype prefix are returned in input order.
        """
        anchor = None
        if self.gfa_ref:
            if self.gfa_ref in paths:
                anchor = self.gfa_ref
            else:
                # Permit an exact sample/haplotype prefix such as Col-0_T2T#1.
                for name in paths:
                    if self._pggb_path_group(name) == self.gfa_ref:
                        anchor = name
                        break
                if anchor is None:
                    for name in paths:
                        if self.gfa_ref in name:
                            anchor = name
                            break
            if anchor is None:
                raise ValueError(
                    f"--gfa-ref '{self.gfa_ref}' not found in paths: "
                    f"{list(paths.keys())}")
        else:
            for name in paths:
                if 't2t' in name.lower():
                    anchor = name
                    break
            if anchor is None:
                anchor = next(iter(paths))

        group = self._pggb_path_group(anchor)
        selected = [name for name in paths
                    if self._pggb_path_group(name) == group]
        return selected

    def _reconstruct_pggb_references(
            self, path_names: List[str], paths: dict,
            seg_seqs: dict) -> int:
        """Reconstruct all selected reference paths as separate FASTA records."""
        ref_fa = os.path.join(self.output_dir, 'gfa_reference.fa')
        line_width = 80
        total_bp = 0
        used_headers = set()

        with open(ref_fa, 'w') as out:
            for path_name in path_names:
                parts = path_name.split('#')
                header = parts[-1] if len(parts) >= 3 else path_name
                if header in used_headers:
                    header = path_name.replace('#', '_')
                used_headers.add(header)

                seq_parts = []
                for seg_id, orient in paths[path_name]:
                    seq = seg_seqs.get(seg_id, '')
                    if orient == '-':
                        seq = _reverse_complement(seq)
                    seq_parts.append(seq)
                full_seq = ''.join(seq_parts)
                total_bp += len(full_seq)

                out.write(f'>{header}\n')
                for i in range(0, len(full_seq), line_width):
                    out.write(full_seq[i:i + line_width] + '\n')
                logger.info(f"Reconstructed PGGB reference: {header}, "
                            f"{len(full_seq) / 1e6:.1f} Mb")

        logger.info(f"Reconstructed {len(path_names)} PGGB reference paths, "
                    f"{total_bp / 1e6:.1f} Mb total")
        return total_bp

    def _extract_pggb_chains(self, paths: dict,
                              ref_path_names: Set[str],
                              ref_seg_set: set, seg_seqs: dict
                              ) -> List[Tuple[str, str]]:
        """Walk non-ref paths, extracting chains of non-ref segments.

        At each ref segment, flush the current chain. At each non-ref
        segment, extend the current chain. Chains in TE size range are
        emitted.
        """
        candidates = []
        candidate_idx = 0

        for path_name, steps in paths.items():
            if path_name in ref_path_names:
                continue

            chain_segs = []  # [(seg_id, orient)]
            chain_len = 0

            for seg_id, orient in steps:
                if seg_id in ref_seg_set:
                    # Flush current chain
                    if chain_len >= self.min_insertion_size and chain_len <= self.max_insertion_size:
                        candidate_idx += 1
                        seq = self._build_chain_seq(chain_segs, seg_seqs)
                        header = (f">gfa_ins_{candidate_idx} "
                                  f"chain_len={chain_len} "
                                  f"n_segs={len(chain_segs)} "
                                  f"path={path_name} source=pggb_chain")
                        candidates.append((header, seq))
                    chain_segs = []
                    chain_len = 0
                else:
                    seg_len = len(seg_seqs.get(seg_id, ''))
                    new_len = chain_len + seg_len
                    if new_len > self.max_insertion_size:
                        # Abort chain, reset
                        chain_segs = []
                        chain_len = 0
                    else:
                        chain_segs.append((seg_id, orient))
                        chain_len = new_len

            # Flush trailing chain
            if chain_len >= self.min_insertion_size and chain_len <= self.max_insertion_size:
                candidate_idx += 1
                seq = self._build_chain_seq(chain_segs, seg_seqs)
                header = (f">gfa_ins_{candidate_idx} "
                          f"chain_len={chain_len} "
                          f"n_segs={len(chain_segs)} "
                          f"path={path_name} source=pggb_chain")
                candidates.append((header, seq))

        nonref_path_count = len(paths) - len(ref_path_names)
        logger.info(f"PGGB chain extraction: {candidate_idx} candidates "
                    f"from {nonref_path_count} non-ref paths")
        self._stats['candidates_raw'] = candidate_idx
        self._stats['candidates_merged'] = candidate_idx  # no merge step for GFA
        return candidates

    def _build_chain_seq(self, chain_segs: List[Tuple[str, str]],
                         seg_seqs: dict) -> str:
        """Concatenate chain segment sequences respecting orientation."""
        parts = []
        for seg_id, orient in chain_segs:
            seq = seg_seqs.get(seg_id, '')
            if orient == '-':
                seq = _reverse_complement(seq)
            parts.append(seq)
        return ''.join(parts)

    def _cleanup_candidate_intermediates(self) -> None:
        if self.keep_intermediates:
            return
        for tmp_name in [
                'gfa_raw_insertions.fa', 'gfa_unique_locus_alleles.fa',
                'gfa_filtered_insertions.fa', 'gfa_te_insertions.fa']:
            tmp_path = os.path.join(self.output_dir, tmp_name)
            if os.path.exists(tmp_path):
                os.remove(tmp_path)

    def _write_allele_candidates(
            self, records: List[AlleleCandidate], output_path: str) -> str:
        """Write representative allele candidates with support metadata."""
        with open(output_path, 'w') as out:
            for record in records:
                out.write(record.fasta_header() + '\n')
                sequence = record.sequence
                for start in range(0, len(sequence), 80):
                    out.write(sequence[start:start + 80] + '\n')
        logger.info(
            f"Wrote {len(records)} unique locus-aware allele candidates")
        return output_path

    def _mark_allele_qc(self, filtered_fasta: str) -> None:
        """Propagate representative QC status to every allele occurrence."""
        passed = set()
        with open(filtered_fasta) as handle:
            for line in handle:
                if line.startswith('>'):
                    passed.add(line[1:].split()[0])
        if self._allele_store is not None:
            self._allele_store.set_qc(passed)
            return
        for occurrence in self._allele_occurrences:
            occurrence.qc_status = (
                'passed' if occurrence.representative_id in passed
                else 'filtered')

    def _write_allele_manifest(self) -> None:
        """Write occurrence provenance for locus-aware PGGB or rGFA+GAF alleles."""
        if not self._allele_mode_active:
            return
        manifest = os.path.join(self.output_dir, 'gfa_allele_manifest.tsv')
        if self._allele_store is not None:
            self._allele_store.write_manifest(manifest)
            self._allele_store.close(remove=True)
            self._allele_store = None
        else:
            write_manifest(manifest, self._allele_occurrences)
        self._stats['allele_manifest'] = manifest
        if self._rgfa_call_store is not None:
            presence = os.path.join(self.output_dir, 'gfa_rgfa_locus_presence.tsv')
            samples = self._rgfa_call_store.write_presence_matrix(
                presence, (record.locus_key for record in self._allele_representatives),
                reference_samples=self._rgfa_reference_samples)
            self._rgfa_call_store.close(remove=True)
            self._rgfa_call_store = None
            self._stats['rgfa_locus_presence'] = presence
            self._stats['rgfa_gaf_samples'] = samples

    def _cleanup_open_allele_stores(self) -> None:
        """Close/remove temporary stores after any failed GFA run."""
        for attribute in ('_allele_store', '_rgfa_call_store'):
            store = getattr(self, attribute, None)
            if store is None:
                continue
            try:
                store.close(remove=True)
            finally:
                setattr(self, attribute, None)

    # ------------------------------------------------------------------
    # Candidate output
    # ------------------------------------------------------------------
    def _write_raw_candidates(self, candidates: List[Tuple[str, str]]) -> str:
        """Write raw candidate sequences to FASTA."""
        out_path = os.path.join(self.output_dir, 'gfa_raw_insertions.fa')
        with open(out_path, 'w') as out:
            for header, seq in candidates:
                out.write(header + '\n')
                out.write(seq + '\n')
        n = len(candidates)
        logger.info(f"Wrote {n} raw candidate sequences")
        self._stats['extracted_sequences'] = n
        return out_path

    # ------------------------------------------------------------------
    # QC filter
    # ------------------------------------------------------------------
    def _qc_filter(self, fasta_path: str) -> str:
        """Filter by entropy, DUST, N-content, length."""
        out_path = os.path.join(self.output_dir, 'gfa_filtered_insertions.fa')
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
    # CD-HIT-EST dedup
    # ------------------------------------------------------------------
    def _cdhit_dedup(self, fasta_path: str) -> str:
        """Deduplicate with cd-hit-est."""
        n_input = count_fasta_sequences(fasta_path)
        if n_input <= 1:
            final = os.path.join(self.output_dir, 'gfa_te_insertions.fa')
            shutil.copy2(fasta_path, final)
            return final

        out_path = os.path.join(self.output_dir, 'gfa_te_insertions.fa')
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

        clstr = out_path + '.clstr'
        if os.path.exists(clstr):
            os.remove(clstr)

        return out_path

    # ------------------------------------------------------------------
    # Scaffold construction
    # ------------------------------------------------------------------
    def _build_synthetic_scaffolds(self, fasta_path: str, output_path: str):
        """Pack candidates, or emit independent records when spacer length is 0."""
        if self.n_spacer_length == 0:
            shutil.copyfile(fasta_path, output_path)
            count = count_fasta_sequences(output_path)
            total_bp = 0
            with open(output_path) as handle:
                for line in handle:
                    if not line.startswith('>'):
                        total_bp += len(line.strip())
            self._stats['scaffold_count'] = count
            self._stats['scaffold_total_bp'] = total_bp
            self._stats['artificial_spacer_bp'] = 0
            logger.info(
                f"Wrote {count} independent candidate contigs, "
                f"{total_bp / 1e6:.1f} Mb total")
            return

        spacer = 'N' * self.n_spacer_length
        line_width = 80

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

        with open(output_path, 'w') as out:
            for i, (scaffold_seq, n_parts) in enumerate(scaffolds, 1):
                out.write(f">gfa_scaffold_{i} "
                          f"n_insertions={n_parts} "
                          f"len={len(scaffold_seq)}\n")
                for j in range(0, len(scaffold_seq), line_width):
                    out.write(scaffold_seq[j:j + line_width] + '\n')

        total_bp = sum(len(seq) for seq, _ in scaffolds)
        self._stats['scaffold_count'] = len(scaffolds)
        self._stats['scaffold_total_bp'] = total_bp
        self._stats['artificial_spacer_bp'] = max(
            0, (len(sequences) - len(scaffolds)) * self.n_spacer_length)
        logger.info(f"Built {len(scaffolds)} synthetic scaffold(s), "
                    f"total {total_bp / 1e6:.1f} Mb")

    # ------------------------------------------------------------------
    # Utilities
    # ------------------------------------------------------------------
    def _log_config(self):
        """Log current configuration for reproducibility."""
        logger.info(
            f"GFA Processor v1.0 | min_ins={self.min_insertion_size} "
            f"max_ins={self.max_insertion_size} "
            f"entropy={self.entropy_threshold_long}/{self.entropy_threshold_short} "
            f"dust={self.dust_threshold} gfa_ref={self.gfa_ref} "
            f"gfa_mode={self.gfa_mode} overlap={self.window_overlap} "
            f"gaf={self.gfa_gaf} gaf_fasta={self.gfa_gaf_fasta} "
            f"gaf_min_mapq={self.gfa_gaf_min_mapq} "
            f"gaf_min_identity={self.gfa_gaf_min_identity} "
            f"spacer={self.n_spacer_length} keep_intermediates="
            f"{self.keep_intermediates}")

    def _write_stats(self):
        """Write processing statistics to JSON."""
        stats = dict(self._stats)
        stats['gfa_format'] = self.gfa_format
        stats_path = os.path.join(self.output_dir, 'gfa_processing_stats.json')
        with open(stats_path, 'w') as f:
            json.dump(stats, f, indent=2)


# ------------------------------------------------------------------
# Standalone CLI
# ------------------------------------------------------------------
def main():
    """Standalone CLI for GFA TE insertion extraction."""
    import argparse
    parser = argparse.ArgumentParser(
        description='Extract TE insertion sequences from pangenome GFA graphs')
    parser.add_argument('--gfa', required=True, help='GFA file (rGFA or PGGB)')
    parser.add_argument('--ref', help='Reference genome FASTA (enrichment mode)')
    parser.add_argument('--gfa-ref',
                        help='Explicit reference path name for PGGB')
    parser.add_argument('--gfa-mode', choices=('alleles', 'legacy'),
                        default='alleles')
    parser.add_argument('--window-overlap', type=int, default=1000)
    parser.add_argument('--gfa-gaf', help='Primary GAF used to recover rGFA sample paths')
    parser.add_argument('--gfa-gaf-fasta', help='Query FASTA used to validate GAF names and lengths')
    parser.add_argument('--gfa-gaf-min-mapq', type=int, default=20)
    parser.add_argument('--gfa-gaf-min-identity', type=float, default=0.90)
    parser.add_argument('--spacer-length', type=int, default=0)
    parser.add_argument('--clean-intermediates', action='store_true')
    parser.add_argument('--outdir', required=True, help='Output directory')
    parser.add_argument('--threads', type=int, default=4)
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO,
                        format='[%(asctime)s] %(levelname)s: %(message)s')

    proc = GfaProcessor(
        gfa_file=args.gfa,
        output_dir=args.outdir,
        reference_genome=args.ref,
        threads=args.threads,
        gfa_ref=args.gfa_ref,
        gfa_mode=args.gfa_mode,
        window_overlap=args.window_overlap,
        spacer_length=args.spacer_length,
        keep_intermediates=not args.clean_intermediates,
        gfa_gaf=args.gfa_gaf,
        gfa_gaf_fasta=args.gfa_gaf_fasta,
        gfa_gaf_min_mapq=args.gfa_gaf_min_mapq,
        gfa_gaf_min_identity=args.gfa_gaf_min_identity
    )
    result = proc.process()
    print(f"Output: {result}")


if __name__ == '__main__':
    main()
