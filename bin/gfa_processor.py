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

    For PGGB (smoothxg), a reference path is selected and all segments
    NOT in that path are considered non-reference. Consecutive non-ref
    segments along each non-ref path form chains; chains in the TE size
    range are emitted.

    In both cases, extracted TE candidates undergo QC filtering (entropy,
    DUST, N-content), CD-HIT-EST deduplication, and scaffold construction.
    """

    def __init__(self, gfa_file: str, output_dir: str,
                 reference_genome: str = None, threads: int = 4,
                 gfa_ref: str = None):
        """
        Args:
            gfa_file: Path to GFA file (rGFA or PGGB).
            output_dir: Working directory for intermediate and output files.
            reference_genome: If provided, enrichment mode; None means
                reference-free mode.
            threads: Thread count for cd-hit-est.
            gfa_ref: Optional explicit reference path/sample name for PGGB.
        """
        self.gfa_file = gfa_file
        self.output_dir = output_dir
        self.reference_genome = reference_genome
        self.threads = threads
        self.gfa_ref = gfa_ref

        # Size parameters
        self.min_insertion_size = 100
        self.max_insertion_size = 30000

        # QC parameters (same as PAFProcessor)
        self.entropy_threshold_long = 1.0
        self.entropy_threshold_short = 0.8
        self.dust_threshold = 0.7
        self.max_n_percent = 0.20
        self.cdhit_identity = 0.8

        # Scaffold parameters
        self.n_spacer_length = 5000
        self.max_scaffold_length = 50000000

        # Runtime state
        self.gfa_format: str = 'unknown'
        self._stats: Dict = defaultdict(int)

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
        self._log_config()
        self.gfa_format = self._detect_format()
        logger.info(f"GFA format detected: {self.gfa_format}")

        if self.gfa_format == 'rGFA':
            candidates, ref_bp = self._process_rgfa()
        else:
            candidates, ref_bp = self._process_pggb()

        self._stats['ref_total_bp'] = ref_bp

        scaffold_path = os.path.join(self.output_dir, 'gfa_scaffolds.fa')

        if not candidates:
            logger.warning("No TE insertion candidates found in GFA")
            with open(scaffold_path, 'w'):
                pass
            self._write_stats()
            if not self.reference_genome:
                ref_fa = os.path.join(self.output_dir, 'gfa_reference.fa')
                if os.path.exists(ref_fa):
                    return ref_fa
            return scaffold_path

        raw_fasta = self._write_raw_candidates(candidates)
        filtered_fasta = self._qc_filter(raw_fasta)
        deduped_fasta = self._cdhit_dedup(filtered_fasta)

        n_seqs = count_fasta_sequences(deduped_fasta)
        self._stats['final_sequences'] = n_seqs

        if n_seqs == 0:
            with open(scaffold_path, 'w'):
                pass
            self._write_stats()
            if not self.reference_genome:
                ref_fa = os.path.join(self.output_dir, 'gfa_reference.fa')
                if os.path.exists(ref_fa):
                    return ref_fa
            return scaffold_path

        self._build_synthetic_scaffolds(deduped_fasta, scaffold_path)
        logger.info(f"GFA enrichment complete: {n_seqs} TE insertions "
                    f"assembled into scaffolds")

        # Cleanup intermediate files
        for tmp_name in ['gfa_raw_insertions.fa', 'gfa_filtered_insertions.fa',
                         'gfa_te_insertions.fa']:
            tmp_path = os.path.join(self.output_dir, tmp_name)
            if os.path.exists(tmp_path):
                os.remove(tmp_path)

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
                self._write_stats()
                return combined

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
        """Process PGGB: path-based chain extraction.

        Returns:
            (candidates, ref_total_bp)
        """
        logger.info("Parsing PGGB segments...")
        seg_seqs = {}  # seg_id -> sequence

        with open(self.gfa_file) as fh:
            for line in fh:
                if line.startswith('S\t'):
                    parts = line.rstrip('\n').split('\t')
                    seg_seqs[parts[1]] = parts[2]

        total_segs = len(seg_seqs)
        logger.info(f"PGGB: {total_segs} segments loaded")

        # Parse P lines
        logger.info("Parsing PGGB paths...")
        paths: Dict[str, List[Tuple[str, str]]] = {}  # name -> [(seg_id, orient)]
        with open(self.gfa_file) as fh:
            for line in fh:
                if line.startswith('P\t'):
                    parts = line.rstrip('\n').split('\t')
                    path_name = parts[1]
                    steps = []
                    for step in parts[2].split(','):
                        seg_id = step[:-1]
                        orient = step[-1]
                        steps.append((seg_id, orient))
                    paths[path_name] = steps

        if not paths:
            raise ValueError("PGGB GFA has no P lines")

        logger.info(f"PGGB: {len(paths)} paths")
        samples = sorted(paths.keys())
        self._stats['samples'] = samples

        # Select reference path
        ref_path_name = self._select_pggb_reference(paths)
        logger.info(f"PGGB reference path: {ref_path_name}")
        self._stats['ref_sample'] = ref_path_name

        ref_steps = paths[ref_path_name]
        ref_seg_set: Set[str] = set()
        for seg_id, _ in ref_steps:
            ref_seg_set.add(seg_id)

        n_ref = len(ref_seg_set)
        n_nonref = total_segs - n_ref
        self._stats['total_segments'] = total_segs
        self._stats['ref_segments'] = n_ref
        self._stats['nonref_segments'] = n_nonref

        # Reconstruct reference if reference-free mode
        ref_bp = 0
        if not self.reference_genome:
            ref_bp = self._reconstruct_pggb_reference(
                ref_path_name, ref_steps, seg_seqs)
        else:
            ref_bp = sum(len(seg_seqs[sid]) for sid, _ in ref_steps
                         if sid in seg_seqs)

        # Extract candidates from non-ref paths
        candidates = self._extract_pggb_chains(
            paths, ref_path_name, ref_seg_set, seg_seqs)

        return candidates, ref_bp

    def _select_pggb_reference(self, paths: dict) -> str:
        """Select the reference path for PGGB.

        Priority: explicit --gfa-ref > T2T pattern match > first path.
        """
        if self.gfa_ref:
            # Exact match first
            if self.gfa_ref in paths:
                return self.gfa_ref
            # Substring match
            for name in paths:
                if self.gfa_ref in name:
                    return name
            raise ValueError(
                f"--gfa-ref '{self.gfa_ref}' not found in paths: "
                f"{list(paths.keys())}")

        # Auto-detect: prefer T2T
        for name in paths:
            if 'T2T' in name or 't2t' in name:
                return name

        # Fall back to first path
        return next(iter(paths))

    def _reconstruct_pggb_reference(self, path_name: str,
                                     steps: List[Tuple[str, str]],
                                     seg_seqs: dict) -> int:
        """Reconstruct reference by walking the reference P line."""
        ref_fa = os.path.join(self.output_dir, 'gfa_reference.fa')
        line_width = 80

        # Extract chromosome name from path name
        parts = path_name.split('#')
        chr_name = parts[-1] if '#' in path_name else path_name

        # Build sequence
        seq_parts = []
        for seg_id, orient in steps:
            seq = seg_seqs.get(seg_id, '')
            if orient == '-':
                seq = _reverse_complement(seq)
            seq_parts.append(seq)
        full_seq = ''.join(seq_parts)
        total_bp = len(full_seq)

        with open(ref_fa, 'w') as out:
            out.write(f'>{chr_name}\n')
            for i in range(0, len(full_seq), line_width):
                out.write(full_seq[i:i + line_width] + '\n')

        logger.info(f"Reconstructed PGGB reference: {chr_name}, "
                    f"{total_bp / 1e6:.1f} Mb")
        return total_bp

    def _extract_pggb_chains(self, paths: dict, ref_path_name: str,
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
            if path_name == ref_path_name:
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

        logger.info(f"PGGB chain extraction: {candidate_idx} candidates "
                    f"from {len(paths) - 1} non-ref paths")
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
        """Concatenate TE insertions with N-spacers into synthetic scaffolds."""
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
            f"dust={self.dust_threshold} gfa_ref={self.gfa_ref}")

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
        gfa_ref=args.gfa_ref
    )
    result = proc.process()
    print(f"Output: {result}")


if __name__ == '__main__':
    main()
