import logging
import shutil
from typing import Dict, List, Tuple, Optional
from Bio import SeqIO
from Bio.Seq import Seq
import subprocess
import tempfile
import os

logger = logging.getLogger(__name__)

# Module-level cache for pysam.FastaFile handles (one per genome path)
_fasta_handles = {}

def _get_fasta_handle(genome_file: str):
    """Get or create a pysam.FastaFile handle for indexed O(1) extraction.
    Returns None if pysam is unavailable or the file can't be opened."""
    global _fasta_handles
    if genome_file in _fasta_handles:
        return _fasta_handles[genome_file]
    try:
        import pysam
        fai_path = genome_file + '.fai'
        if not os.path.exists(fai_path):
            # Try to create index
            subprocess.run(['samtools', 'faidx', genome_file],
                           capture_output=True, check=True, timeout=120)
        if os.path.exists(fai_path):
            handle = pysam.FastaFile(genome_file)
            _fasta_handles[genome_file] = handle
            logger.info(f"Opened pysam indexed FASTA: {genome_file}")
            return handle
    except ImportError:
        logger.debug("pysam not available, will try samtools CLI")
    except Exception as e:
        logger.debug(f"pysam open failed: {e}")
    _fasta_handles[genome_file] = None
    return None

def _has_samtools() -> bool:
    """Check if samtools is available on PATH."""
    return shutil.which('samtools') is not None

def _samtools_faidx_extract(genome_file: str, chrom: str, start: int, end: int) -> str:
    """Extract sequence using samtools faidx (0-based start, 0-based exclusive end).
    samtools uses 1-based inclusive coordinates, so convert."""
    region = f"{chrom}:{start+1}-{end}"
    try:
        result = subprocess.run(
            ['samtools', 'faidx', genome_file, region],
            capture_output=True, text=True, timeout=30
        )
        if result.returncode != 0:
            return ""
        lines = []
        for line in result.stdout.split('\n'):
            if not line.startswith('>') and line.strip():
                lines.append(line.strip())
        return ''.join(lines)
    except (subprocess.TimeoutExpired, subprocess.SubprocessError, FileNotFoundError):
        return ""

def reverse_complement(seq: str) -> str:
    """获取序列的反向互补"""
    return str(Seq(seq).reverse_complement())

def calculate_identity(seq1: str, seq2: str) -> float:
    """计算两个序列的相似度（简单版本）"""
    if not seq1 or not seq2:
        return 0.0
    
    # 使用较短序列的长度作为基准
    min_len = min(len(seq1), len(seq2))
    matches = sum(1 for i in range(min_len) if seq1[i] == seq2[i])
    
    return matches / min_len if min_len > 0 else 0.0

def calculate_sequence_identity(seq1: str, seq2: str) -> float:
    """计算两个序列的相似度（与calculate_identity相同的实现）"""
    return calculate_identity(seq1, seq2)

def extract_sequence_from_genome(genome_file: str, chrom: str,
                               start: int, end: int, strand: str = '+') -> str:
    """Extract a genomic region using indexed access (O(1) per call).

    Uses pysam if available, otherwise falls back to samtools faidx CLI,
    and finally to BioPython sequential scan as last resort.
    Coordinates are 0-based, half-open [start, end).
    """
    logger.debug(f"Extracting sequence from {genome_file}: {chrom}:{start}-{end} ({strand})")

    if not os.path.exists(genome_file):
        logger.error(f"Genome file not found: {genome_file}")
        return ""

    start = max(0, start)
    if start >= end:
        logger.warning(f"Invalid coordinates: {chrom}:{start}-{end}")
        return ""

    extracted_seq = ""

    # Strategy 1: pysam (fastest — memory-mapped indexed access)
    handle = _get_fasta_handle(genome_file)
    if handle is not None:
        try:
            chrom_len = handle.get_reference_length(chrom)
            adj_end = min(end, chrom_len)
            if start >= adj_end:
                logger.warning(f"Coordinates out of range: {chrom}:{start}-{end} (chrom_len={chrom_len})")
                return ""
            extracted_seq = handle.fetch(chrom, start, adj_end)
        except (KeyError, ValueError):
            # Chromosome not found in index
            logger.debug(f"pysam: chromosome '{chrom}' not in index, trying samtools CLI")
            extracted_seq = ""

    # Strategy 2: samtools faidx CLI
    if not extracted_seq and _has_samtools() and os.path.exists(genome_file + '.fai'):
        extracted_seq = _samtools_faidx_extract(genome_file, chrom, start, end)

    # Strategy 3: BioPython sequential scan (slowest fallback)
    if not extracted_seq:
        logger.debug(f"Falling back to BioPython for {chrom}:{start}-{end}")
        try:
            with open(genome_file, 'r') as f:
                for record in SeqIO.parse(f, 'fasta'):
                    if record.id == chrom:
                        genome_seq = str(record.seq)
                        adj_end = min(end, len(genome_seq))
                        if start < adj_end:
                            extracted_seq = genome_seq[start:adj_end]
                        break
        except Exception as e:
            logger.error(f"BioPython extraction failed: {e}")
            return ""

    if not extracted_seq:
        logger.warning(f"Chromosome '{chrom}' not found in genome file")
        return ""

    if strand == '-':
        extracted_seq = reverse_complement(extracted_seq)

    return extracted_seq

def extract_sequence_with_flanking(genome_file: str, chrom: str,
                                  start: int, end: int,
                                  flanking: int = 50) -> Dict:
    """Extract target sequence plus flanking regions using indexed access."""
    extract_start = max(0, start - flanking)
    extract_end = end + flanking

    full_seq = extract_sequence_from_genome(genome_file, chrom, extract_start, extract_end)
    if not full_seq:
        return {'sequence': '', 'left_flank': '', 'right_flank': '', 'full_sequence': ''}

    # Calculate offsets within extracted sequence
    left_len = start - extract_start
    right_start = left_len + (end - start)

    left_flank = full_seq[:left_len]
    target_seq = full_seq[left_len:right_start]
    right_flank = full_seq[right_start:]

    return {
        'sequence': target_seq,
        'left_flank': left_flank,
        'right_flank': right_flank,
        'full_sequence': full_seq
    }

def find_best_genome_match(consensus_seq: str, genome_file: str, config) -> Dict:
    """找到共识序列在基因组中的最佳匹配位置"""
    # 创建临时文件
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as query_file:
        query_file.write(f">query\n{consensus_seq}\n")
        query_file_path = query_file.name
    
    try:
        # 运行BLAST搜索
        cmd = [
            config.blastn_exe,
            '-query', query_file_path,
            '-subject', genome_file,
            '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore',
            '-max_target_seqs', '1',
            '-max_hsps', '1'
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        if result.stdout.strip():
            # 解析BLAST结果
            fields = result.stdout.strip().split('\t')
            return {
                'chrom': fields[1],
                'identity': float(fields[2]),
                'length': int(fields[3]),
                'start': int(fields[8]),
                'end': int(fields[9]),
                'evalue': float(fields[10]),
                'score': float(fields[11])
            }
        else:
            logger.warning("No BLAST matches found")
            return {}
            
    except subprocess.CalledProcessError as e:
        logger.error(f"BLAST search failed: {e}")
        return {}
    finally:
        # 清理临时文件
        if os.path.exists(query_file_path):
            os.unlink(query_file_path)

def detect_tsd_from_sequence(sequence: str, flanking_size: int = 30, 
                           min_length: int = 4, max_length: int = 20, 
                           max_mismatches: int = 1) -> Optional[Dict]:
    """从序列中检测目标位点重复（TSD）"""
    if not sequence or len(sequence) < 2 * flanking_size + min_length:
        return None
    
    # 从序列两端提取侧翼序列
    left_flank = sequence[:flanking_size]
    right_flank = sequence[-flanking_size:]
    
    best_tsd = None
    best_score = 0
    
    # 搜索可能的TSD
    for tsd_len in range(min_length, min(max_length + 1, len(left_flank) + 1)):
        # 获取左侧潜在TSD（从左侧侧翼的末尾）
        left_tsd = left_flank[-tsd_len:] if len(left_flank) >= tsd_len else ""
        
        if not left_tsd:
            continue
        
        # 在右侧搜索匹配（从右侧侧翼的开头）
        for i in range(min(tsd_len, len(right_flank))):
            right_tsd = right_flank[i:i+tsd_len]
            
            if len(right_tsd) != tsd_len:
                continue
            
            # 计算匹配度
            mismatches = sum(1 for j in range(tsd_len) if left_tsd[j] != right_tsd[j])
            
            if mismatches <= max_mismatches:
                # 计算得分（长度和匹配度的组合）
                score = tsd_len * (1 - mismatches / tsd_len)
                
                if score > best_score:
                    best_score = score
                    best_tsd = {
                        'sequence': left_tsd,
                        'length': tsd_len,
                        'mismatches': mismatches,
                        'type': 'perfect' if mismatches == 0 else 'imperfect',
                        'confidence': score / max_length,
                        'left_pos': flanking_size - tsd_len,
                        'right_pos': len(sequence) - flanking_size + i,
                        'start_in_seq': flanking_size - tsd_len,
                        'end_in_seq': len(sequence) - flanking_size + i + tsd_len
                    }
    
    return best_tsd

def detect_tsd(genome_file: str, chrom: str, start: int, end: int,
              search_range: int = 30, min_length: int = 4, 
              max_length: int = 20, max_mismatches: int = 1) -> Optional[Dict]:
    """检测目标位点重复（TSD）"""
    # 提取目标序列及侧翼
    extracted = extract_sequence_with_flanking(
        genome_file, chrom, start, end, flanking=search_range
    )
    
    if not extracted['sequence']:
        return None
    
    left_flank = extracted['left_flank']
    right_flank = extracted['right_flank']
    
    best_tsd = None
    best_score = 0
    
    # 搜索可能的TSD
    for tsd_len in range(min_length, min(max_length + 1, len(left_flank) + 1)):
        # 获取左侧潜在TSD
        left_tsd = left_flank[-tsd_len:] if len(left_flank) >= tsd_len else ""
        
        if not left_tsd:
            continue
        
        # 在右侧搜索匹配
        for i in range(min(tsd_len, len(right_flank))):
            right_tsd = right_flank[i:i+tsd_len]
            
            if len(right_tsd) != tsd_len:
                continue
            
            # 计算匹配度
            mismatches = sum(1 for j in range(tsd_len) if left_tsd[j] != right_tsd[j])
            
            if mismatches <= max_mismatches:
                # 计算得分（长度和匹配度的组合）
                score = tsd_len * (1 - mismatches / tsd_len)
                
                if score > best_score:
                    best_score = score
                    best_tsd = {
                        'sequence': left_tsd,
                        'length': tsd_len,
                        'mismatches': mismatches,
                        'type': 'perfect' if mismatches == 0 else 'imperfect',
                        'confidence': score / max_length,
                        'left_pos': len(left_flank) - tsd_len,
                        'right_pos': i,
                        'start_in_seq': -1,  # 需要根据具体情况调整
                        'end_in_seq': -1
                    }
    
    return best_tsd

def build_consensus_from_alignment(alignment, min_coverage: float = 0.6) -> str:
    """从多序列比对构建共识序列"""
    from collections import Counter
    
    if not alignment or len(alignment) == 0:
        return ""
    
    # 获取比对长度
    alignment_length = len(alignment[0])
    consensus = []
    
    for pos in range(alignment_length):
        # 获取该位置的所有碱基
        column = [seq[pos] for seq in alignment if pos < len(seq)]
        
        # 统计非gap字符
        non_gap_bases = [base for base in column if base != '-']
        
        # 计算覆盖度
        coverage = len(non_gap_bases) / len(column) if column else 0
        
        if coverage < min_coverage:
            # 覆盖度太低，跳过该位置
            continue
        
        if non_gap_bases:
            # 找到最常见的碱基
            base_counts = Counter(non_gap_bases)
            most_common_base = base_counts.most_common(1)[0][0]
            consensus.append(most_common_base)
    
    return ''.join(consensus)

def build_consensus_from_msa(msa_result, min_coverage: float = 0.3, use_iupac: bool = True,
                             quality_threshold: float = 0.5,
                             boundary_trim_coverage: float = 0.30) -> str:
    """Build consensus from MSA with coverage-based boundary trimming.

    After building the raw consensus, trims boundaries to the region where
    at least ``boundary_trim_coverage`` fraction of aligned copies have
    non-gap characters.  This prevents truncated copies from eroding
    consensus boundaries with gap-dominated columns.

    Args:
        msa_result: BioPython MultipleSeqAlignment or list of aligned strings.
        min_coverage: Minimum per-position coverage to include in consensus.
        use_iupac: Use IUPAC ambiguity codes for polymorphic sites.
        quality_threshold: Minimum dominant-base frequency for clean call.
        boundary_trim_coverage: Minimum coverage fraction to retain at
            consensus boundaries (default 0.30 = 30% of copies must have
            non-gap bases).

    Returns:
        Consensus sequence string.
    """
    from collections import Counter

    if not msa_result or len(msa_result) == 0:
        return ""

    # Extract sequences from BioPython alignment
    if hasattr(msa_result, '__iter__'):
        sequences = [str(record.seq) for record in msa_result]
    else:
        sequences = msa_result

    if not sequences:
        return ""

    n_seqs = len(sequences)
    alignment_length = len(sequences[0])

    # Step 1: Compute per-position coverage for boundary trimming
    coverage_profile = []
    for pos in range(alignment_length):
        column = [seq[pos].upper() for seq in sequences if pos < len(seq)]
        non_gap = sum(1 for base in column if base != '-' and base in 'ACGTN')
        cov_frac = non_gap / len(column) if column else 0
        coverage_profile.append(cov_frac)

    # Step 2: Find trim boundaries (region where coverage >= boundary_trim_coverage)
    trim_start = 0
    trim_end = alignment_length
    for i in range(alignment_length):
        if coverage_profile[i] >= boundary_trim_coverage:
            trim_start = i
            break
    for i in range(alignment_length - 1, -1, -1):
        if coverage_profile[i] >= boundary_trim_coverage:
            trim_end = i + 1
            break

    if trim_start >= trim_end:
        # Fallback: use full alignment
        trim_start = 0
        trim_end = alignment_length
    else:
        trimmed_cols = (trim_start) + (alignment_length - trim_end)
        if trimmed_cols > 0:
            logger.debug(f"Boundary trimming: removed {trim_start} left + "
                         f"{alignment_length - trim_end} right columns "
                         f"(coverage < {boundary_trim_coverage:.0%})")

    # Step 3: Log copy type diagnostics
    full_length_count = 0
    five_prime_trunc = 0
    three_prime_trunc = 0
    for seq in sequences:
        # A copy is "full length" if it has non-gap bases at both boundaries
        has_left = any(seq[i] != '-' and seq[i].upper() in 'ACGTN'
                       for i in range(min(20, len(seq))))
        has_right = any(seq[i] != '-' and seq[i].upper() in 'ACGTN'
                        for i in range(max(0, len(seq) - 20), len(seq)))
        if has_left and has_right:
            full_length_count += 1
        elif has_left and not has_right:
            five_prime_trunc += 1
        elif has_right and not has_left:
            three_prime_trunc += 1
        # else: both-end truncated (very short fragment)

    if n_seqs > 1:
        logger.debug(f"Copy types in MSA ({n_seqs} total): "
                     f"{full_length_count} full-length, "
                     f"{five_prime_trunc} 5'-truncated, "
                     f"{three_prime_trunc} 3'-truncated")

    # IUPAC ambiguous nucleotide codes
    iupac_codes = {
        frozenset(['A', 'G']): 'R',
        frozenset(['C', 'T']): 'Y',
        frozenset(['G', 'C']): 'S',
        frozenset(['A', 'T']): 'W',
        frozenset(['G', 'T']): 'K',
        frozenset(['A', 'C']): 'M',
        frozenset(['C', 'G', 'T']): 'B',
        frozenset(['A', 'G', 'T']): 'D',
        frozenset(['A', 'C', 'T']): 'H',
        frozenset(['A', 'C', 'G']): 'V',
        frozenset(['A', 'C', 'G', 'T']): 'N'
    }

    # Step 4: Build consensus within trimmed boundaries
    consensus = []
    for pos in range(trim_start, trim_end):
        column = [seq[pos].upper() for seq in sequences if pos < len(seq)]
        non_gap_bases = [base for base in column if base != '-' and base in 'ACGTN']
        coverage = len(non_gap_bases) / len(column) if column else 0

        if coverage < min_coverage:
            continue
        if not non_gap_bases:
            continue

        base_counts = Counter(non_gap_bases)
        total_bases = len(non_gap_bases)
        most_common = base_counts.most_common()
        dominant_base = most_common[0][0]
        dominant_freq = most_common[0][1] / total_bases

        if dominant_freq >= quality_threshold:
            consensus.append(dominant_base)
        elif use_iupac and len(most_common) > 1:
            significant_bases = set()
            cumulative_freq = 0
            for base, count in most_common:
                freq = count / total_bases
                if freq >= 0.2:
                    significant_bases.add(base)
                    cumulative_freq += freq
                if cumulative_freq >= 0.8:
                    break
            if len(significant_bases) > 1:
                iupac_key = frozenset(significant_bases)
                iupac_code = iupac_codes.get(iupac_key, 'N')
                consensus.append(iupac_code)
            else:
                consensus.append(dominant_base)
        else:
            consensus.append(dominant_base)

    consensus_seq = ''.join(consensus)
    consensus_seq = consensus_seq.strip('N')

    return consensus_seq

def extend_consensus_boundaries(consensus_seq: str, extended_sequence: str, aggressive: bool = False) -> str:
    """Extend consensus boundaries using a reference genomic copy.

    Biology-aware boundary extension:
    - Requires 80% similarity (not 60%) to accept extensions, preventing
      host DNA flanking TE insertions from being incorporated.
    - No aggressive fallback that blindly uses longer sequences.
    - Maximum extension capped at 500bp per side.

    Args:
        consensus_seq: Current consensus sequence.
        extended_sequence: Longer reference sequence (genomic copy with flanks).
        aggressive: Deprecated, kept for API compatibility but ignored.

    Returns:
        Extended consensus or original if extension is unsupported.
    """
    if not consensus_seq or not extended_sequence:
        return consensus_seq

    # Skip if extended_sequence is not meaningfully longer
    if len(extended_sequence) <= len(consensus_seq):
        return consensus_seq

    # IUPAC-aware base matching
    def matches_iupac(base1, base2):
        base1, base2 = base1.upper(), base2.upper()
        if base1 == base2 or base1 == 'N' or base2 == 'N':
            return True
        iupac_matches = {
            'R': 'AG', 'Y': 'CT', 'S': 'GC', 'W': 'AT',
            'K': 'GT', 'M': 'AC', 'B': 'CGT', 'D': 'AGT',
            'H': 'ACT', 'V': 'ACG', 'N': 'ACGT'
        }
        if base1 in iupac_matches:
            return base2 in iupac_matches[base1]
        if base2 in iupac_matches:
            return base1 in iupac_matches[base2]
        return False

    # Find best alignment position using sliding window
    search_range = len(extended_sequence) - len(consensus_seq) + 1
    if search_range <= 0:
        return consensus_seq

    best_match_pos = -1
    best_score = 0.0
    step = 1 if search_range < 1000 else max(1, search_range // 500)

    for i in range(0, search_range, step):
        # Quick screen on first 100 bases
        quick_checks = min(100, len(consensus_seq))
        quick_matches = sum(
            1 for j in range(quick_checks)
            if i + j < len(extended_sequence) and matches_iupac(consensus_seq[j], extended_sequence[i + j])
        )
        if quick_checks > 0 and quick_matches / quick_checks < 0.6:
            continue

        # Full comparison
        matches = 0
        comparisons = 0
        for j in range(len(consensus_seq)):
            if i + j < len(extended_sequence):
                comparisons += 1
                if matches_iupac(consensus_seq[j], extended_sequence[i + j]):
                    matches += 1

        if comparisons > 0:
            score = matches / comparisons
            if score > best_score:
                best_score = score
                best_match_pos = i

    # Require 80% similarity to accept extension (prevents host DNA contamination)
    if best_match_pos < 0 or best_score < 0.80:
        logger.debug(f"Boundary extension rejected: best_score={best_score:.2f} < 0.80 threshold")
        return consensus_seq

    # Cap extension at 500bp per side
    max_extension = 500

    # Left extension
    left_avail = best_match_pos
    left_ext_len = min(left_avail, max_extension)
    left_extension = extended_sequence[best_match_pos - left_ext_len:best_match_pos]

    # Right extension
    right_start = best_match_pos + len(consensus_seq)
    right_avail = len(extended_sequence) - right_start
    right_ext_len = min(right_avail, max_extension)
    right_extension = extended_sequence[right_start:right_start + right_ext_len]

    extended_consensus = left_extension + consensus_seq + right_extension
    extended_consensus = extended_consensus.strip('N')

    if len(extended_consensus) > len(consensus_seq):
        logger.debug(f"Extended consensus from {len(consensus_seq)}bp to {len(extended_consensus)}bp "
                     f"(match={best_score:.2f}, left=+{left_ext_len}bp, right=+{right_ext_len}bp)")

    return extended_consensus

def calculate_quick_similarity(seq1: str, seq2: str) -> float:
    """快速计算两个序列的相似度（采样方法）"""
    if not seq1 or not seq2:
        return 0.0
    
    # 采样比较，每100bp取一个点
    sample_size = min(100, min(len(seq1), len(seq2)))
    step = max(1, min(len(seq1), len(seq2)) // sample_size)
    
    matches = 0
    comparisons = 0
    
    for i in range(0, min(len(seq1), len(seq2)), step):
        if i < len(seq1) and i < len(seq2):
            comparisons += 1
            if seq1[i].upper() == seq2[i].upper() or seq1[i] == 'N' or seq2[i] == 'N':
                matches += 1
    
    return matches / comparisons if comparisons > 0 else 0.0



def run_cdhit(sequences: List[Dict], identity_threshold: float, config) -> List[Dict]:
    """使用CD-HIT-EST进行序列去冗余"""
    # 创建临时输入文件
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as input_file:
        for i, seq_data in enumerate(sequences):
            input_file.write(f">{seq_data['id']}\n{seq_data['sequence']}\n")
        input_file_path = input_file.name
    
    # 创建临时输出文件
    output_file_path = input_file_path + '.cdhit'
    
    try:
        # 运行CD-HIT-EST
        cmd = [
            config.cdhit_exe,
            '-i', input_file_path,
            '-o', output_file_path,
            '-c', str(identity_threshold),
            '-n', '10' if identity_threshold >= 0.95 else '8',
            '-T', str(config.threads),
            '-M', '0',  # 不限制内存
            '-d', '0'   # 完整描述
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        # 读取去冗余后的序列
        non_redundant = []
        with open(output_file_path, 'r') as f:
            for record in SeqIO.parse(f, 'fasta'):
                # 找到原始序列数据
                for seq_data in sequences:
                    if seq_data['id'] == record.id:
                        non_redundant.append(seq_data)
                        break
        
        return non_redundant
        
    except subprocess.CalledProcessError as e:
        logger.error(f"CD-HIT-EST failed: {e}")
        return sequences  # 返回原始序列
    finally:
        # 清理临时文件
        for file_path in [input_file_path, output_file_path, output_file_path + '.clstr']:
            if os.path.exists(file_path):
                os.unlink(file_path)