import logging
from typing import Dict, List, Tuple, Optional
from Bio import SeqIO
from Bio.Seq import Seq
import subprocess
import tempfile
import os

logger = logging.getLogger(__name__)

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
    """从基因组中提取指定区域的序列"""
    logger.debug(f"Extracting sequence from {genome_file}: {chrom}:{start}-{end} ({strand})")
    
    # 检查输入参数
    if not os.path.exists(genome_file):
        logger.error(f"Genome file not found: {genome_file}")
        return ""
    
    if start < 0 or end < 0 or start >= end:
        logger.warning(f"Invalid coordinates: {chrom}:{start}-{end}")
        return ""
    
    # 首先收集所有可用的染色体名称
    available_chroms = []
    try:
        with open(genome_file, 'r') as f:
            for record in SeqIO.parse(f, 'fasta'):
                available_chroms.append(record.id)
                if record.id == chrom:
                    genome_seq = str(record.seq)
                    chrom_length = len(genome_seq)
                    logger.debug(f"Found chromosome {chrom} with length {chrom_length}")
                    
                    # 确保坐标在合理范围内
                    original_start, original_end = start, end
                    start = max(0, start)
                    end = min(len(genome_seq), end)
                    
                    if original_start != start or original_end != end:
                        # 只在调整幅度较大时发出警告（超过100bp）
                        adjustment = abs(original_start - start) + abs(original_end - end)
                        if adjustment > 100:
                            logger.warning(f"Large coordinate adjustment from {original_start}-{original_end} to {start}-{end} (adjusted by {adjustment}bp)")
                        else:
                            logger.debug(f"Minor coordinate adjustment from {original_start}-{original_end} to {start}-{end}")
                    
                    if start >= end:
                        logger.warning(f"Invalid adjusted coordinates: {start}>={end}")
                        return ""
                    
                    # 提取序列
                    extracted_seq = genome_seq[start:end]
                    logger.debug(f"Extracted sequence length: {len(extracted_seq)}")
                    
                    # 如果是负链，获取反向互补
                    if strand == '-':
                        extracted_seq = reverse_complement(extracted_seq)
                        logger.debug(f"Applied reverse complement for negative strand")
                    
                    # 检查序列质量
                    n_count = extracted_seq.count('N')
                    n_percent = n_count / len(extracted_seq) if len(extracted_seq) > 0 else 1.0
                    logger.debug(f"Extracted sequence: length={len(extracted_seq)}, N_content={n_percent:.3f}")
                    
                    return extracted_seq
    except Exception as e:
        logger.error(f"Error reading genome file {genome_file}: {e}")
        return ""
    
    # 如果没有找到染色体，提供更详细的错误信息
    logger.warning(f"Chromosome '{chrom}' not found in genome file")
    logger.debug(f"Available chromosomes: {available_chroms[:10]}{'...' if len(available_chroms) > 10 else ''}")
    
    # 尝试模糊匹配
    possible_matches = []
    for avail_chrom in available_chroms:
        if chrom in avail_chrom or avail_chrom in chrom:
            possible_matches.append(avail_chrom)
    
    if possible_matches:
        logger.info(f"Possible chromosome name matches: {possible_matches}")
    
    return ""

def extract_sequence_with_flanking(genome_file: str, chrom: str, 
                                  start: int, end: int, 
                                  flanking: int = 50) -> Dict:
    """从基因组提取序列及其侧翼序列"""
    sequences = {}
    with open(genome_file, 'r') as f:
        for record in SeqIO.parse(f, 'fasta'):
            if record.id == chrom:
                genome_seq = str(record.seq)
                
                # 计算提取范围
                extract_start = max(0, start - flanking)
                extract_end = min(len(genome_seq), end + flanking)
                
                # 提取序列
                left_flank = genome_seq[extract_start:start]
                target_seq = genome_seq[start:end]
                right_flank = genome_seq[end:extract_end]
                
                return {
                    'sequence': target_seq,
                    'left_flank': left_flank,
                    'right_flank': right_flank,
                    'full_sequence': left_flank + target_seq + right_flank
                }
    
    logger.warning(f"Chromosome {chrom} not found in genome file")
    return {'sequence': '', 'left_flank': '', 'right_flank': '', 'full_sequence': ''}

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

def build_consensus_from_msa(msa_result, min_coverage: float = 0.3, use_iupac: bool = True, quality_threshold: float = 0.5) -> str:
    """从多序列比对结果构建共识序列"""
    from collections import Counter
    
    if not msa_result or len(msa_result) == 0:
        return ""
    
    # 从BioPython MultipleSeqAlignment提取序列
    if hasattr(msa_result, '__iter__'):
        sequences = [str(record.seq) for record in msa_result]
    else:
        sequences = msa_result
    
    if not sequences:
        return ""
    
    alignment_length = len(sequences[0])
    consensus = []
    
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
    
    for pos in range(alignment_length):
        # 获取该位置的所有碱基
        column = [seq[pos].upper() for seq in sequences if pos < len(seq)]
        
        # 统计非gap字符
        non_gap_bases = [base for base in column if base != '-' and base in 'ACGTN']
        
        # 计算覆盖度
        coverage = len(non_gap_bases) / len(column) if column else 0
        
        if coverage < min_coverage:
            continue
        
        if not non_gap_bases:
            continue
            
        # 统计碱基频率
        base_counts = Counter(non_gap_bases)
        total_bases = len(non_gap_bases)
        
        # 找到最频繁的碱基
        most_common = base_counts.most_common()
        dominant_base = most_common[0][0]
        dominant_freq = most_common[0][1] / total_bases
        
        if dominant_freq >= quality_threshold:
            # 主导碱基频率足够高
            consensus.append(dominant_base)
        elif use_iupac and len(most_common) > 1:
            # 使用IUPAC模糊碱基
            significant_bases = set()
            cumulative_freq = 0
            
            for base, count in most_common:
                freq = count / total_bases
                if freq >= 0.2:  # 至少20%频率才考虑
                    significant_bases.add(base)
                    cumulative_freq += freq
                if cumulative_freq >= 0.8:  # 累计80%频率
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
    
    # 清理首尾的N
    consensus_seq = consensus_seq.strip('N')
    
    return consensus_seq

def extend_consensus_boundaries(consensus_seq: str, extended_sequence: str, aggressive: bool = True) -> str:
    """使用扩展序列来改进共识序列的边界 - 优化版本，更积极地扩展"""
    if not consensus_seq or not extended_sequence:
        return consensus_seq
    
    # 如果扩展序列本身就比共识序列长，直接返回扩展序列（积极策略）
    if aggressive and len(extended_sequence) > len(consensus_seq) * 1.2:
        # 检查是否有基本的相似性
        similarity = calculate_quick_similarity(consensus_seq, extended_sequence)
        if similarity > 0.6:  # 60%相似度
            logger.debug(f"Using full extended sequence ({len(extended_sequence)}bp vs {len(consensus_seq)}bp consensus)")
            return extended_sequence.strip('N')
    
    # 找到共识序列在扩展序列中的最佳匹配位置
    best_match_pos = -1
    best_score = 0
    
    # 简化的相似度匹配，支持IUPAC codes
    def matches_iupac(base1, base2):
        """检查两个碱基是否匹配（支持IUPAC codes）"""
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
    
    # 使用滑动窗口找到最佳匹配位置 - 扩大搜索范围
    search_range = len(extended_sequence) - len(consensus_seq) + 1
    
    # 使用采样搜索以提高效率
    step = 1 if search_range < 1000 else max(1, search_range // 500)
    
    for i in range(0, search_range, step):
        # 快速检查前100个碱基
        quick_matches = 0
        quick_checks = min(100, len(consensus_seq))
        
        for j in range(quick_checks):
            if i + j < len(extended_sequence):
                if matches_iupac(consensus_seq[j], extended_sequence[i + j]):
                    quick_matches += 1
        
        quick_score = quick_matches / quick_checks if quick_checks > 0 else 0
        
        # 如果快速检查分数太低，跳过详细检查
        if quick_score < 0.5:
            continue
        
        # 详细检查
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
    
    # 降低阈值以更积极地扩展 - 从70%降到60%
    if best_match_pos >= 0 and best_score > 0.6:
        # 获取左侧扩展 - 限制最大扩展长度避免过度扩展
        max_left_extension = min(best_match_pos, len(consensus_seq) // 2)  # 最多扩展50%长度
        left_extension = extended_sequence[max(0, best_match_pos - max_left_extension):best_match_pos]
        
        # 获取右侧扩展
        right_start = best_match_pos + len(consensus_seq)
        max_right_extension = min(len(extended_sequence) - right_start, len(consensus_seq) // 2)
        right_extension = extended_sequence[right_start:right_start + max_right_extension]
        
        # 组合扩展序列
        extended_consensus = left_extension + consensus_seq + right_extension
        
        # 清理首尾的N和低质量区域，但保留内部的N
        extended_consensus = extended_consensus.strip('N')
        
        if len(extended_consensus) > len(consensus_seq):
            logger.debug(f"Extended consensus from {len(consensus_seq)}bp to {len(extended_consensus)}bp ")
        
        return extended_consensus
    
    # 如果没有找到好的匹配，但扩展序列更长，考虑直接使用
    if aggressive and len(extended_sequence) > len(consensus_seq) * 1.1:
        logger.debug(f"No good match found, but using longer extended sequence anyway ({len(extended_sequence)}bp)")
        return extended_sequence.strip('N')
    
    # 返回原序列
    return consensus_seq

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