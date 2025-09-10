"""
ID格式检查和诊断工具
用于检测RepeatMasker结果与基因组文件间的ID不匹配问题
"""

import logging
from typing import Dict, List, Set, Tuple
from collections import Counter
import re

logger = logging.getLogger(__name__)

def analyze_id_mismatch(hits: List[Dict], genome_file: str) -> Dict:
    """
    分析RepeatMasker hits与基因组文件间的ID不匹配问题
    
    Args:
        hits: RepeatMasker hits列表
        genome_file: 基因组文件路径
        
    Returns:
        分析结果字典
    """
    # 从hits中提取染色体ID
    hit_chroms = set()
    for hit in hits:
        chrom = hit.get('chrom', '')
        if chrom:
            # 处理可能包含坐标的ID
            if ':' in chrom:
                chrom = chrom.split(':')[0]
            hit_chroms.add(chrom)
    
    # 从基因组文件中提取染色体ID
    genome_chroms = set()
    try:
        from Bio import SeqIO
        for record in SeqIO.parse(genome_file, "fasta"):
            genome_chroms.add(record.id)
    except Exception as e:
        logger.error(f"Failed to read genome file {genome_file}: {e}")
        return {'error': str(e)}
    
    # 分析匹配情况
    exact_matches = hit_chroms & genome_chroms
    hit_only = hit_chroms - genome_chroms
    genome_only = genome_chroms - hit_chroms
    
    # 分析ID格式模式
    hit_patterns = _analyze_id_patterns(hit_chroms)
    genome_patterns = _analyze_id_patterns(genome_chroms)
    
    # 尝试模糊匹配
    fuzzy_matches = _find_fuzzy_matches(hit_only, genome_chroms)
    
    # 生成报告
    report = {
        'total_hits': len(hits),
        'unique_hit_chroms': len(hit_chroms),
        'unique_genome_chroms': len(genome_chroms),
        'exact_matches': len(exact_matches),
        'hit_only': len(hit_only),
        'genome_only': len(genome_only),
        'match_rate': len(exact_matches) / len(hit_chroms) if hit_chroms else 0,
        'hit_patterns': hit_patterns,
        'genome_patterns': genome_patterns,
        'fuzzy_matches': fuzzy_matches,
        'sample_hit_chroms': list(hit_chroms)[:10],
        'sample_genome_chroms': list(genome_chroms)[:10],
        'sample_unmatched_hits': list(hit_only)[:10],
        'sample_unmatched_genome': list(genome_only)[:10]
    }
    
    return report

def _analyze_id_patterns(ids: Set[str]) -> Dict:
    """分析ID的格式模式"""
    patterns = {
        'has_lcl_prefix': 0,
        'has_chr_prefix': 0,
        'has_ref_prefix': 0,
        'has_gi_prefix': 0,
        'has_underscore': 0,
        'has_colon': 0,
        'is_numeric': 0,
        'has_pipe': 0,
        'length_distribution': Counter()
    }
    
    for id_str in ids:
        if id_str.startswith('lcl|'):
            patterns['has_lcl_prefix'] += 1
        if id_str.startswith('chr'):
            patterns['has_chr_prefix'] += 1
        if id_str.startswith('ref|'):
            patterns['has_ref_prefix'] += 1
        if id_str.startswith('gi|'):
            patterns['has_gi_prefix'] += 1
        if '_' in id_str:
            patterns['has_underscore'] += 1
        if ':' in id_str:
            patterns['has_colon'] += 1
        if '|' in id_str:
            patterns['has_pipe'] += 1
        if id_str.isdigit():
            patterns['is_numeric'] += 1
        
        # 长度分布
        length_range = f"{len(id_str)//10*10}-{len(id_str)//10*10+9}"
        patterns['length_distribution'][length_range] += 1
    
    return patterns

def _find_fuzzy_matches(unmatched_hits: Set[str], genome_chroms: Set[str]) -> List[Tuple[str, str]]:
    """寻找模糊匹配的染色体ID对"""
    fuzzy_matches = []
    
    for hit_id in list(unmatched_hits)[:20]:  # 限制检查数量
        for genome_id in genome_chroms:
            if _ids_potentially_match(hit_id, genome_id):
                fuzzy_matches.append((hit_id, genome_id))
                if len(fuzzy_matches) >= 20:  # 限制结果数量
                    break
        if len(fuzzy_matches) >= 20:
            break
    
    return fuzzy_matches

def _ids_potentially_match(id1: str, id2: str) -> bool:
    """检查两个ID是否可能匹配"""
    # 提取数字部分
    nums1 = re.findall(r'\d+', id1)
    nums2 = re.findall(r'\d+', id2)
    
    # 如果数字部分相同，可能匹配
    if nums1 and nums2 and nums1 == nums2:
        return True
    
    # 检查移除前缀后是否匹配
    for prefix in ['lcl|', 'ref|', 'gi|', 'chr']:
        clean_id1 = id1[len(prefix):] if id1.startswith(prefix) else id1
        clean_id2 = id2[len(prefix):] if id2.startswith(prefix) else id2
        
        if clean_id1 == clean_id2:
            return True
    
    return False

def print_id_mismatch_report(report: Dict):
    """打印ID不匹配分析报告"""
    if 'error' in report:
        logger.error(f"Analysis failed: {report['error']}")
        return
    
    print("\n" + "="*60)
    print("Chromosome ID Mismatch Analysis Report")
    print("="*60)
    
    print(f"Total RepeatMasker hits: {report['total_hits']}")
    print(f"Unique chromosome IDs in hits: {report['unique_hit_chroms']}")
    print(f"Unique chromosome IDs in genome: {report['unique_genome_chroms']}")
    print(f"Exact matches: {report['exact_matches']}")
    print(f"Match rate: {report['match_rate']:.1%}")
    print(f"Unmatched hit IDs: {report['hit_only']}")
    print(f"Unused genome IDs: {report['genome_only']}")
    
    print("\nID Format Patterns in RepeatMasker hits:")
    for pattern, count in report['hit_patterns'].items():
        if isinstance(count, int) and count > 0:
            print(f"  {pattern}: {count}")
    
    print("\nID Format Patterns in Genome file:")
    for pattern, count in report['genome_patterns'].items():
        if isinstance(count, int) and count > 0:
            print(f"  {pattern}: {count}")
    
    if report['sample_unmatched_hits']:
        print("\nSample unmatched hit IDs:")
        for hit_id in report['sample_unmatched_hits']:
            print(f"  {hit_id}")
    
    if report['sample_unmatched_genome']:
        print("\nSample unmatched genome IDs:")
        for genome_id in report['sample_unmatched_genome']:
            print(f"  {genome_id}")
    
    if report['fuzzy_matches']:
        print("\nPotential fuzzy matches:")
        for hit_id, genome_id in report['fuzzy_matches'][:10]:
            print(f"  {hit_id} -> {genome_id}")
    
    print("\nRecommendations:")
    if report['match_rate'] < 0.1:
        print("⚠️  Very low match rate! Consider:")
        print("   1. Check if RepeatMasker used a different genome file")
        print("   2. Use ID mapping to convert between formats")
        print("   3. Verify genome file preprocessing steps")
    elif report['match_rate'] < 0.5:
        print("⚠️  Low match rate. ID mapping recommended.")
    else:
        print("✅ Match rate is acceptable.")
    
    print("="*60)

def quick_id_check(hits: List[Dict], genome_file: str):
    """快速ID检查并输出简要报告"""
    report = analyze_id_mismatch(hits, genome_file)
    
    if 'error' in report:
        logger.error(f"ID check failed: {report['error']}")
        return
    
    match_rate = report['match_rate']
    if match_rate < 0.1:
        logger.error(f"Severe ID mismatch: {match_rate:.1%} match rate")
        logger.error("RepeatMasker hits and genome file use incompatible ID formats")
        logger.info("Using ID mapping to resolve mismatches...")
    elif match_rate < 0.5:
        logger.warning(f"ID mismatch detected: {match_rate:.1%} match rate") 
        logger.info("Using ID mapping to improve compatibility...")
    else:
        logger.info(f"ID compatibility good: {match_rate:.1%} match rate")