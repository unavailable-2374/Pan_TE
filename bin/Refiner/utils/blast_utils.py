import logging
import subprocess
import tempfile
import os
from typing import List, Dict, Set
from collections import defaultdict

logger = logging.getLogger(__name__)

def run_all_vs_all_blast(sequences: List[Dict], config) -> List[Dict]:
    """运行all-vs-all BLAST比较"""
    # 创建临时FASTA文件
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as fasta_file:
        for seq_data in sequences:
            fasta_file.write(f">{seq_data['id']}\n{seq_data['sequence']}\n")
        fasta_file_path = fasta_file.name
    
    # 创建BLAST数据库
    db_path = fasta_file_path + '.db'
    try:
        # 构建数据库
        makedb_cmd = [
            config.makeblastdb_exe,
            '-in', fasta_file_path,
            '-dbtype', 'nucl',
            '-out', db_path
        ]
        subprocess.run(makedb_cmd, capture_output=True, text=True, check=True)
        
        # 运行BLAST
        blast_cmd = [
            config.blastn_exe,
            '-query', fasta_file_path,
            '-db', db_path,
            '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen',
            '-num_threads', str(config.threads),
            '-max_target_seqs', '1000'
        ]
        
        result = subprocess.run(blast_cmd, capture_output=True, text=True, check=True)
        
        # 解析结果
        blast_results = []
        for line in result.stdout.strip().split('\n'):
            if not line:
                continue
            
            fields = line.split('\t')
            if len(fields) >= 14:
                query_id = fields[0]
                subject_id = fields[1]
                
                # 跳过自我比对
                if query_id == subject_id:
                    continue
                
                identity = float(fields[2])
                align_length = int(fields[3])
                query_length = int(fields[12])
                subject_length = int(fields[13])
                
                # 计算覆盖度
                query_coverage = align_length / query_length if query_length > 0 else 0
                subject_coverage = align_length / subject_length if subject_length > 0 else 0
                coverage = min(query_coverage, subject_coverage)
                
                blast_results.append({
                    'query': query_id,
                    'subject': subject_id,
                    'identity': identity,
                    'coverage': coverage,
                    'align_length': align_length,
                    'evalue': float(fields[10]),
                    'score': float(fields[11])
                })
        
        return blast_results
        
    except subprocess.CalledProcessError as e:
        logger.error(f"BLAST all-vs-all failed: {e}")
        return []
    finally:
        # 清理临时文件
        for file_path in [fasta_file_path, db_path + '.nhr', db_path + '.nin', db_path + '.nsq']:
            if os.path.exists(file_path):
                os.unlink(file_path)

def extract_connected_components(similarity_graph: Dict[str, List[str]]) -> List[Set[str]]:
    """从相似性图中提取连通分量"""
    visited = set()
    components = []
    
    def dfs(node: str, component: Set[str]):
        """深度优先搜索"""
        if node in visited:
            return
        visited.add(node)
        component.add(node)
        
        for neighbor in similarity_graph.get(node, []):
            if neighbor not in visited:
                dfs(neighbor, component)
    
    # 找到所有连通分量
    for node in similarity_graph:
        if node not in visited:
            component = set()
            dfs(node, component)
            if component:
                components.append(component)
    
    return components

def run_all_vs_all_blast_chunked(query_sequences: List[Dict], 
                                 db_sequences: List[Dict], 
                                 config) -> List[Dict]:
    """分块运行BLAST比较以减少内存使用"""
    # 创建数据库文件
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as db_file:
        for seq_data in db_sequences:
            db_file.write(f">{seq_data['id']}\n{seq_data['sequence']}\n")
        db_file_path = db_file.name
    
    # 创建查询文件
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as query_file:
        for seq_data in query_sequences:
            query_file.write(f">{seq_data['id']}\n{seq_data['sequence']}\n")
        query_file_path = query_file.name
    
    # 创建BLAST数据库
    db_path = db_file_path + '.db'
    try:
        # 构建数据库
        makedb_cmd = [
            config.makeblastdb_exe,
            '-in', db_file_path,
            '-dbtype', 'nucl',
            '-out', db_path
        ]
        subprocess.run(makedb_cmd, capture_output=True, text=True, check=True)
        
        # 运行BLAST - 使用更少的线程以减少内存
        blast_cmd = [
            config.blastn_exe,
            '-query', query_file_path,
            '-db', db_path,
            '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen',
            '-num_threads', str(min(4, config.threads)),  # 限制线程数
            '-max_target_seqs', '500',  # 减少目标序列数
            '-evalue', '1e-5'  # 更严格的E值阈值
        ]
        
        result = subprocess.run(blast_cmd, capture_output=True, text=True, check=True)
        
        # 解析结果
        blast_results = []
        for line in result.stdout.strip().split('\n'):
            if not line:
                continue
            
            fields = line.split('\t')
            if len(fields) >= 14:
                query_id = fields[0]
                subject_id = fields[1]
                
                # 跳过自我比对
                if query_id == subject_id:
                    continue
                
                identity = float(fields[2])
                align_length = int(fields[3])
                query_length = int(fields[12])
                subject_length = int(fields[13])
                
                # 计算覆盖度
                query_coverage = align_length / query_length if query_length > 0 else 0
                subject_coverage = align_length / subject_length if subject_length > 0 else 0
                coverage = min(query_coverage, subject_coverage)
                
                blast_results.append({
                    'query': query_id,
                    'subject': subject_id,
                    'identity': identity,
                    'coverage': coverage,
                    'align_length': align_length,
                    'evalue': float(fields[10]),
                    'score': float(fields[11])
                })
        
        return blast_results
        
    except subprocess.CalledProcessError as e:
        logger.error(f"BLAST chunked failed: {e}")
        return []
    finally:
        # 清理临时文件
        for file_path in [query_file_path, db_file_path, 
                         db_path + '.nhr', db_path + '.nin', db_path + '.nsq']:
            if os.path.exists(file_path):
                os.unlink(file_path)


def run_blast_segments(segment: str, genome_file: str, config) -> List[Dict]:
    """对序列片段进行BLAST搜索"""
    # 创建临时查询文件
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as query_file:
        query_file.write(f">segment\n{segment}\n")
        query_file_path = query_file.name
    
    try:
        # 运行BLAST
        cmd = [
            config.blastn_exe,
            '-query', query_file_path,
            '-subject', genome_file,
            '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore',
            '-max_target_seqs', '10',
            '-evalue', '1e-10'
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        # 解析结果
        hits = []
        for line in result.stdout.strip().split('\n'):
            if not line:
                continue
            
            fields = line.split('\t')
            if len(fields) >= 12:
                hits.append({
                    'chrom': fields[1],
                    'identity': float(fields[2]),
                    'length': int(fields[3]),
                    'start': int(fields[8]),
                    'end': int(fields[9]),
                    'evalue': float(fields[10]),
                    'score': float(fields[11])
                })
        
        return hits
        
    except subprocess.CalledProcessError as e:
        logger.error(f"BLAST segment search failed: {e}")
        return []
    finally:
        # 清理临时文件
        if os.path.exists(query_file_path):
            os.unlink(query_file_path)