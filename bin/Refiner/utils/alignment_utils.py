import logging
import subprocess
import tempfile
import os
from typing import List, Dict, Any
from Bio import SeqIO, AlignIO
from Bio.Align import MultipleSeqAlignment

logger = logging.getLogger(__name__)

def run_repeatmasker_batch_detailed(sequences: List[Dict], genome_file: str, 
                                   params: Dict, config) -> Dict[str, Dict]:
    """批量运行RepeatMasker并返回详细结果（包括所有hits）"""
    results = {}
    
    # 创建临时库文件
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as lib_file:
        for seq_data in sequences:
            lib_file.write(f">{seq_data['id']}\n{seq_data['sequence']}\n")
        lib_file_path = lib_file.name
    
    # 创建临时输出目录
    temp_dir = tempfile.mkdtemp(prefix='rm_batch_')
    
    try:
        # 构建RepeatMasker命令
        cmd = [config.repeatmasker_exe]
        
        # 添加快速模式参数（如果启用）
        if getattr(config, 'repeatmasker_quick', False):
            cmd.append('-q')
        
        # 添加参数
        if params.get('s'):
            cmd.append('-s')
        if params.get('no_is'):
            cmd.append('-no_is')
        if params.get('nolow'):
            cmd.append('-nolow')
        if 'cutoff' in params:
            cmd.extend(['-cutoff', str(params['cutoff'])])
        if 'pa' in params:
            cmd.extend(['-pa', str(params['pa'])])
        
        cmd.extend([
            '-lib', lib_file_path,
            '-dir', temp_dir,
            genome_file
        ])
        
        # 运行RepeatMasker
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            logger.error(f"RepeatMasker failed: {result.stderr}")
            return results
        
        # 解析输出文件获取详细hits
        out_file = os.path.join(temp_dir, os.path.basename(genome_file) + '.out')
        if os.path.exists(out_file):
            with open(out_file, 'r') as f:
                lines = f.readlines()
                
            # 初始化每个序列的结果
            for seq_data in sequences:
                results[seq_data['id']] = {
                    'copy_number': 0,
                    'hits': [],  # 详细的hit列表
                    'total_identity': 0,
                    'avg_identity': 0
                }
            
            # 解析每个hit
            for line in lines[3:]:  # 跳过header
                if line.strip():
                    parts = line.split()
                    if len(parts) >= 15:
                        repeat_name = parts[9]
                        
                        # 找到对应的序列
                        if repeat_name in results:
                            hit = {
                                'score': int(parts[0]) if parts[0].isdigit() else 0,
                                'div': float(parts[1]) if parts[1].replace('.','').isdigit() else 0,
                                'del': float(parts[2]) if parts[2].replace('.','').isdigit() else 0,
                                'ins': float(parts[3]) if parts[3].replace('.','').isdigit() else 0,
                                'chrom': parts[4],
                                'start': int(parts[5]) if parts[5].isdigit() else 0,
                                'end': int(parts[6]) if parts[6].isdigit() else 0,
                                'strand': parts[8],
                                'identity': 100 - float(parts[1]) if parts[1].replace('.','').isdigit() else 0
                            }
                            
                            results[repeat_name]['hits'].append(hit)
                            results[repeat_name]['copy_number'] += 1
                            results[repeat_name]['total_identity'] += hit['identity']
            
            # 计算平均相似度
            for seq_id in results:
                if results[seq_id]['copy_number'] > 0:
                    results[seq_id]['avg_identity'] = (
                        results[seq_id]['total_identity'] / results[seq_id]['copy_number']
                    )
        
        return results
        
    except Exception as e:
        logger.error(f"RepeatMasker batch processing failed: {e}")
        return results
    finally:
        # 清理临时文件
        if os.path.exists(lib_file_path):
            os.unlink(lib_file_path)
        # 清理临时目录
        import shutil
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)

def run_repeatmasker_batch(sequences: List[Dict], genome_file: str, 
                          params: Dict, config) -> Dict[str, Dict]:
    """批量运行RepeatMasker"""
    results = {}
    
    # 创建临时库文件
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as lib_file:
        for seq_data in sequences:
            lib_file.write(f">{seq_data['id']}\n{seq_data['sequence']}\n")
        lib_file_path = lib_file.name
    
    # 创建临时输出目录
    temp_dir = tempfile.mkdtemp(prefix='rm_batch_')
    
    try:
        # 构建RepeatMasker命令
        cmd = [config.repeatmasker_exe]
        
        # 添加快速模式参数（如果启用）
        if getattr(config, 'repeatmasker_quick', False):
            cmd.append('-q')
        
        # 添加参数
        if params.get('s'):
            cmd.append('-s')
        if params.get('no_is'):
            cmd.append('-no_is')
        if params.get('nolow'):
            cmd.append('-nolow')
        if 'cutoff' in params:
            cmd.extend(['-cutoff', str(params['cutoff'])])
        if 'pa' in params:
            cmd.extend(['-pa', str(params['pa'])])
        
        cmd.extend([
            '-lib', lib_file_path,
            '-dir', temp_dir,
            genome_file
        ])
        
        # 运行RepeatMasker
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        # 解析输出文件
        out_file = os.path.join(temp_dir, os.path.basename(genome_file) + '.out')
        if os.path.exists(out_file):
            results = parse_repeatmasker_output(out_file, sequences)
        
    except Exception as e:
        logger.error(f"RepeatMasker batch run failed: {e}")
    finally:
        # 清理临时文件
        if os.path.exists(lib_file_path):
            os.unlink(lib_file_path)
        # 清理临时目录
        import shutil
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)
    
    return results

def run_repeatmasker_single(sequence: Dict, genome_file: str, 
                           params: Dict, config) -> List[Dict]:
    """对单个序列运行RepeatMasker"""
    # 创建临时库文件
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as lib_file:
        lib_file.write(f">{sequence['id']}\n{sequence['sequence']}\n")
        lib_file_path = lib_file.name
    
    # 创建临时输出目录
    temp_dir = tempfile.mkdtemp(prefix='rm_single_')
    
    hits = []
    
    try:
        # 构建RepeatMasker命令
        cmd = [config.repeatmasker_exe]
        
        # 添加快速模式参数（如果启用）
        if getattr(config, 'repeatmasker_quick', False):
            cmd.append('-q')
        
        # 添加参数
        if params.get('s'):
            cmd.append('-s')
        if params.get('no_is'):
            cmd.append('-no_is')
        if params.get('nolow'):
            cmd.append('-nolow')
        if 'cutoff' in params:
            cmd.extend(['-cutoff', str(params['cutoff'])])
        
        cmd.extend([
            '-lib', lib_file_path,
            '-dir', temp_dir,
            genome_file
        ])
        
        # 运行RepeatMasker
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        # 解析输出文件
        out_file = os.path.join(temp_dir, os.path.basename(genome_file) + '.out')
        if os.path.exists(out_file):
            hits = parse_repeatmasker_hits(out_file)
        
    except Exception as e:
        logger.error(f"RepeatMasker single run failed: {e}")
    finally:
        # 清理临时文件
        if os.path.exists(lib_file_path):
            os.unlink(lib_file_path)
        # 清理临时目录
        import shutil
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)
    
    return hits

def parse_repeatmasker_output(out_file: str, sequences: List[Dict]) -> Dict[str, Dict]:
    """解析RepeatMasker输出文件"""
    results = {}
    
    # 初始化结果
    for seq_data in sequences:
        results[seq_data['id']] = {
            'copy_number': 0,
            'total_length': 0,
            'avg_identity': 0,
            'hits': []
        }
    
    try:
        with open(out_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('SW') or line.startswith('score'):
                    continue
                
                fields = line.split()
                if len(fields) >= 15:
                    # 解析字段
                    repeat_name = fields[9]
                    
                    # 找到对应的序列
                    for seq_data in sequences:
                        if seq_data['id'] == repeat_name:
                            identity = 100 - float(fields[1])  # 转换为相似度
                            length = int(fields[6]) - int(fields[5]) + 1
                            
                            # 处理染色体名称中可能包含的坐标信息
                            chrom = fields[4]
                            if ':' in chrom:
                                chrom = chrom.split(':')[0]
                            
                            results[repeat_name]['copy_number'] += 1
                            results[repeat_name]['total_length'] += length
                            results[repeat_name]['hits'].append({
                                'chrom': chrom,
                                'start': int(fields[5]),
                                'end': int(fields[6]),
                                'identity': identity,
                                'length': length
                            })
                            break
        
        # 计算平均相似度
        for seq_id, data in results.items():
            if data['hits']:
                identities = [hit['identity'] for hit in data['hits']]
                data['avg_identity'] = sum(identities) / len(identities)
    
    except Exception as e:
        logger.error(f"Failed to parse RepeatMasker output: {e}")
    
    return results

def parse_repeatmasker_hits(out_file: str) -> List[Dict]:
    """解析RepeatMasker hits"""
    hits = []
    
    try:
        with open(out_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('SW') or line.startswith('score'):
                    continue
                
                fields = line.split()
                if len(fields) >= 15:
                    score = int(fields[0])
                    identity = 100 - float(fields[1])  # 转换为相似度
                    chrom = fields[4]
                    
                    # 处理染色体名称中可能包含的坐标信息
                    # 例如 "lcl|34_8:121301-121525" 应该只提取 "lcl|34_8"
                    if ':' in chrom:
                        chrom = chrom.split(':')[0]
                    
                    start = int(fields[5])
                    end = int(fields[6])
                    
                    hits.append({
                        'score': score,
                        'identity': identity,
                        'chrom': chrom,
                        'start': start,
                        'end': end
                    })
    
    except Exception as e:
        logger.error(f"Failed to parse RepeatMasker hits: {e}")
    
    return hits

def run_mafft(sequences: List[Dict], algorithm: str = 'localpair',
             maxiterate: int = 1000, adjustdirection: bool = True,
             thread: int = None, config = None) -> MultipleSeqAlignment:
    """运行MAFFT多序列比对"""
    # 创建临时输入文件
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as input_file:
        for seq_data in sequences:
            input_file.write(f">{seq_data['id']}\n{seq_data['sequence']}\n")
        input_file_path = input_file.name
    
    # 创建临时输出文件
    output_file_path = input_file_path + '.aln'
    
    try:
        # 构建MAFFT命令
        cmd = [config.mafft_exe if config else 'mafft']
        
        # 算法选择和TE特异参数
        if algorithm == 'localpair':
            cmd.extend(['--localpair'])
            # L-INS-i: 对于高相似度TE，使用更精确的参数
            cmd.extend(['--maxiterate', str(min(maxiterate, 2000))])
        elif algorithm == 'globalpair':
            cmd.extend(['--globalpair'])
        elif algorithm == 'genafpair':
            cmd.extend(['--genafpair'])
            # E-INS-i: 对于有大插入/缺失的TE，允许更多迭代
            cmd.extend(['--maxiterate', str(min(maxiterate, 5000))])
        elif algorithm == 'auto':
            cmd.extend(['--auto'])
            # 自动模式：快速但质量较低
            cmd.extend(['--maxiterate', '100'])
        else:
            # 默认参数
            if maxiterate:
                cmd.extend(['--maxiterate', str(maxiterate)])
        
        # TE序列的方向性处理：始终启用
        if adjustdirection:
            cmd.append('--adjustdirection')
        
        # TE特异的额外参数
        if len(sequences) > 50:
            # 对于大数据集，使用更宽松的gap penalty
            cmd.extend(['--op', '2.0'])  # 增加gap opening penalty
            cmd.extend(['--ep', '0.5'])  # 增加gap extension penalty
        
        # 使用传入的thread参数或配置中的线程数（限制MAFFT线程数避免线程爆炸）
        if thread is not None:
            threads = thread
        else:
            # 限制MAFFT线程数，避免线程爆炸
            # 当有多个工作进程时，每个MAFFT应该使用较少的线程
            max_mafft_threads = max(1, config.threads // 8) if config else 1
            threads = min(max_mafft_threads, 8)  # 每个MAFFT最多8个线程
        cmd.extend(['--thread', str(threads)])
        
        # 对于大数据集，使用更快的算法
        if len(sequences) > 100:
            if algorithm == 'localpair':
                cmd.append('--retree')
                cmd.append('2')
            # 对于超大数据集，强制使用最快的算法
            if len(sequences) > 200:
                cmd = [config.mafft_exe if config else 'mafft']
                cmd.extend(['--auto', '--thread', str(threads), '--quiet'])
                logger.warning(f"Using MAFFT auto mode for {len(sequences)} sequences")
        
        cmd.extend([input_file_path])
        
        # 运行MAFFT
        with open(output_file_path, 'w') as out_file:
            result = subprocess.run(cmd, stdout=out_file, stderr=subprocess.PIPE, text=True)
        
        if result.returncode != 0:
            logger.error(f"MAFFT failed: {result.stderr}")
            return None
        
        # 读取比对结果
        alignment = AlignIO.read(output_file_path, 'fasta')
        return alignment
        
    except Exception as e:
        logger.error(f"MAFFT execution failed: {e}")
        return None
    finally:
        # 清理临时文件
        for file_path in [input_file_path, output_file_path]:
            if os.path.exists(file_path):
                os.unlink(file_path)