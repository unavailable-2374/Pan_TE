# CRITICAL FIX - Phase 2 运行太快的问题

## 问题描述
Phase 2运行太快（仅15秒），没有进行扩展：
- 输入: 2174个序列
- 输出: 2174个序列（扩展比例1.00x）
- 全部是"Single"类型，没有"Multiple"

## 根本原因
Phase 2 `get_genome_copies`函数中的字段名错误：
- 使用了`hit['chr']`
- 但RepeatMasker结果中实际字段名是`hit['chrom']`
- 导致所有序列提取失败，没有找到任何拷贝

## 修复内容
文件：`phase2_consensus_expansion.py`

### 修改1（第283行）：
```python
# 错误：
hit['chr']
# 修正为：
hit['chrom']
```

### 修改2（第305行）：
```python
# 错误：
'chr': hit['chr'],
# 修正为：
'chr': hit['chrom'],
```

## 重新运行步骤

1. **上传修复后的文件**：
```bash
scp phase2_consensus_expansion.py user@server:/path/to/TCB/god/
```

2. **清除Phase 2的checkpoint**（必须！）：
```bash
rm -f checkpoints/phase2_complete.pkl
```

3. **重新运行Phase 2和后续步骤**：
```bash
# 如果Phase 1已完成，会自动从Phase 2开始
python main.py \
  -r /path/to/repeatscout_output.fa \
  -g /path/to/genome.fa \
  -o output_dir \
  -t 32 \
  --resume \
  --keep-checkpoints
```

## 预期结果

修复后Phase 2应该：
1. **运行时间更长**（几分钟到几十分钟，取决于数据量）
2. **产生更多共识序列**（扩展比例>1.5x）
3. **显示"Multiple"类型的统计**（不只是"Single"）

示例输出：
```
Phase 2 complete: Generated 3500 consensus sequences from 2174 inputs
Expansion ratio: 1.61x
Single: 1200, Multiple: 800, Failed: 174
```

## 验证修复

查看日志中应该出现：
```bash
grep "Found.*copies" pipeline.log
# 应该看到类似：
# DEBUG - seq_1: Found 15 copies
# DEBUG - seq_2: Found 8 copies
```

而不是：
```bash
# ERROR - Error extracting copy from hit: 'chr'
```

## Phase 3速度快的原因

Phase 3运行快（5秒）是正常的，因为：
1. 它主要使用CD-HIT进行序列聚类（已优化的C++程序）
2. 在relaxed模式下跳过了chimera检测
3. 使用并行处理（8个workers）
4. 2174个序列对CD-HIT来说是中等规模，可以快速处理

## 总结

这个修复至关重要，因为：
1. Phase 2是整个pipeline的核心，负责扩展和改进序列
2. 没有正确的Phase 2，最终TE库质量会很差
3. 字段名错误导致Phase 2完全失效，只是简单复制了输入序列