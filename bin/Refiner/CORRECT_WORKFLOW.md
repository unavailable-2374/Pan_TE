# 正确的TE共识序列构建流程

## 问题背景

原始实现存在逻辑错误：对RepeatScout序列本身进行聚类，而不是对它们在基因组中的拷贝进行聚类。

## 正确的流程

### Phase 1: 筛选与评分
- 加载RepeatScout输出序列
- 计算复杂度分数（DUST、Shannon熵等）
- 初步RepeatMasker评估覆盖度
- 分类为A/B/C三个等级

### Phase 2: 基于拷贝的共识构建（正确逻辑）

对每个高质量RepeatScout序列：

1. **找到基因组拷贝**
   - 使用RepeatMasker搜索基因组
   - 提取所有匹配的基因组片段
   - 记录位置、方向、相似度信息

2. **拷贝聚类识别亚家族**
   - 计算拷贝之间的相似度矩阵
   - 层次聚类识别亚家族
   - 动态阈值（基于距离分布）

3. **亚家族共识构建**
   - 对每个亚家族执行多序列比对(MAFFT)
   - 构建加权共识序列
   - 检测TSD和边界特征

### Phase 3: 质量控制（宽松版）
- 宽松的嵌合检测（只移除明显嵌合）
- 自适应质量筛选（多轨道评估）
- 两档冗余去除（95%和90%）

## 关键改进

1. **正确的聚类对象**：聚类基因组拷贝而不是RepeatScout序列
2. **亚家族识别**：一个RepeatScout序列可能产生多个亚家族
3. **更准确的共识**：基于实际基因组拷贝而不是种子序列

## 使用方法

### 使用正确逻辑（默认）
```bash
python main.py -r repeatscout.fa -g genome.fa -o output_dir -t 16
```

### 使用旧逻辑（对比测试）
```bash
python main.py -r repeatscout.fa -g genome.fa -o output_dir -t 16 --use-old-logic
```

### 清除缓存重新运行
```bash
rm -rf checkpoints/phase2_* cache/
python main.py -r repeatscout.fa -g genome.fa -o output_dir -t 16
```

## 预期结果

使用正确逻辑后：
- 更多的共识序列（包含亚家族）
- 更准确的边界定义
- 更好的拷贝数估计
- 减少假阳性嵌合检测

## 参数调优建议

对于不同类型的基因组：

### 高重复基因组
```python
config.min_copy_number = 5
config.identity_threshold = 0.80
```

### 低重复基因组
```python
config.min_copy_number = 2
config.identity_threshold = 0.70
```

### 古老/分化的TEs
```python
config.identity_threshold = 0.65
config.coverage_threshold = 0.40
```