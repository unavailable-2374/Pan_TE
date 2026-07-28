---
title: "Pan_TE完整工作流、核心算法与论文方法学说明"
subtitle: "基于当前可执行实现、数据契约与拟南芥验证结果"
lang: zh-CN
toc-title: "目录"
---

# 文档定位与方法学边界

本文档描述 Pan_TE 当前可执行代码所采用的数据处理流程、算法、参数、质量控制和输出契约，可作为论文 Materials and Methods、软件方法说明和补充方法的事实来源。所有算法描述均以仓库中的 `bin/Pan_TE`、`bin/paf_processor.py`、`bin/gfa_processor.py`、`bin/rgfa_gaf.py`、`bin/gfa_annotator.py`、`bin/Refiner_mdl/` 及所调用的 `dtr` 命令为准。

Pan_TE 的核心任务是从一个基因组或多个组装中建立非冗余重复序列库。生产流程由三个互补发现轨道组成：Look4LTRs 识别具有典型长末端重复结构的候选；mdl-repeat 利用重复序列的压缩收益发现广谱重复家族，并由 Refiner_mdl 校正边界和共识；TE-looker 利用高频种子、成员聚类、片段拼接和家族门控补充远缘或破碎重复。三个轨道的结果经跨轨去冗余后形成最终重复库。分类是可选步骤；未提供分类模型时，产物是未分类的从头重复库，不自动等同于已经完成生物学分类的 TE 库。

TE-looker 的生产输出由 `dtr run` 及其门控契约决定。只有声明 `contract=te-looker-gated-v1`、`mode=gated` 且 `final_library=families.fasta` 的 provenance 才被直接接纳。没有该契约时，Pan_TE 对发现结果执行片段拼接、长度和成员数门控以及轨道内去冗余。Track 2 的窗口图、LSH、社区、POA、non-TE 守卫、copy catalog、结构审计和 G1–G6 审计输出属于部分证据链，`family_call=false`，不进入生产最终库。

# 1. 软件输入、运行模式与总体流程

## 1.1 输入模式

Pan_TE 接受四类互斥或组合输入。

1. 单基因组模式使用 `--genome genome.fa`，三个发现轨道在规范化后的基因组上运行。
2. PAF 泛基因组模式使用 `--paf all_vs_all.paf --paf-fasta assembly1.fa assembly2.fa ...`，从装配间共线比对的不对称缺口和 CIGAR 插入中提取候选序列。
3. GFA 泛基因组模式使用 `--gfa graph.gfa`。PGGB/smoothxg 图由内嵌 P/W 路径恢复样本等位序列；minigraph rGFA 由 `SR/SN/SO` 标签和 primary GAF 恢复样本 traversal。
4. 富集模式同时提供参考基因组和 PAF、GFA 或 VCF。候选插入被追加到规范化参考序列，三个发现轨道在该富集序列集上运行。

GFA 注释是独立于重复库发现的投影功能。`--anno-only --gfa graph.gfa --anno-manifest annotations.tsv` 将所有组装已有的 RepeatMasker、BED 或 GFF3 注释投影到图节点，不调用 RepeatMasker，也不依赖 `combine/TEs.fa`。

## 1.2 阶段依赖关系

工作流按以下有向依赖执行：

```text
输入合法性和环境预检
  -> 基因组规范化及 PAF/GFA/VCF 富集
  -> Look4LTRs 与 mdl-repeat/Refiner_mdl 并行运行
  -> TE-looker（可使用 mdl 共识构建的 HMM）
  -> 三轨合并和去冗余
  -> 可选分类
  -> 可选 GFA 全基因组注释投影
```

Look4LTRs 与 mdl-repeat 由两个任务并行运行。TE-looker 位于二者之后，因为其 HMM 输入可以由 Refiner_mdl 的 masking consensus 构建。阶段完成状态由 `genome.ok`、`ltr.ok`、`mdl_repeat.ok`、`te_looker.ok`、`combine.ok` 和 `anno.ok` 表示。GFA 发现模式额外保存输入路径、文件大小、mtime、图模式、GAF 路径和阈值指纹；指纹不一致时不复用旧的基因组阶段检查点。

## 1.3 坐标与标识符约定

内部 BED、PAF、GAF、occurrence 和投影 TSV 均采用 0-based、half-open 区间 $[s,e)$，长度为 $e-s$。GFF3 输入由 1-based inclusive 转为内部区间 $[start-1,end)$，输出 node-local GFF3 时再转换回 1-based inclusive。RepeatMasker `.out` 的 `query_begin/query_end` 同样转换为 $[begin-1,end)$。调用 `samtools faidx` 时使用 1-based inclusive 区间 `start+1-end`。

泛基因组样本名称优先采用 PanSN：`sample#haplotype#contig`。manifest 中的 `(sample,haplotype)` 是样本主键；`seqid_map` 用于显式连接注释 seqid 和图路径 contig。所有自动名称匹配都保留 unresolved 审计，歧义匹配不作猜测。

# 2. 环境预检、资源控制与输入规范化

## 2.1 环境预检

预检同时检查外部程序、Python/Perl 模块、分类模型、RepeatMasker/Dfam 资产和 TE-looker 内存。检查结果集中报告，使缺失依赖在长任务开始前暴露。`--anno-only` 只读取已有注释，因此不要求 RepeatMasker 可执行文件。

可用内存取以下上限中的最小值：命令行 `-M`、cgroup 限制、SLURM 分配和 `/proc/meminfo`。这种取最小值策略反映进程实际受到的最严格资源约束。

TE-looker seed harvest 对每个保留的 k-mer occurrence 使用 16-byte 打包条目。峰值内存估算为：

$$
M_{peak}=1.15\left[G+16G\rho(w)c\right]
$$

其中 $G$ 为基因组 bp 数，$w$ 为 minimizer window，$\rho(w)=1$（$w\leq1$）或 $2/(w+1)$（$w>1$），$c$ 为内存副本数。`dtr` 版本无法识别或版本低于 0.3.1 时取 $c=2$，版本不低于 0.3.1 时取 $c=1$。估算值超过可用上限时流程在 seed harvest 前终止，并报告所需内存和显式可选的 $w$；流程不自动改变 $w$，因为改变 $w$ 会改变种子集合。

## 2.2 FASTA 规范化

`clean_seq_fast` 以流式方式读取 FASTA，将序列转为大写，并把非 `A/C/G/T/N` 的退化碱基替换为 `N`。内部 contig 统一命名为 `chr1...chrN`，随后使用 `samtools faidx` 和项目索引程序建立随机访问索引。原始标识符与内部标识符的关系通过输入映射和 manifest 保存。

## 2.3 VCF 插入序列富集

VCF 模式只提取具有显式核苷酸序列的插入。非符号 ALT 直接提供插入序列；符号 `<INS>` 仅在 INFO 中含 `SVINSSEQ`、`INSSEQ`、`SEQ` 或 `CONSENSUS` 时可提取。没有序列的结构变异只表达事件存在，不能用于从头建立重复共识，因此不加入发现基因组。

# 3. PAF 泛基因组候选提取算法

## 3.1 输入验证和记录模型

PAF 记录表示为：

```text
PAFRecord(qname, qlen, qstart, qend, strand,
          tname, tlen, tstart, tend,
          matches, block_len, mapq, cigar)
```

解析器要求至少 12 个字段；长度、起止坐标、匹配数、block length 和 MAPQ 必须为合法整数；strand 必须为 `+` 或 `-`；query 和 target 名称必须存在于所给 FASTA 索引且记录长度一致。若提供 `--paf-mapping`，同一 assembly 内的旁系比对和不同染色体之间的比对被排除。记录按 `(qname,tname)` 分组，正负链分别建链。

默认保留条件为比对长度至少 5,000 bp、MAPQ 至少 20。比对器类型根据 `tp` 或 `cg` 标签识别为 minimap2，否则按不含 CIGAR 的 fastGA 类记录处理。

## 3.2 链构建

每个 query-target-strand 组按 `tstart` 升序排序。贪婪算法从每个未使用记录启动一条链，并扫描后续记录。对相邻记录 $a,b$，target gap 为：

$$
g_t=b.tstart-a.tend
$$

正链 query gap 为：

$$
g_q=b.qstart-a.qend
$$

负链随 target 前进时 query 坐标递减，因此：

$$
g_q=a.qstart-b.qend
$$

只有 $g_t$ 和 $g_q$ 均不小于 `-overlap_tolerance` 且不大于 `max_chain_gap` 时，$b$ 才加入该链。默认 overlap tolerance 为 200 bp，最大链间距为 500 kb。由于记录按 target 排序，一旦 $g_t>500$ kb 即终止该内层扫描。保留链至少含 3 条 alignment，且 target aligned length 总和至少为 $\max(50\text{ kb},2L_{max})$；默认候选最大长度 $L_{max}=30$ kb，因此默认总锚定长度门槛为 60 kb。

该贪婪算法最坏复杂度为 $O(n^2)$，排序后大 target gap 的提前终止限制了典型数据的有效比较数。链是候选提取的共线背景，不被解释为全局最优全基因组比对。

## 3.3 链间不对称缺口

对链内相邻 alignment 计算 $g_q$ 和 $g_t$，定义：

$$
g_{max}=\max(|g_q|,|g_t|),\quad
g_{min}=\min(|g_q|,|g_t|),\quad
r=g_{min}/g_{max}
$$

候选满足 $100\le g_{max}\le30,000$ bp 且 $r\le0.30$。若 $|g_q|>|g_t|$，插入位于 query；若 $|g_t|>|g_q|$，插入位于 target。负链 query 候选的端点按原 query forward 坐标重新排序。零长或反向区间被丢弃。

不对称条件的目的不是直接判断 TE，而是识别“一侧显著多出序列”的局部事件。候选仍需经过序列 QC 和三个重复发现轨道。

## 3.4 CIGAR 内部插入

minimap2 `cg:Z` 被完整重构，操作含义如下：`M/=/X` 同时推进 query 和 target；`I` 只推进 query 并产生 query-side insertion；`D` 只推进 target 并产生 target-side insertion；`N` 只推进 target；`S/H` 使该记录不进入 CIGAR 候选提取；`P` 不推进坐标。负链 query cursor 从 `qend` 向 `qstart` 递减。

解析后重新计算 query 和 target 消耗长度，并分别与 PAF 的 `qend-qstart`、`tend-tstart` 比较。CIGAR 字符串不完整、存在未知操作或坐标消耗不一致时，整条记录的 CIGAR 候选被拒绝，避免一个错误操作使后续插入坐标整体偏移。

## 3.5 候选合并、序列提取与 QC

来自链间缺口和 CIGAR 的候选统一转为 `InsertionCandidate(seqname,start,end,strand,source,partner)`。同一序列上重叠或近邻的候选按坐标合并，并保留来源证据。所有区间再次验证边界不越过 FASTA 长度。

序列通过 `samtools faidx` 批量提取；负链候选取反向互补。默认序列 QC 包括 100–30,000 bp 长度、N 比例、单核苷酸熵、二核苷酸熵、低复杂度窗口比例和 DUST。通过 QC 的候选使用 CD-HIT-EST 去冗余，并以独立 FASTA record 追加到分析基因组；`--gfa-spacer-length` 仅在下游需要单 scaffold 时引入人工 N 间隔。

# 4. GFA 到发现序列的转换

## 4.1 图格式判别

图解析器检查 S 记录的 `SR:i` 标签。具有 rGFA rank 语义的图按 minigraph rGFA 处理；具有内嵌 P/W sample paths 的图按 PGGB/smoothxg 处理。S 记录必须提供序列或正的 `LN`；路径引用不存在的节点、无法解释的 overlap 和不一致坐标均构成硬错误或 unresolved 记录。

## 4.2 PGGB/smoothxg 参考路径和等位区间

P 路径的节点步表示为 `node+` 或 `node-`，W walk 使用 `>node` 或 `<node`。节点序列按方向拼接。相邻节点 overlap 来自 P overlap 字段或 bidirected L link；只有完全由 `M/=/X` 组成的 overlap CIGAR可转换为确定的 overlap 长度。路径坐标前缀为：

$$
p_0=0,\qquad p_{i+1}=p_i+|S_i|-o_i
$$

其中 $o_i$ 是节点 $i$ 与 $i+1$ 的 overlap。bidirected 边同时检查正向键 `(A,o_A,B,o_B)` 与反向互补键 `(B,flip(o_B),A,flip(o_A))`。

参考样本由 `--gfa-ref` 或路径支持自动确定。参考路径的节点顺序构成 anchor rank。对每条非参考路径，提取方向一致、rank 单调递增的参考节点作为锚点，并取最大共线锚链。相邻锚点之间的非参考 node walk 构成 alternate allele；同一区间的 reference internal nodes 构成 reference allele。两者序列不同且 alternate 长度在 100–30,000 bp 范围时生成候选。

长于 30 kb 的 allele 按窗口拆分，窗口默认有 1,000 bp overlap；末端窗口向前移动以避免产生小于 100 bp 的尾片段。候选记录 reference path、左右锚点、样本路径、方向、node IDs、reference/alternate span、序列 SHA-256 和 QC 状态。

等位序列使用磁盘 SQLite 按 locus 和 sequence SHA-256 聚合。相同位点和相同序列只保留一个代表，同时累计 occurrence 数、支持样本数和样本集合。该设计避免完整图的所有 candidate 常驻内存。

## 4.3 minigraph rGFA 与 GAF 等位提取

严格 rGFA 后端要求每个 S 记录含内嵌序列以及 `SR/SN/SO`，L overlap 为 `0M`。解析器建立五类索引：node ID 到 segment；双向 link；`(SN,start,end)` 到 node；每个 source path 的有序区间；SR=0 reference contig 的 node rank 和 prefix length。SR=0 节点必须在各 contig 上无重叠且坐标连续。

GAF 只接受 primary 记录，并验证 mandatory fields、query/path 范围、`cg` 消耗、resolved path length 和每条相邻图边。path 支持三种表示：直接 node walk、`source:start-end` 坐标 token walk、命名 source path。坐标 token 只有在精确命中一个节点或能展开为完整连续节点序列时才被接受。

对每条高质量记录，只使用被 `path_start/path_end` 完整覆盖的节点作等位序列锚定，部分终端节点不用于构造完整 allele。负链 traversal 反转节点顺序并翻转节点方向，标准化到 query-forward。

SR=0 节点按 reference rank 构建最大递增锚链。每对相邻参考锚点之间拼接 alternate complete nodes，并与 reference internal sequence 比较。变化且长度满足候选范围时写入 locus-aware allele。presence 矩阵第二次扫描 GAF：只有同一条高质量记录完整跨越一个 locus 的两个外侧锚点时，样本才得到 0 或 1；没有完整 spanning record 时为 `NA`，不把未观察等同于缺失。

## 4.4 legacy 拓扑 walker 的解释范围

`--gfa-mode legacy` 从参考节点出发遍历非参考分支并寻找回接点，可用于没有 GAF 的 rGFA 候选提取。该模式不具有可靠的样本路径、位点频率和等位关系语义，因此输出只表示拓扑候选，不用于样本 presence/absence 推断。

# 5. 基于所有组装注释的 GFA 节点身份判定

## 5.1 manifest 数据契约

manifest 为制表符分隔文件，核心列如下。

| 字段 | 含义 |
|---|---|
| `sample` | 必填；与 PanSN 样本名一致 |
| `haplotype` | 单倍型，默认 `0` |
| `te_annotation` | RepeatMasker `.out`、BED 或 TE GFF3 |
| `gff3` | 可同时含 gene、transcript、exon、CDS、UTR 和 TE |
| `seqid_map` | 注释 seqid 到 graph contig 的显式映射 |
| `path_prefix` | 非 PanSN path 的唯一前缀 |
| `gaf` | 该条目专用 GAF；也可使用全局 `--gfa-gaf` |
| `fasta` | 输入溯源和 checksum |
| `annotation_complete` | 该组装是否具有完整 TE annotation |

同一 `(sample,haplotype)` 不得重复。相对路径以 manifest 所在目录解析。所有输入文件写入 provenance，并计算 SHA-256。

## 5.2 注释解析与 feature 分类

GFF3 的 `transposable_element`、`repeat_region`、`mobile_element` 及包含 transposon/retroelement 语义的 feature 被归为 TE；gene、pseudogene、mRNA、transcript、exon、CDS、intron 和 UTR 被归为 gene structure。若 transcript 只有 exon，按同一 parent、同一 seqid 排序相邻 exon，并以 $[exon_i.end,exon_{i+1}.start)$ 派生 intron。

RepeatMasker `.out` 默认将 repeat class 归为 TE，但 `Simple_repeat`、`Low_complexity`、`Satellite`、`rRNA`、`tRNA`、`snRNA`、`scRNA` 和 `srpRNA` 归为 other。无效坐标、空区间和 `end<start` 记录不进入 feature 集。BED 按 0-based half-open 读取。

feature 使用 1 kb bin 建立索引。区间 $[s,e)$ 被加入从整数商 `s // 1000` 到整数商 `(e - 1) // 1000` 的 bin。查询 occurrence 时只检查相交 bin 内的候选 feature，并以精确区间交集复核。

## 5.3 PGGB occurrence 恢复

P 路径优先按 `sample#haplotype#contig` 解析，也支持 manifest `path_prefix`。每个节点 occurrence 保存 genome interval、node-local interval、方向、path ID 和 evidence。若节点路径起点为 $g_i$、节点长度为 $l_i$、与前一节点 overlap 为 $o_{i-1}$，则：

$$
start_i=g_i-o_{i-1},\qquad end_i=start_i+l_i
$$

W walk 使用 `seqStart` 作为初始基因组坐标且相邻 overlap 为 0。PGGB 的 P/W 表示完整组装路径；某样本的 `Absent` 只有在其所有 manifest haplotype 都具有同 contig 的完整内嵌路径且该节点未出现时才成立，否则为 `NA`。

## 5.4 minigraph occurrence 与 gapped CIGAR 变换

rGFA 的 `SN/SO` 为节点来源组装提供 origin occurrence。GAF 进一步提供同一节点在其他组装 traversal 中的坐标。默认使用 MAPQ 不低于 20、`matches/block_length` 不低于 0.90 的 primary GAF。

每条 resolved GAF path 先建立 node layout $[p_i,p_i+l_i)$。CIGAR 从 `path_start` 开始扫描：

- `M/=/X` 同时推进 query 和 path，并将该 match block 与 node layout 求交；
- `I/S` 只推进 query；
- `D/N` 只推进 path；
- `H/P` 不推进。

对正链记录，match block 内 path 偏移 $d$、交集长度 $h$ 映射到：

$$
[qstart+q_{consumed}+d,\ qstart+q_{consumed}+d+h)
$$

对负链记录映射到：

$$
[qend-q_{consumed}-d-h,\ qend-q_{consumed}-d)
$$

若 node step 为正向，node-local 区间是 `[overlap_start-node_left, overlap_end-node_left)`；若 step 为反向，则为 `[node_right-overlap_end, node_right-overlap_start)`。负链 query 还会翻转 occurrence 方向。连续 `M/=/X` 均是 1:1 坐标变换，因此合并为一个 match run；indel 仍作为线性块边界保留。

一个真实 traversal 中的同一 node 可能由 indel 拆为多个 occurrence block。节点支持率按 `(entry,seqid,node_id,node_step,path_id)` context 聚合，而不是按 block 数投票，使 CIGAR 的表示碎片数不影响 TE 支持率和 disagreement。

## 5.5 feature 到 node 的方向变换

对 feature $[f_0,f_1)$ 与 occurrence genome 区间 $[g_0,g_1)$，交集为：

$$
[u,v)=[\max(f_0,g_0),\min(f_1,g_1))
$$

若 $v\le u$ 则无投影。正向 occurrence 的 node-local 坐标为：

$$
[n_0+(u-g_0),\ n_0+(v-g_0))
$$

反向 occurrence 为：

$$
[n_1-(v-g_0),\ n_1-(u-g_0))
$$

反向 occurrence 同时翻转 feature strand。每条 projected feature 保留原 sample、haplotype、seqid、genome interval、feature ID、parent、family、class、source、strand 和 phase。

## 5.6 occurrence、context 与 node 身份

对单个 occurrence block，TE covered bp 是所有 TE node-local 区间的并集长度，fraction 为 covered bp 除以该 block 的 node span。默认 fraction 不低于 0.80 为 TE；不低于 0.10 为 PartialTE；完整 TE annotation 中没有 TE 为 NonTE；不完整 annotation 中没有 TE 为 Unknown。

节点汇总先对所有 projected TE 的 node-local 区间求并集。令并集覆盖的碱基数为 $C_{TE}$，节点长度为 $L_{node}$，则：

$$
F_{TE}=\frac{C_{TE}}{L_{node}}
$$

TE context support 为含 TE 的 callable traversal context 数除以全部 callable context 数；TE genome support 为含 TE 的样本数除以 callable 样本数。节点状态按以下顺序确定：

1. 同一 occurrence block 的 TE 与 CDS 有碱基重叠时为 `Conflict`。
2. 节点同时具有 TE 与 gene structure，或完整注释的 traversal context 对 TE 与非 TE 判断不一致时为 `Mixed`。
3. $F_{TE}\ge0.80$ 且 TE context support 不低于 0.5 时为 `TE`。
4. $F_{TE}\ge0.10$ 或存在任意 TE 投影时为 `PartialTE`。
5. 无 TE 且存在 gene structure 时为 `GeneStructure`。
6. 无 TE/gene structure 且存在 callable context 时为 `NonTE`。
7. 无 callable 证据时为 `Unknown`。

样本状态使用 `TE/PartialTE/NonTE/Unknown/Absent/NA`，这些状态不能合并。特别地，rGFA 中未观察到高质量 traversal 只得到 `NA`，不推断 absence。

## 5.7 输出文件

| 文件 | 数据层级和用途 |
|---|---|
| `gfa_node_annotation_intervals.tsv` | feature 投影明细及原基因组坐标 |
| `gfa_node_occurrence_annotations.tsv` | 线性 occurrence block 的 callability 和 TE coverage |
| `gfa_node_annotations.tsv` | node 身份、覆盖率、支持度、family/class 和宽格式样本状态 |
| `gfa_node_sample_matrix.tsv` | node × sample 长格式状态矩阵 |
| `gfa_node_annotation.gff3` | node-local GFF3 |
| `gfa_node_annotation_conflicts.tsv` | TE–CDS overlap |
| `gfa_annotation_unresolved.tsv` | 无法精确解释的 path/GAF/seqid |
| `pangenome.annotated.gfa` | S 行附加身份、TE fraction、class、family、gene structure 和 confidence |
| `gfa_annotation_provenance.json` | 输入 checksum、阈值、图模式、计数和 GAF 统计 |

# 6. Look4LTRs 结构轨道

`LTR_detect` 在规范化基因组上识别具有长末端重复结构的候选，并输出 `.rtr`。`build_ltr_library.py` 根据图分组和家族归属建立共识。分块结果通过 `cd-hit-est -c 0.8 -aS 0.8` 合并为 `Look4LTRs/consensus_lib/ltr_library.fa`。

该轨道的证据来自 LTR 边界和内部结构，因此优先覆盖结构保存良好且相对年轻的 LTR retrotransposon。结构破损、嵌套或高度分化的候选由 mdl-repeat 和 TE-looker 补充。80% identity 和 80% shorter-sequence coverage 用于轨道内去除明显重复共识，同时保留分化家族。

# 7. mdl-repeat 的 MDL 重复发现

## 7.1 最小描述长度目标

mdl-repeat 使用 Minimum Description Length 原则比较候选家族模型 $H$ 对观测序列 $D$ 的压缩能力：

$$
L_{total}=L(H)+L(D\mid H)
$$

$L(H)$ 描述共识、家族参数和模型本身；$L(D\mid H)$ 描述所有拷贝相对共识的差异、位置和未解释序列。真实重复家族使“一个模型加多个差异编码”短于逐段原样编码；过多的偶然家族会增大 $L(H)$，过度合并分化序列会增大 $L(D\mid H)$。

Pan_TE 使用 mdl-repeat 输出 header 中的 `R`、`length`、`copies`、`mdl`、`tier` 和 `accept`，并使用 instances BED 与 stats TSV 作为坐标和家族证据。mdl-repeat 内部优化器属于独立可执行程序，本文不对未暴露的内部搜索顺序作推断。

## 7.2 低复杂度预掩蔽

`build_mdl` 先运行 `dustmasker`，再以 TRF 参数 `2 7 7 80 10 50 500 -m -h` 掩蔽串联重复。参数依次表示 match reward、mismatch penalty、indel penalty、match probability、indel probability、minimum score 和 maximum period。预掩蔽降低微卫星、同聚体和简单重复对 MDL 压缩收益的支配。

mdl-repeat 在掩蔽序列上输出 raw family FASTA、instances BED 和 stats TSV，三者同时进入 Refiner_mdl。

## 7.3 大基因组确定性采样

输入超过 2 Gb 时，发现和 Phase 1 在确定性窗口样本上运行，Phase 2b 使用完整基因组重新计算 recurrence/copy evidence。采样坐标属于样本 contig；full genome 与 sample 由 genome stamp 区分，防止复用坐标体系不同的掩蔽和检查点。采样提高高拷贝家族发现的可计算性，但发现范围受样本覆盖限制，因此低拷贝家族的未观察不解释为不存在。

# 8. Refiner_mdl 家族精修

## 8.1 Phase 0：解析、硬过滤和分层

Phase 0 读取 MDL header 和辅助文件，应用以下默认门槛。

| 指标 | 默认值 | 数学或生物学作用 |
|---|---:|---|
| length | 至少 50 bp | 排除信息量过低的片段 |
| N fraction | 至多 0.20 | 保证可比较碱基比例 |
| k=1 Shannon entropy | 长序列至少 1.0；短序列至少 0.8 | 排除单碱基偏置 |
| k=2 Shannon entropy | 长序列至少 2.0；短序列至少 1.6 | 捕捉二核苷酸重复 |
| low-complexity window fraction | 至多 0.60 | 排除多数窗口为简单重复的序列 |
| DUST score | 至多 0.70 | 抑制 SSR/simple repeat |
| MDL score | 大于 0 | 要求正压缩收益 |

Shannon entropy 以 $H_k=-\sum_i p_i\log_2p_i$ 计算。k=1 最大为 2 bit；k=2 最大为 4 bit。长度低于 100 bp 的序列使用较低门槛，以容忍有限样本对 k-mer 频率估计的方差。

`mdl_per_copy` 的 P25/P75 将家族分为 T1/T2/T3。分层用于分配计算预算和输出用途，不表示 TE superfamily 分类。

## 8.2 Phase 0 片段共现图与拼接

instances BED 提供家族在基因组上的位置和方向。算法在同 chromosome、同 strand、允许间距内的实例之间统计有向共现。负链实例先翻转家族顺序，使两条链对同一生物学顺序投票。候选家族边必须同时满足：至少 3 个独立共现位点；共现数除以两家族较小 copy 数至少 0.30；主方向一致性至少 0.80；较小/较大 copy 数比至少 0.30；overlap 不超过较短序列的 0.50。

对通过门槛的有向边按证据强度排序，贪婪选择每个节点最多一个左邻和一个右邻，形成低分支链。stitch gap 默认 300 bp，look-ahead window 为 4，最大 stitched length 为 25 kb。链生成 scaffold，Phase 1 从真实跨片段实例重新提取和建共识。无 BED 时可使用 BLASTN 坐标共现作为替代证据。

Phase 0 使用 near-exact CD-HIT-EST identity 0.98 去除重复描述。更松的边界收敛合并在 Phase 1 共识建立之后执行。

## 8.3 Phase 0.5 copy gate

copy gate 位于 MSA 之前。启用 copy recruitment 时，rmblastn HSP 按 query、chromosome 和 strand 合并，实例需要 consensus coverage 至少 0.5。`hard_min_copies` 大于 0 时构成全家族 recurrence floor。结构信号只作为 rescue，不作为所有家族必须满足的条件。

## 8.4 Phase 1 实例提取和边界校正

实例按 strand 统一方向，并在两端加入：

$$
padding=\min(500,\max(50,0.2\times family\ length))
$$

FASTA 提取按最多 200 个 region 的批次执行。每家族至少需要 5 个 copy 进入 MSA，否则保留 MDL seed。确定性 subsampling seed 为 42；MAFFT 路径最多 100 个 copy，abPOA 路径最多 300 个 copy。

MSA 引擎根据 copy 数、分化和总 bp 选择。低分化小家族可使用 MAFFT L-INS-i；常规家族使用 MAFFT auto；copy 数或总 MSA bp 较大时使用 abPOA。单序列超过 50 kb、总输入超过 8 Mb 或每次 MSA 超过 60 s 时不生成推测性共识，而是保留原 seed。subfamily consensus length 大于 2,500 bp 或 `copy_count × consensus_length` 超过 500,000 时不进入昂贵的 subfamily MSA。

对齐列 occupancy 为该列非 gap copy 数。两端连续低 occupancy 列被裁剪；默认相对 floor 为 0.30，且至少需 3 个 copy 支持。共识由保留列的多数碱基形成。若对齐、occupancy 或边界证据不足，输出原 seed，形成 never-regress fallback。

## 8.5 嵌合与亚家族

MSA 嵌合检测寻找内部 occupancy valley。valley occupancy 低于中位 occupancy 的 0.5、连续至少 3 列，并要求左右 copy span group 各占至少 0.20；允许的 full-span copy fraction 至多 0.50。该 AND gate 区分内部连接断点与随机低覆盖列。

copy clustering 位于实例提取和共识建立之间。局部 BLASTN identity 构成主要 pairwise distance；大 copy set 可采用 k-mer Jaccard 粗筛。满足 cluster 最小规模和 coherence 的子簇分别生成 consensus，均质家族产生一个 cluster。嵌合拆分解决 A+B 错接，亚家族聚类解决同一元件内部的进化分支，二者不互相替代。

## 8.6 Phase 2 质量控制和双库输出

Phase 2 对 refined consensus 再次执行长度、N fraction、熵、低复杂度和 DUST 检查，因为边界延伸和拼接可重新引入低复杂度尾部。

TE protein、LTR/TIR 和已知 TE homology 默认作为注释，不作为普遍硬删除条件。MDL quality gate 保留满足任一条件的家族：tier 为 `core`；存在 TE signal；copies 至少 50。只有低 tier、无 TE signal 且低 copy 的交集被排除。缺少 tier 字段或无法获得结构信号时，该门控不利用缺失证据作负判断。

masking library 使用长度依赖 copy floor：长度至少 500 bp 需 2 copies；200–499 bp 需 3；100–199 bp 需 5；低于 100 bp 需 10。analysis library 包含 masking library，并加入 copies 至少 3 的 T3。输出分别为 `consensus_masking.fa` 和 `phase3_analysis_library.fa`。

# 9. TE-looker 生产轨道

## 9.1 capability probe 与自适应 seed floor

Pan_TE 从 `--te-looker-bin`、`TE_LOOKER_BIN`、PATH 和仓库 `bin/dtr` 解析可执行文件，并查询 `--version` 与 `dtr run --help`。capability probe 决定 seed-harvest 内存副本数、`--final-gated` 是否为显式能力以及 provenance 的接纳方式。

canonical 16-mer 的基础 `--min-count` 随基因组大小取值：小于 300 Mb 为 50；300 Mb–1 Gb 为 100；1–3 Gb 为 200；3–10 Gb 为 300；大于 10 Gb 为 500。若用户显式使用 $w>1$，门槛按 minimizer 密度缩放：

$$
C(w)=\max(5,\operatorname{round}[C(1)\rho(w)])
$$

该缩放保持“一个真实家族需要约多少观察种子”的语义。`RUST_MIN_STACK` 设为 64 MiB，降低并行 seed harvest 的线程栈溢出风险。

## 9.2 HMM 种子

Refiner_mdl 的 masking consensus 可通过 MAFFT 对齐和 `hmmbuild` 构建 nucleotide HMM。该 HMM 作为已知/高置信家族的 profile 证据输入 TE-looker Track 1。HMM 构建失败不生成伪模型；TE-looker 仍可使用 de novo discovery 轨道。

## 9.3 gated contract 和宿主门控

满足 `te-looker-gated-v1` 的 dtr provenance 表示 `families.fasta` 包含 fragment stitch、copy/length gate 和跨 seed merge 的门控结果，Pan_TE 对文件和契约作一致性验证后使用该库。

不满足 gated contract 时，Pan_TE 对 discovery family 执行：

1. 使用 `disc.members.bed` 统计不同 fragment family 在多个实例中的共线邻接，重提取支持链并建立 stitched consensus；拼接无可信证据时保留原 family。
2. 从 header 读取 `members=N`，保留成员数至少 5 且长度至少 100 bp 的 family；缺少 members 字段时按 5 处理。
3. 使用 `cd-hit-est -c 0.80 -aS 0.60 -n 5 -r 1 -g 1` 合并跨 seed 的部分重叠 family。
4. 最终 family 数不得超过 `15 × genome_Mb`；超过该上限表示 discovery inflation 未被门控控制，流程停止。

`pan_te_handoff.json` 记录输入、发现库、拼接库、门控库和最终轨道库的 SHA-256、序列数、总 bp、阈值和拼接状态。

## 9.4 嵌套实例分解与自适应边界重估

该 refinement 以显式命令作用于 gated family library。`dtr-copy` 首先生成全基因组 copy BED。`dtr-nested --min-flank 80` 在同 chromosome 上识别一个 family 的 inner instance 被另一 family 的 outer instance 完整包含、且 outer 两侧各保留至少 80 bp 的事件。含 inner 的 outer instance 不参与 outer family 边界估计；inner family 仍由自身 copy 表示。

`dtr-boundary` 从 clean instances 提取默认 100 bp flanks，统一方向后以 POA 生成 proposal。family 至少需要 5 个成员；proposal length 必须位于原长度的 $[1/1.5,1.5]$ 倍范围；单 family POA 输入默认最多 4,000 bp。预算不足或证据不足时原 consensus 原样保留。

proposal library 重新运行 `dtr-copy`。`dtr-boundary-gate` 对每个 family/chromosome 的实例区间求并集，只有 proposal coverage 相对原库增加至少 1% 且序列发生改变时采用 proposal，其余 family 回退到原 consensus。coverage 是边界 proposal 的选择代理，最终库效果由匹配参数的独立 RepeatMasker 验证。

## 9.5 Track 2 部分证据链的当前输出语义

`dtr-window` 从重复 canonical spaced seed occurrence 中按全局 hash-bottom 选择候选窗口；默认 seed pattern 为 `11101110111011101111`，span 20、weight 16。`dtr-lsh` 对窗口计算 canonical 15-mer bottom-k MinHash（128 hashes），按 band 生成候选对，每个 bucket 最多 128 个节点，精确 Jaccard 至少 0.08，每节点最多保留 16 条边。

`dtr-community` 输出连通分量容器，`dtr-community-refine` 执行确定性 first-level Louvain local moving。该社区结果不是 Leiden 多层优化，也不直接构成 TE family。`dtr-component-poa` 只对成员数满足门槛的社区生成 provisional consensus。

`dtr-copy` 使用 BLASTN 建 copy catalog：HSP 按 query/chromosome/strand 分组，query gap 至多 200 bp、genome gap 至多 500 bp 时合并；合并实例要求 query coverage 至少 0.60、identity 至少 75%；同 query/chromosome/strand 上起点相距不超过 100 bp 的冗余实例只保留一个。

`dtr-nonte-guard` 使用用户明确提供的 non-TE reference。identity 至少 80% 的 HSP 在 query 坐标上求并集；high-identity query coverage 至少 0.80 时，候选被解释为 non-TE。`dtr-structure` 审计 copy flanks 中 3–15 bp TSD、至少 12 bp terminal direct/inverted seed 和 terminal poly-A/T；少于 20 copies 时标记 boundary uncertain。

`dtr-accept-audit` 分列报告 copy、structure、boundary、non-TE、segmental-duplication 和外部守卫状态。Track 2 provenance 固定 `mode=partial`、`family_call=false`、`acceptance_complete=false`、`final_library=null`，因此这些 provisional consensus 不进入 Pan_TE 的生产最终库。

# 10. 三轨合并、去冗余与可选分类

合并输入为 Look4LTRs `ltr_library.fa`、Refiner_mdl `phase3_analysis_library.fa` 和 TE-looker `consensus.fa`。三者直接串接为 `raw_te_library.all.fa`，再执行：

```text
cd-hit-est -c 0.90 -aS 0.80 -aL 0.80 -G 0 -n 8 -r 1 -mask NX
```

`-c 0.90` 要求高序列 identity；`-aS/-aL 0.80` 要求两条序列均有充分覆盖；`-G 0` 使用局部 identity 定义；`-r 1` 同时考虑反向互补。跨轨阈值高于 TE-looker 轨道内 0.80，是因为跨轨合并只消除近重复，不把相关但独立的家族合并。

去冗余未分类库为 `combine/raw_te_library.fa`。提供 `--model-dir` 时，RepeatClassifier 与 ClassifyTE 在隔离环境中运行，并由 `Combine_for_Two` 合并标签，输出 `combine/TEs.fa`。分类结果与发现证据分层保存；分类缺失不改变序列库的重复家族身份，但不能把 Unknown family 自动解释为特定 TE class。

# 11. 参数选择原则与敏感性

| 参数 | 默认值 | 直接影响 | 调整后的解释 |
|---|---:|---|---|
| PAF minimum alignment | 5 kb | 共线锚点稳定性 | 降低会增加局部重复 HSP |
| PAF MAPQ | 20 | 多重比对过滤 | 降低会增加错误链 |
| PAF insertion length | 100–30 kb | 候选尺度 | 范围外事件不进入发现集 |
| PAF gap ratio | 0.30 | 缺口不对称程度 | 增大会提高召回并增加 SV 噪声 |
| GAF MAPQ | 20 | rGFA traversal | 降低会扩大错误 occurrence/presence |
| GAF identity | 0.90 | rGFA path 可信度 | 降低适用于远缘样本但需审计 |
| GFA long-allele overlap | 1 kb | 长候选窗口连续性 | 增大增加重复窗口 |
| annotation TE fraction | 0.80 | TE node/occurrence 判定 | 代表节点大部分被 TE 覆盖 |
| annotation partial fraction | 0.10 | PartialTE 判定 | 捕捉跨节点 TE 边界 |
| TE-looker w | 1 | seed density和内存 | 大于 1 改变 seed set |
| TE-looker min members | 5 | family recurrence | 提高增加特异性、降低低拷贝召回 |
| TE-looker min length | 100 bp | family信息量 | 提高会排除短元件 |
| final CD-HIT identity | 0.90 | 跨轨 family 合并 | 降低可能吞并亚家族 |

参数敏感性分析应一次只改变一个具有明确统计含义的参数，并保留输入 checksum、工具版本、完整命令、provenance 和输出计数。`w`、MAPQ、identity、copy floor 和 family merge identity 会改变搜索空间或家族定义，不能仅视为计算性能参数。

# 12. 质量控制、溯源与复现

## 12.1 必需审计

每次正式运行至少保存：原始 FASTA/PAF/GFA/GAF/VCF 的 SHA-256；manifest；Pan_TE 和外部工具版本；完整命令；线程和内存上限；每阶段 provenance/stats；最终 FASTA sequence count 和总 bp；CD-HIT cluster 文件；GFA unresolved/conflict 表；分类输入与输出。

FASTA 最低 QC 包括空文件、重复 ID、非法字符、长度范围、N fraction、序列数和总 bp。坐标文件最低 QC 包括 0-based/1-based 契约、start/end 合法性、contig 是否存在、长度是否越界和 strand 值。

## 12.2 GFA 注释验证

投影验证应在 genome-coordinate 层完成。对每个 `sample#hap#contig`，先合并 occurrence 表的 callable intervals，再独立合并源 TE intervals，并计算：

$$
E=(source\ TE\ union)\cap(callable\ union)
$$

将 $E$ 与 interval TSV 中 projected TE 的 genome-coordinate union 逐区间比较。该比较排除 GAF 未覆盖区，适用于 PGGB 和 minigraph。node-local TE fraction 是图中唯一节点序列的覆盖率，不等于任何线性基因组的 RepeatMasker mask rate。

## 12.3 重复库验证

库版本比较使用同一基因组、同一 RepeatMasker 版本、同一参数和同一线程配置。报告 masked bases、genome mask fraction、相对变化、绝对百分点变化、family 数、库总 bp 和运行时间。mask rate 增加只表示该库解释了更多基因组区间，不单独证明新增区间均为 TE；需要结合 family recurrence、non-TE 守卫和编码区 overlap 审计。

# 13. 拟南芥验证数据摘要

## 13.1 GFA 注释

PGGB Chr1 图含 2,270,366 nodes、3,144,225 个 `0M` links 和 7 条 P paths。7 套 RepeatMasker 输入共 800,587 features，恢复 17,011,115 occurrences，投影 12,150,687 intervals，unresolved 为 0。节点状态为 TE 853,457、PartialTE 10,872、Mixed 33,237、NonTE 1,372,800。7 个样本的 projected TE genome union 与源 Chr1 TE union 逐 bp 一致。

minigraph `Ara.gfa` 含 74,660 nodes、总 node sequence 166,090,619 bp。1,376 条 primary GAF 均通过结构验证；默认质量门槛投影 521 条记录，其中 295 条 coordinate walk、226 条 named path，855 条质量不足记录不进入 occurrence。结果含 868,506 CIGAR match runs、1,179,408 occurrence blocks、1,133,307 projected intervals，unresolved 为 0。节点状态为 TE 24,596、PartialTE 17,621、Mixed 5,383、NonTE 27,060。node-local TE union 为 66,307,931 bp，占唯一 node sequence 的 39.9227%。

minigraph 独立 genome-coordinate 审计得到可映射 TE：Col-0 T2T 41,719,175 bp、TAIR10 28,571,472 bp、Cvi-0 25,520,542 bp、Ey15-2 26,638,897 bp、IP-Met-6 26,671,318 bp、Rabacal-1 26,251,204 bp、Tsu-0 27,629,146 bp，总计 203,001,754 bp。`source TE union ∩ callable union` 与 projected TE union 完全一致。

## 13.2 TE-looker 生产库和边界 refinement

拟南芥 gated TE-looker 运行产生 337 个发现/精修 families，共 673,452 bp。fragment stitch 得到 289 families、673,719 bp；成员数和长度门控保留 211 families、326,477 bp；轨道内 CD-HIT-EST 形成 69 families、167,509 bp。

69-family library 的 copy catalog 含 2,307 instances；嵌套分解识别 176 个跨 family nested events，clean copy 集含 2,242 instances。边界 proposal gate 接纳 2 个 family 的 proposal。匹配参数 RepeatMasker 的 masked bases 为 23,018,188 bp，mask fraction 为 17.188446%；对应未采用 proposal 的 69-family library 为 23,012,842 bp、17.184454%，绝对增加 5,346 bp 和 0.003992 percentage points。

# 14. 结果解释边界

1. 从头重复 family 表示可重复序列模型，不自动给出 TE superfamily。分类模型缺失时应报告为 de novo repeat library。
2. PAF 不对称 gap 是候选插入证据，不是 TE 判定；复杂重排、assembly gap 和重复比对可产生候选。
3. PGGB allele 依赖完整 P/W path 和稳定参考锚链；rGFA absence 依赖完整 spanning GAF，未观察为 `NA`。
4. GFA node 是共享图序列单位，同一 node 可处于多个基因组和多种注释 context；node identity 是 context 汇总，不替代 occurrence-level 明细。
5. query-only GAF insertion 没有 graph node，不能产生 node annotation；它仍可由 PAF/GFA allele extraction 作为发现候选处理。
6. 大基因组采样保留高拷贝家族的发现能力，但低拷贝 family 的发现范围受样本覆盖限制。
7. TE structure、protein domain 和已知 homology 在 Refiner_mdl 中默认是证据标记或 rescue，不构成所有 family 的必要条件。
8. Track 2 输出是部分证据审计，`family_call=false`，不与生产最终库混用。
9. mask rate 是库对基因组的覆盖指标，不单独等于 TE precision、family classification accuracy 或生物学完整性。

# 15. 论文 Methods 写作映射

| 论文 Methods 小节 | 本文档对应内容 | 建议报告的可复现信息 |
|---|---|---|
| Input assemblies and preprocessing | 第 1–2 节 | assembly version、FASTA checksum、清洗规则、坐标体系 |
| Pangenome insertion recovery | 第 3–4 节 | PAF/GFA/GAF 生成命令、MAPQ、identity、length 和 gap thresholds |
| De novo repeat discovery | 第 6–7 节 | Look4LTRs、Dust/TRF、mdl-repeat 版本和参数 |
| Consensus refinement | 第 8 节 | copy source、padding、MSA engine、occupancy、chimera 和 copy thresholds |
| Twilight-zone repeat discovery | 第 9 节 | dtr version、w、min-count、memory、gated contract 和 family gates |
| Library integration and classification | 第 10 节 | 三轨输入计数、CD-HIT 参数、分类模型版本 |
| GFA functional annotation | 第 5 节 | manifest、annotation completeness、CIGAR transform 和 node-state thresholds |
| Validation and reproducibility | 第 11–13 节 | provenance、unresolved、mask rate、callable-intersection audit |

# 16. 代码与输出溯源索引

| 功能 | 实现位置 |
|---|---|
| 主调度、资源模型、TE-looker handoff、三轨合并 | `bin/Pan_TE` |
| PAF 解析、链、gap/CIGAR 候选和 QC | `bin/paf_processor.py` |
| PGGB/rGFA 图发现入口 | `bin/gfa_processor.py` |
| PGGB locus-aware allele 聚合 | `bin/gfa_alleles.py` |
| strict rGFA/GAF 锚链、allele 和 presence | `bin/rgfa_gaf.py` |
| 全基因组注释到节点的投影 | `bin/gfa_annotator.py` |
| mdl-repeat wrapper 和预掩蔽 | `bin/build_mdl` |
| Refiner_mdl phase orchestration | `bin/Refiner_mdl/main.py` |
| fragment triage/stitch | `bin/Refiner_mdl/phase0_triage.py`、`bin/fragment_stitcher.py` |
| copy gate | `bin/Refiner_mdl/phase05_copy_gate.py` |
| copy extraction、boundary、MSA、consensus、chimera | `bin/Refiner_mdl/phase1_*.py` |
| final QC、copy recruit、non-TE split、双库 | `bin/Refiner_mdl/phase2_*.py` |
| 输入环境预检 | `bin/preflight.py` |
