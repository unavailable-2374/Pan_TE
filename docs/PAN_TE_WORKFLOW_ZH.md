# Pan_TE 完整工作流、算法与参数说明

> 适用代码：本仓库当前 `bin/Pan_TE` 及其调用脚本。本文以**当前可执行实现**为准；对尚未由当前 `dtr` 二进制完整实现的 TE-looker 设计，会明确标为“目标架构”，避免把设计文档当作已执行事实。

## 1. 目的、输入模式与总体策略

Pan_TE 是泛基因组转座元件（transposable element，TE）从头发现与可选分类流程。它把三个互补的发现轨道合并为一个非冗余文库：

| 轨道 | 覆盖重点 | 主要产物 |
|---|---|---|
| Look4LTRs | 保守结构明显、尤其较年轻的 LTR 逆转座子 | LTR 家族共识 |
| mdl-repeat + Refiner_mdl | 广谱、拷贝驱动的从头重复序列发现 | 经拷贝和序列质量校正的共识 |
| TE-looker (`dtr`) | 与已知家族远缘、破碎或处于“twilight zone”的重复序列 | 发现级或全门控 TE 家族 |

三条轨道不是互相替代的关系。Look4LTRs 用结构约束提高年轻 LTR 的特异性；mdl-repeat 用压缩/重复性信号提供广泛种子；TE-looker 试图补足较深分化、在传统相似性阈值附近变得难以聚类的家族。最终以较严格的跨轨道去冗余合并，避免将相关但不同的家族过早压成一个。

支持四类输入情景：

1. 单基因组：`--genome genome.fa`。
2. PAF 泛基因组：`--paf all_vs_all.paf --paf-fasta asm1.fa asm2.fa ...`。从装配间不对称缺口提取插入序列，构造供三个发现轨道使用的富集序列集。
3. GFA 泛基因组：`--gfa graph.gfa`。PGGB 与 rGFA 走不同的、坐标可追溯的等位序列提取路径，详见第 7 节。
4. 富集模式：在 `--genome` 已给定参考基因组时同时给 PAF、GFA 或 `--vcf-dir`；发现引擎分析“清洗后的参考 + 新发现的候选插入”而非仅参考。

入口程序为 `bin/Pan_TE`。默认输出目录为当前目录，也可用 `--out-dir` 指定。所有坐标性内部处理使用 0-based、half-open 区间 `[start,end)`；仅调用 `samtools faidx` 时转换为其 1-based、闭区间写法 `start+1-end`。

## 2. 端到端执行图

```text
输入 FASTA / PAF / GFA / VCF
          │
          ├─ 环境预检（工具、模块、资产、dtr 内存）
          │
          ├─ 基因组准备
          │    ├─ 清洗、重命名、faidx 与索引
          │    └─ PAF/GFA/VCF 候选序列富集（可选）
          │
          ├──────────── 并行 ────────────┐
          │                              │
          │  Look4LTRs                    │  mdl-repeat
          │  LTR_detect → 共识构建         │  Dust/TRF → MDL 发现 → Refiner_mdl
          │                              │
          └───────────────┬──────────────┘
                          │
                    TE-looker / dtr
                    （可由 mdl 共识构建 HMM 种子）
                          │
                  片段拼接、门控、轨道内去冗余
                          │
                  三源合并 + CD-HIT-EST 0.90
                          │
             RepeatClassifier + ClassifyTE（可选）
                          │
                 GFA 节点注释与样本存在矩阵（可选）
```

主程序只在环境预检通过后建立工作目录。`ltr` 与 `mdl_repeat` 使用两个 `ThreadPoolExecutor` 任务并行；TE-looker 依赖 mdl-repeat 的高质量共识以构建可选 HMM，因此在二者完成后执行；随后才合并与分类。阶段检查点由 `.ok` 文件控制：`genome.ok`、`ltr.ok`、`mdl_repeat.ok`、`te_looker.ok`、`combine.ok` 和可选 `anno.ok`。GFA 模式另以 `gfa_run_fingerprint.json` 记录输入文件大小、mtime、模式、GAF 及阈值；配置或输入改变时拒绝复用旧 `genome.ok`，防止混用不兼容候选集。

## 3. 阶段 0：环境预检与基因组规范化

### 3.1 预检

`--check-env` 只执行检查而不写分析结果。预检汇总而不是遇到首个错误即退出，覆盖：外部程序、Perl/Python 导入模块、可选 ClassifyTE 模型和 Dfam/RepeatMasker 资产，以及 TE-looker 的内存可行性。这样可避免运行数小时后才发现一个环境缺口。

TE-looker 的发现核心在默认精确模式 `w=1` 下为每个 k-mer 位置保留一个 16-byte 打包条目；估算峰值为：

```text
peak ≈ 1.15 × [genome_bp + genome_bp × 16 × density(w) × copies]
density(w) = 1                         (w ≤ 1)
           = 2 / (w + 1)               (w > 1，期望 minimizer 密度)
```

旧版（无可靠 `--version` 或 < 0.3.1）在拼接排序时可能同时保留两份条目，`copies=2`；新版为 1。可用内存取 `-M`、cgroup、SLURM 和 `/proc/meminfo` 中最小的硬上限。流程不会因内存较小而暗中调大 `w`，因为这会改变种子集合，使结果依赖于机器；用户必须显式使用 `--te-looker-w`。

### 3.2 基因组清洗

`clean_seq_fast` 流式处理 FASTA：序列转大写，并将退化碱基（如 `W,R,M,K,Y,S,H,B,V,D,X`）替换为 `N`；contig 随后统一为 `chr1...chrN`。清洗后的文件写为 `genome/genome.fa`，然后执行 `samtools faidx` 和项目索引程序 `index`。

其理由是让坐标、后续 `faidx` 提取、不同输入来源的字符集和重复检测行为一致。代价是原始 contig 名称不再可作为内部名称使用，因此 PAF 参考名称在清洗前会被捕获，供富集阶段排除参考候选；GFA/PAF 的原始映射均应保存在各自的 manifest 或统计文件中。

### 3.3 VCF 富集（可选）

`--vcf-dir` 接受 `.vcf` 和 `.vcf.gz`。仅将可得到实际核苷酸序列的插入加入发现输入：显式 ALT 序列，或符号插入中 `SVINSSEQ`、`INSSEQ`、`SEQ`、`CONSENSUS` 标签提供的序列。没有序列的 `<INS>` 不可用于从头发现，故会跳过而不是猜测；这样保证 TE 库始终来自可审核的碱基序列。

## 4. Look4LTRs 轨道

`process_ltr` 对 `genome/genome.fa` 调用 `LTR_detect`。后者封装 Look4LTRs C++ 检测器；输入文件大于 5 GB 时按既定规模策略分块。所有 `.rtr` 在一次 `build_ltr_library.py` 调用中联合处理，文件来源加入内部家族键，使分块内局部 element ID 和 GraphGroup ID 不会碰撞，同时避免对同一完整基因组串行重复执行 FastGA。

`.rtr` 完整元件按 LTR identity、PPT 和 TSD 分层：identity≥0.85 且同时有 PPT/TSD 为 Tier 1；identity≥0.70 且有 PPT 为 Tier 2；其余完整元件为 Tier 3；缺少右 LTR 的记录为 solo。GraphGroup 是局部候选家族，不被解释为跨染色体的全局家族。builder 以基因组 FASTA 文件≥250 MB或原始 GraphGroup≥2500作为 scalable mode 触发条件；阈值分别可由 `--scale-genome-mb` 和 `--scale-family-count` 调整。小输入使用 sensitive mode，保留全部 GraphGroup；scalable mode 保留含至少两个完整元件、或至少一个 Tier 1/2完整元件的 GraphGroup。solo仍可在合格组中贡献LTR序列，但solo-only组不能定义内部区；缺乏第二拷贝和高等级结构证据的Tier-3单例只在scalable mode中门控。`--include-low-confidence-singletons` 可在scalable mode恢复全组输入。

序列提取分别处理左 LTR、右 LTR和内部区。已知嵌套元件的坐标从内部区中扣除，剩余片段按原顺序拼接，因此内部共识不会直接吸收已报告的嵌套插入。每个家族最多使用 identity 最高的 20 个元件（最多 40 条左右 LTR 序列）及 10 条内部序列；MAFFT 后按 LTR identity 加权逐列投票，非 gap 权重大于总权重的 0.5 才输出该列。初始完整共识为 `LTR_consensus + interior_consensus + LTR_consensus`；没有可用内部区时仅输出 LTR 共识。

scalable mode默认在全基因组校正前以`cd-hit-est -c 0.8 -aS 0.8 -d 0`对初始共识去冗余；`-d 0`保留完整家族ID，聚类成员与FASTA ID必须精确闭合。sensitive mode默认跳过这一固定开销；`--precluster`和`--no-precluster`只覆盖预聚类选择，不改变scale mode、证据门控或FastGA查询选择。每个代表首先继承Look4LTRs已报告的完整元件坐标。sensitive mode将全部代表送入FastGA，以保持小基因组中的跨GraphGroup连接；scalable mode只搜索少于两个非重叠已知loci的代表。sensitive mode的FastGA阈值为identity≥0.60、alignment length≥200 bp、query coverage≥0.30；scalable mode相应为0.70、300 bp和0.50。三项均可由命令行覆盖。每个query只保留按`identity × query_coverage`排序的前30个命中，PAF在流入时即限流，因而内存、磁盘和后续区间复杂度均受`30 × query数`约束。只有FastGA序列比对命中参与50%重叠的family union；`.rtr`完整区间可因嵌套或相邻而重叠，只用于copy support与序列提取，不参与union。默认硬超时为7200 s；运行期间每60 s输出心跳。发生超时时，流程保留`.rtr`已知拷贝证据继续完成，无法取得第二个非重叠locus的低支持单例在copy QC中淘汰，而不是让整个LTR轨道失败。

共享基因组区间重叠达到较短区间的 50% 时，代表通过 union-find 合并。每个合并组最多保留 30 个去重基因组 loci 作为 copy support；校正 MSA 从这些 loci 的坐标序列中确定性均匀抽取最多 10 条，避免长 LTR 家族使 MAFFT 规模失控，而 header 中的 copy 数仍按全部 loci 计算。可用拷贝不足两个时回退到组内最长初始共识。最终过滤长度<200 bp、N比例>0.25或基因组拷贝数<2的序列，再以 80% identity/80% shorter coverage 去冗余，输出 `Look4LTRs/consensus_lib/ltr_library.fa`。`--no-extend` 使用 `.rtr` 已知完整元件数执行相同 copy gate，不会因缺少 Phase 6 header 而把全部家族误删。

FastGA PAF、参数指纹和基因组索引指纹默认保存在 `work/`；只有 query 哈希、基因组路径/大小/mtime及阈值完全一致时才复用。大型逐家族 MSA 目录在成功后清除；`--clean-work` 可显式删除所有断点。Pan_TE 将 builder 的阶段日志实时转发到主日志，因此索引、比对、合并和校正进度均可观察。

这一轨道依赖 LTR 边界及内部结构信号，因而优先覆盖完整、相对年轻的 LTR 元件；结构腐蚀、严重截断或高度分化家族由 mdl-repeat 和 TE-looker 补充。

## 5. mdl-repeat：MDL 核心原理与运行过程

### 5.1 Minimum Description Length（最小描述长度）原理

mdl-repeat 的核心是 MDL 原则：对候选模型/家族 (H) 和观测序列 (D)，比较总编码长度：

```text
L_total = L(H) + L(D | H)
```

其中 `L(H)` 是描述一个重复家族模型、其共识和必要参数的成本，`L(D|H)` 是用该模型解释所有拷贝后剩余误差、突变、未解释序列和位置的成本。若把若干相似片段表示为“一个家族 + 拷贝差异”比逐段原样编码更短，该模式就有重复家族证据。它自然体现两项相互制衡：

- 把每个偶然相似片段各建一个家族会使 `L(H)` 很大，抑制过拟合；
- 把真实但有一定分化的拷贝纳入同一模型会明显降低 `L(D|H)`，支持发现广谱重复。

本仓库只包装独立的 `mdl-repeat` 可执行文件，而不包含其内部优化器源码；因此本文不虚构其未公开的确切搜索顺序或打分公式。Pan_TE 可直接利用的是它输出 FASTA header 中的 `R`、`length`、`copies`、`mdl`，以及较新版本的 `tier=core/warn/reject`、`accept=...`；实例 BED 和统计 TSV 则给出坐标和家族层面的辅助证据。

### 5.2 预掩蔽与调用

`bin/build_mdl` 的顺序为：

1. 对分析基因组运行 `dustmasker -outfmt fasta -out ...`，硬掩蔽局部低复杂度。
2. 将 FASTA 依线程分块，分别运行 `trf input 2 7 7 80 10 50 500 -m -h`，合并为 TRF 掩蔽的 FASTA。
3. 在掩蔽结果上运行 `mdl-repeat`，同时要求 `-instances mdl_repeat.instances.bed` 和 `-stats mdl_repeat.stats.tsv`。
4. 将原始 `mdl_repeat.raw.fa`、BED、stats 交给 `Refiner_mdl/main.py`。

Dust 与 TRF 的作用不是删除 TE，而是尽早压制多聚体、微卫星和简单重复造成的廉价“伪压缩”。这些模式可在 MDL 中产生大量短而高拷贝的假家族，后续再过滤成本更高、也会污染拷贝招募。TRF 参数采用经典重复模式设定：匹配奖励 2、错配/indel 罚分 7、匹配概率 80、indel 概率 10、最小评分 50、最大周期 500；对应重点是串联重复而非一般散在 TE。

### 5.3 大基因组采样路径

当输入大于 2 Gb，`build_mdl` 切到确定性窗口采样：mdl-repeat 和 Refiner Phase 1 都在样本上运行，完整基因组仅在最终 Phase 2b 上用于全基因组复制数/复发门槛。设计依据是 mdl-repeat 面向高拷贝家族：例如全基因组有 100 个拷贝的家族，在约 5% 样本中仍期望留下约 5 个观测。采样的 BED 坐标属于样本 contig，不能与全基因组坐标混用；脚本以 `basename|size` genome stamp 识别全基因组/样本切换并使旧掩蔽、MDL 结果和检查点失效。

这是一种计算折衷，不代表样本中未见的低拷贝家族不存在。大基因组结果应将“发现范围由样本限制”作为解释限制，并检查 `phase2b` 的全基因组计数审计。

## 6. Refiner_mdl：为什么必须 refine、每相的算法与门槛

原始从头重复发现的“家族”通常仍混合短片段、低复杂度、嵌合、边界过伸和不同亚家族；直接拿它掩蔽会使重复库膨胀，并将许多主机序列误标为 TE。Refiner_mdl 的总目标是：保留有可重复支持的真实家族，同时把共识的边界和代表性提升到可用于后续 HMM、RepeatMasker 与分类的程度。

当前调度为：**Phase 0 → Phase 0.5 copy gate → Phase 1 → Phase 2 →（大基因组）Phase 2b**。每相有 pickle 检查点，临时文件默认清理，可用 `--keep-temp` 保留。

### 6.1 Phase 0：解析、硬过滤、分层、片段重组与近似去冗余

`phase0_triage.py` 解析 MDL header 并保留 `R/copies/mdl/length/tier/accept` 元数据。硬过滤包括：

| 规则 | 默认值 | 设计原因 |
|---|---:|---|
| 最短长度 | 50 bp | 低于此值难形成可靠家族模型，随机/低复杂度命中占比高 |
| 最大 N 比例 | 0.20 | 大量未知碱基不可提供同源或结构信息 |
| 单核苷酸熵 | 长序列 1.0；短序列 0.8 | 过滤单一碱基偏置；短序列采样误差更大，阈值稍放宽 |
| 二核苷酸熵 | 长序列 2.0；短序列 1.6 | 捕捉单碱基熵看似正常的 `(CA)n` 等二核苷酸重复 |
| 低复杂度窗口比例 | ≤0.60 | 防止平均 DUST 掩盖“多数为简单重复”的序列 |
| DUST | ≤0.70 | 将 SSR/simple repeat 从散在 TE 候选中排除 |
| MDL 分值 | >0 | 无压缩收益的候选没有 MDL 支持 |

接着以 `mdl_per_copy` 的分位数（P25/P75）把家族分为 T1/T2/T3。分层不是生物学分类，而是分配计算预算：证据更强的 T1 可投入更完整 MSA，T2 标准处理，T3 采用轻量方案且只在满足分析库门槛时输出。

**片段重组。** mdl-repeat 可能将同一真实元件的保守块分别报为不同家族。Phase 0 首选 BED 共现：在同染色体、同链、有限间距的实例间建立有向边；负链时翻转元素顺序，使双链拷贝投票给相同方向。边需要同时通过共现次数/比例、方向一致性和两家族拷贝数比门槛；再按证据强度贪婪选择每节点最多一个左邻和一个右邻，构成无过度分支的链。若无 BED 可退回 BLASTN 坐标共现。对每条可信链构建 scaffold，Phase 1 会重新提取跨片段的真实拷贝并重建全长共识。

最后进行两轮 CD-HIT-EST（95% 后 90%）。先紧后松的目的，是先消掉几乎相同的发现种子，再合并仍明显重叠的家族，且不把远缘亚家族在 refine 前压扁。

### 6.2 Phase 0.5：全家族 copy gate

`phase05_copy_gate.py` 位于昂贵 Phase 1 之前。它可通过 rmblastn 招募真实基因组拷贝（`--enable-copy-recruit`），以覆盖度 `--copy-recruit-cov`（默认 0.5）把 HSP 合并为实例，并可施加 `--hard-min-copies`。在大基因组采样模式，严格全基因组复发地板延后至 Phase 2b，避免样本计数误删真实高拷贝家族。

此相的设计原因是让“是否值得为一个家族做 MSA”取决于基因组复发/连贯性而不只是 MDL header。TE 结构证据可作 rescue，但不应要求每个远缘 TE 都必须带可检测蛋白结构域或 LTR/TIR，因为那会系统性丢失 dark-matter TE。

### 6.3 Phase 1：按拷贝进行边界校正与共识建立

Phase 1 以 BED 实例为主，按 `shard_size` 分片并行；实例索引只流式路由到对应 shard，避免把全部 BED 复制给每个进程。家族按长度 × 拷贝数优先，优先完成计算昂贵且影响最大的任务。每个家族的核心步骤是：

1. 将实例按链方向统一，使用 `samtools faidx` 批量取出带 padding 的序列。padding 由家族长度按比例计算，并受最小值与上限约束；它为边界外延留出证据，而不无限扩张到邻近基因。
2. 实例不足时，可延迟建立一次全基因组 BLAST 数据库并用 `rmblastn` 补招募；默认不会在所有家族上盲目建库。可选 `--enable-rmblastn-copies` 使用 rmblastn 拷贝作为 Phase 1 输入。
3. 根据拷贝数、平均分化和长度选择 MSA 策略；常规路径使用 MAFFT（`--auto` 或 T2 的 `L-INS-i`），过大/不宜 MSA 的情况可用 abPOA。对齐列用 occupancy profile 判断共有边界，低占有端部被裁剪；多数字符形成共识。
4. 若对齐质量不足或低占有，保留原始种子而不是输出貌似精细但不可信的共识。这是显式的保守失败策略。
5. 嵌合检测在共识/对齐深度上寻找内部覆盖谷、跨度/双峰等证据，将“两个不同元件 A+B 被连在一起”的结构嵌合拆分。它不同于亚家族聚类：前者找连接断点，后者找同一元件的进化分支。
6. Phase 1 后再以 `post_refine_merge_identity` 合并因边界收敛而变得近乎相同的结果；失败时保留未合并集合，避免工具故障造成丢失。

`REFINE_STRATEGY_DESIGN_v2.md` 描述了进一步的“representativeness-first”目标架构：对每家族先评估拷贝完整性，再**总是**对真实拷贝做距离聚类、每个合格子簇分别建共识；对缺失/分化信号家族用 sketch 预筛 + 有上限的定向 BLASTN 召回。该文档的重要生物学动机是正确的：单一共识会处在 70–95% 分化谱的中间，反而离两端都远。但在当前代码路径中，选择性召回/子家族功能受配置开关与实施状态制约，不能假定每一次默认运行都输出“每家族多个亚家族共识”。报告结果时应以 `refiner_mdl_stats.json`、shard JSON 和输出 header 为准。

### 6.4 Phase 2：最终 QC、TE 信号标记与双库输出

Phase 2 再执行长度、N%、单/二核苷酸熵、DUST 和低复杂度比例过滤。重复过滤是有意的：Phase 1 的延伸/拼接可能重新引入低复杂度端或 N-rich 片段。

TE 结构检测包括蛋白域、LTR/TIR、已知 TE 同源等。**默认是标记而不是硬删除**（`te_signal_annotation.tsv`）：远缘 TE 不应因没有可检测结构而被误删。`--enable-te-structure-gate` 才恢复旧式硬门控。

若可取得结构/同源信号，默认 MDL 原生质量门控保留满足任一条件的记录：`tier` 在 `mdl_quality_keep_tiers`（默认 `core`），或有 TE 结构信号，或 `copies >= 50`；只丢弃“低 tier 且无信号且低拷贝”的交集，并写入 `mdl_quality_dropped.tsv`。如果无法计算任何结构信号，或 header 不含 `tier=`（老版本 mdl-repeat），门控自动跳过。这一 fail-safe 设计防止只凭 `tier=core` 过度删除远缘、低信号的真实 TE。

若用户提供去除 TE 后的宿主蛋白库，强匹配细胞蛋白且没有 TE 信号的共识可被排除；终端基因样序列可被裁剪。可选非 TE 分离会把细胞器、rRNA、卫星重复写到独立 `nonte_*.fa`，不计入 TE 库。

输出有两个用途不同的库：

- `refiner_output/consensus_masking.fa`：T1 + T2，且满足长度相关的拷贝阈值；用于更保守的掩蔽/HMM 种子。
- `refiner_output/phase3_analysis_library.fa`：前者加上满足 `analysis_t3_min_copies`（默认 3）的 T3；用于最终三源合并和分类，覆盖更广。

长度相关拷贝阈值为：≥500 bp 需 2 拷贝，200–499 bp 需 3，100–199 bp 需 5，<100 bp 需 10。短序列随机匹配和简单重复风险更高，故要求更强复发证据；长序列自身的信息量/特异性更高，可用较少拷贝支持。

## 7. GFA 的两种格式：从图到可发现 FASTA 的转换

GFA 输入不会简单地“取所有非参考节点”。默认 `--gfa-mode alleles` 的目标是从**参考锚点之间的样本等位序列**生成候选，从而保留位点、样本、长度变化和质量来源。最终转换链为：

```text
GFA 图 + 路径/比对证据
  → 每个位点的候选等位序列（≥100 bp；必要时 30 kb 分窗）
  → 同一 locus 内的 exact SHA-256 去重与样本支持聚合
  → FASTA QC（N、熵、DUST、长度）
  → 每条候选独立 FASTA 记录（默认 spacer=0）或 N-spacer scaffold
  → 参考重建 + 候选 scaffold（参考无关模式）/仅 scaffold（富集模式）
```

`GfaProcessor` 用前 50 条 `S` 记录的 `SR:i:` 判断 rGFA；没有 SR 标签而有/预期有 P 路径时视为 PGGB。若图没有嵌入 `S` 序列、路径引用缺失 segment、或 PGGB `P` overlap 非零且无法安全解释，流程会停止，避免构造错误 DNA。

共同候选 QC 是：长度 100–30,000 bp、N 比例 ≤0.20、单核苷酸熵（≥100 bp 为 1.0；更短为 0.8）和 DUST ≤0.7。30 kb 上限匹配常见完整 LTR/LINE 的尺度，也限制极长重排或错误路径对发现器的干扰；超过上限的有效等位序列以默认 1,000 bp overlap 滑窗切分。末窗会后移以避免留下 <100 bp 的无效尾片。`--gfa-spacer-length 0` 是默认值：每个候选独立记录，既不引入人工 N，也避免相邻候选被发现器串接；只有兼容需要时才可打包。

### 7.1 PGGB/smoothxg：P/W 路径等位序列提取

PGGB 的细粒度 segment 通常没有 rGFA 的 `SR/SN/SO` 坐标语义，但有 `P` 路径。算法在 `gfa_alleles.py` 中实现：

1. 选择参考 P 路径。可用 `--gfa-ref` 显式指定；否则根据路径名/样本自动推断。每条参考路径构建 `ReferenceModel`：有向 `(node,orientation)` 的出现 rank、累计 bp prefix 和每种节点出现次数。
2. 出现次数超过 `max_anchor_occurrences`（默认 64）的高重复节点不作为可靠锚点；优先使用只出现一次的锚点。其原因是重复节点不能稳定划定“同一基因组位点”。
3. 对每一条非参考 P 路径，按 contig 名匹配参考路径（必要时允许唯一的大小写折叠匹配）。用 Hunt–Szymanski/LIS 风格的最大共线链算法在 alternate path 与 reference 的所有可用节点出现位置中求递增锚链；不是贪婪的“第一次匹配”，因而能在重复节点和重排中最大化共线支持。
4. 在相邻锚点之间取 alternate 内部 steps 与 reference 内部 steps，按方向拼接 segment DNA；开头/末尾的无锚部分也生成未解析边界候选。若 alt 与 ref 的步序或拼接 DNA 完全相同则跳过。
5. 以 `alt_span - ref_span` 归类：参考跨度 0 为 `insertion`；差值 ≥100 为 `insertion_like`；≤−100 为 `deletion_associated`；其余为 `replacement`。候选仍可作为发现内容，因为替换等位序列也可能携带 TE。
6. 长等位序列按窗输出。每个候选标注 path、sample、contig、两个锚点及其参考 rank、原/替换跨度、窗口坐标、locus key 和 SHA-256。
7. 只对 `(locus_key, sequence_sha256)` 相同的候选去重；相同序列在不同位点**绝不合并**。`SQLiteAlleleStore` 记录每次出现和支持样本，防止大队列把所有复制序列留在内存；最后 manifest 仍保留 occurrence 数和 sample 支持数。

若 alternate path 找不到匹配参考 contig，程序报错并列出示例，而不是静默丢弃；这避免因路径命名规则不一致造成假阴性。

### 7.2 rGFA/minigraph：GAF 支持的样本等位序列提取

rGFA 的 `S` 行包含 `SR:i`（参考 rank）、`SN:Z`（来源）和 `SO:i`（参考偏移），但 L-link 拓扑本身不说明某条非参考支路属于哪个样本或哪个完整位点。因此默认等位模式强制要求 `--gfa-gaf`；`--gfa-gaf-fasta` 可额外校验 query 名称和长度。没有 GAF 时只能显式选 `--gfa-mode legacy`（第 7.3 节），不能把它当作样本感知的等位结果。

GAF 后端（`rgfa_gaf.py`）的关键规则如下：

1. 读取 rGFA，建立有向 segments、GFA links、`(SN,SO,length)→node` 坐标索引、按 SR=0 节点重建的 reference contig/rank/prefix。
2. 流式解析 primary GAF；验证 mandatory fields、路径 walk/`cg`、节点存在性和 GFA 中相邻边。若提供 query FASTA，GAF query 长度也须一致。每条 query 的 sample/contig 由名称约定解析。
3. 只保留被 GAF path interval 完整覆盖的 node；终端部分覆盖节点不以不可靠方式“补全”。负链记录先翻转 path 并互换方向，标准化到 query-forward。
4. 用完整、正向、SR=0 的节点构建最大递增共线锚链；必须有至少两个锚点。参考样本（定义 SR=0 backbone 的 sample）不生成 alternate allele，即使 mapper 可走到一个非参考等价路径。
5. 只接受 `MAPQ ≥ 20`、`identity ≥ 0.90`（`--gfa-gaf-min-mapq`、`--gfa-gaf-min-identity` 可调）的记录。阈值的目的不是提高序列“美观度”，而是避免不确定映射制造虚假的样本存在/缺失。
6. 每对相邻锚点间拼接完整 alternate node，并与参考内部 DNA 比较；变化且长度 ≥100 bp 才成为候选。依 PGGB 同样的长度规则分窗、标记类型、写入图节点序列来源 `graph_complete_nodes`、GAF record id、query 范围、MAPQ、identity。
7. 用磁盘 SQLite 以 `(locus,sha256)` 聚合同位点精确相同等位序列。随后**第二次**流式扫描 GAF：对每个候选 locus，只有同一条高质量完整记录同时跨越两个外侧锚点才做样本 call。这样可容许另一样本在中间含有更多参考锚点，又不把断裂/不确定读段误作 absence。

输出 `gfa_rgfa_locus_presence.tsv` 中：reference sample 固定为 `0`；非参考样本有一个唯一且候选性的等位序列为 `1`，有一个唯一但不满足候选条件的等位为 `0`；没有跨越两锚点证据、或观察到多个冲突 digest 时为 `NA`。`NA` 是故意的保守未知，不可在下游按 `0` 解释。

### 7.3 历史 rGFA legacy walker

`--gfa-mode legacy` 读取 `SR:i:0` 为参考节点，从未访问的非参考节点沿 L links 正/反向走最大非分支线性链；遇参考或分支停止。链总长在 100–30,000 bp 时拼接序列并输出。无参考 FASTA 输入时，参考按 `SN` 分组、`SO` 排序后重建为 `gfa_reference.fa`。

它能为旧工作流提供候选，但不带样本路径语义，也无法可靠界定 locus、样本频率或等位关系。因此仅适合作为兼容模式，不可与 GAF 等位矩阵结果等价解释。

### 7.4 GFA 注释（可选）

当前 `anno` **只保留全基因组既有注释投影方案**。旧实现中“用最终 combined library 对拼接的 GFA path context 再跑一次 RepeatMasker”的路径已经删除；新实现不调用 RepeatMasker，也不读取 `combine/TEs.fa`。其目标是回答两个不同层级的问题：某次样本路径上的 node occurrence 覆盖了哪些 TE/CDS/exon/intron 等特征，以及汇总所有 occurrence 后，一个图节点最合理的身份是什么。

可在完整 Pan_TE 运行末尾使用 `--anno --anno-manifest annotations.tsv`，也可用 `--anno-only` 跳过 TE 库发现，直接注释已有图：

```bash
Pan_TE --gfa graph.gfa --anno-only \
  --anno-manifest annotations.tsv --out-dir graph_annotation

# 等价的独立入口
python bin/gfa_annotator.py --gfa graph.gfa \
  --manifest annotations.tsv --out-dir graph_annotation/anno
```

#### 7.4.1 manifest 契约

manifest 是制表符分隔文件，每个 sample/haplotype 一行。路径可以是绝对路径，也可相对 manifest 所在目录：

| 列 | 必需 | 含义 |
|---|---|---|
| `sample` | 是 | 与 PanSN P/W/GAF 名称中的 sample 一致 |
| `haplotype` | 否 | 默认 `0`；PGGB 常用 `1` |
| `gff3` | 二选一 | 基因结构或同时含 TE 的标准 GFF3 |
| `te_annotation` | 二选一 | RepeatMasker `.out`、BED 或 TE GFF3 |
| `annotation_complete` | 否 | `true` 表示 TE 注释在全基因组完整，没有 TE 命中才可解释为 NonTE；只有 gene GFF3 时不得据此设 true |
| `path_prefix` | 否 | 非 PanSN path/query 名称的样本前缀 |
| `seqid_map` | 否 | 两列 `annotation_seqid -> graph_seqid`，解决 NCBI accession 与 `Chr1` 等命名差异 |
| `gaf` | 否 | 该样本的 rGFA 映射；全局也可用 `--gfa-gaf/--gaf` |
| `fasta` | 否 | 输入溯源与 checksum，并为逐样本 `.tbl` 提供精确的 sequence 数、总长度、N/X-excluded 长度和 GC；投影坐标仍来自图路径/GAF而非重新比对 FASTA |

`sample#haplotype#contig` 形式的注释 seqid 会自动去掉前两段，与 P/W 解出的 contig 对齐。显式 `seqid_map` 优先于自动归一化。重复的 sample/haplotype、冲突映射、缺失文件、非法 GFF3 坐标都会立即失败，而不是静默跳过。

输入 GFF3 统一从 1-based inclusive 转成内部 0-based half-open。保留 `gene/mRNA/transcript/exon/CDS/UTR/ncRNA` 等类型；同一 transcript 的相邻 exon 间空隙会派生为 intron。TE 由 feature type/source/attributes 中的 transposon、retrotransposon、retroelement、repeat_region 或 mobile_element 信号识别。`te_annotation` 的 RepeatMasker `.out` 直接读取原基因组坐标；`Simple_repeat`、`Low_complexity`、`Satellite` 和常见 RNA 类保留为 other feature，但不计作 TE。

#### 7.4.2 PGGB/smoothxg 坐标恢复

对 P 行，优先按 PanSN `sample#haplotype#contig` 解析；也支持 manifest `path_prefix` 和单样本回退。每个 step 的节点长度来自 S sequence 或 `LN`。相邻节点起点为 `previous_end - overlap`；只接受完全由 `M/=/X` 组成的 overlap CIGAR。GFA link 是 bidirected，因此 `A+ -> B+` 与 `B- -> A-` 等价，查边时同时检查正向和反向互补键。缺边、复杂 indel overlap、节点缺失或 overlap 数量不符会使该 traversal 进入 `gfa_annotation_unresolved.tsv`，绝不近似坐标。

W 行直接使用 `seqStart`，按 walk 中 `>`/`<` 恢复方向；W 没有相邻节点重叠。对每个 occurrence 保存 genome interval、node-local interval、方向、path 和 evidence。PGGB 的 P/W 是完整装配路径，但 `Absent` 只在该 sample 的所有 manifest haplotype 都有与该节点其他 occurrence **同 contig** 的完整嵌入路径时成立；缺失该 contig 或任一 haplotype 路径则保持 `NA`。

#### 7.4.3 minigraph rGFA 与 GAF

minigraph 注释采用严格 rGFA 后端。每条 S 必须有内嵌序列以及 `SR/SN/SO`；link overlap 必须是 `0M`。程序建立 node、双向 link、`(SN,SO,end)→node` 精确坐标索引和每个来源 contig 的有序区间索引，并验证 SR=0 参考节点连续。`SN/SO` 生成来源基因组的 origin occurrence；它是该 node 来源坐标的直接证据，但不证明其他样本缺失。

真实 minigraph GAF 不只有 `>node_id`。同一文件可以混合三种 path 表示，程序均作严格解析：

1. `>node1<node2` 形式的直接 node walk；
2. `>sample#hap#contig:start-end` 形式的来源坐标 token。精确命中一个 S 区间时直接解析；跨多个 S 时只允许完整、连续、无 gap 的节点展开；
3. `sample#hap#contig` 命名来源 path，再由 GAF 的 `path_start/path_end` 截取对齐区间。

对每条 primary GAF，先验证 mandatory fields、query/path 范围、`cg` 消耗长度、解析后 path length 以及每一对相邻 GFA edge；query 名再按 PanSN、注释 seqid 或 manifest `path_prefix` 连接到对应注释。默认只投影 `MAPQ >= 20` 且 `matches/block_length >= 0.90` 的记录。

坐标变换按 `cg` CIGAR 精确执行，而不是要求 gapless：`M/=/X` 同时推进 query 和 graph path，并与每个 node 的 path interval 求交；`I/S` 只推进 query；`D/N` 只推进 graph path；`H/P` 不推进。负链从 `query_end` 向前计算并翻转 node 与 feature 方向。因此 indel 之后的注释不会整体错移，终端部分 node 也只记录实际覆盖的 node-local interval。连续的 `M/=/X` 在坐标上都是 1:1，写出前无损合并为 match run，以免 SNP 的 `=/X` 交替制造数百万个等价碎片。

一个 node traversal 仍可能被 indel 拆成多个线性 occurrence block；节点的 TE 支持率和 disagreement 按 `(样本/单倍型、contig、GAF record、node step)` traversal context 聚合，而不是按 block 数量投票。因此 CIGAR 的分块表达不会改变 `TE/PartialTE/Mixed` 身份。明细表保留 block 级精确坐标，provenance 另外记录总/有效/低质量 GAF 数、三种 path 类型、match-run 数和新增 occurrence 数。

没有高质量 GAF traversal 绝不等同于 absence。rGFA 中未观察到某样本时状态保持 `NA`；只有已有 occurrence 且该样本的 TE 注释标为完整时，才可在该 occurrence 上给出 `NonTE`。这是相对具有完整 P/W 装配路径的 PGGB 最大的能力边界。

#### 7.4.4 投影与方向变换

对 genome feature `[f0,f1)` 和 occurrence `[g0,g1)` 先取交集 `[max(f0,g0), min(f1,g1))`。正向 occurrence 的 node 坐标为 `node_start + (overlap-g0)`；反向 occurrence 为 `node_end - (overlap-g0)` 并翻转 feature strand。内部和 TSV 一律为 0-based half-open，只有输出 GFF3 转回 1-based inclusive。

真实 smoothxg 图节点很短，本实现用 1 kb feature bin 查候选，避免每个约 20 bp 节点扫描 100 kb 内所有注释。图 links 在 occurrence 构建后释放；写 occurrence 时不缓存无用行，node summary 只保留写 annotated GFA 所需的紧凑字段。即便如此，PGGB 每条样本路径、minigraph 每个 GAF match-run/node 交集都会产生 occurrence，完整图的 TSV/GFF3 可能达到数 GB，磁盘空间必须按 `路径步数 + CIGAR 线性块 + 投影碎片数` 规划。

#### 7.4.5 occurrence 与 node 身份

occurrence 是权威层。其 TE fraction 为该 occurrence 的 node-local TE 区间并集长度除以实际覆盖的 node span：默认 `>=0.80` 为 TE，`0.10–0.80` 为 PartialTE；有完整 TE 注释但无 TE 为 NonTE，不完整注释无 TE 命中为 Unknown。gene/CDS/exon 命中本身不把 TE call 变成 callable，避免把“有基因注释”错误等价为“已完整搜索且不是 TE”。

节点汇总将不同样本、不同位点的 occurrence 合并：

- `Conflict`：同一 occurrence 的 TE 与 CDS 有碱基重叠；
- `Mixed`：节点既带 TE 又带基因结构，或完整注释的 occurrence 对 TE/非 TE 判断不一致；
- `TE`：节点 TE union fraction 达阈值且至少一半 callable occurrences 支持 TE；
- `PartialTE`：有 TE 证据但未满足完整 TE 判据；
- `GeneStructure`：无 TE 且有基因结构；
- `NonTE`：有 callable occurrence、无 TE/基因结构；
- `Unknown`：没有足够可调用证据。

样本状态为 `TE/PartialTE/NonTE/Unknown/Absent/NA`，不可把 `Unknown`、`Absent` 和 `NA` 合并。`--anno-te-min-fraction` 默认 0.80，容忍边界误差和极短未注释尾部；`--anno-partial-min-fraction` 默认 0.10，用于显式标记跨节点边界的 TE 碎片，避免把任何 1 bp overlap 都提升为 TE。

#### 7.4.6 输出与审计

| 文件 | 含义 |
|---|---|
| `gfa_node_annotation_intervals.tsv` | 所有投影 feature；node、occurrence、原基因组坐标和 family/class 均保留 |
| `gfa_node_occurrence_annotations.tsv` | 每个 occurrence 的坐标、证据、callability、TE coverage 与身份 |
| `gfa_node_annotations.tsv` | 每个节点的汇总身份、支持度、gene structure、conflict、confidence 与宽格式样本状态 |
| `gfa_node_sample_matrix.tsv` | node × sample 长格式状态矩阵 |
| `gfa_annotation.tbl` | RepeatMasker 风格的唯一 GFA node-sequence 空间汇总；共享节点只计一次 |
| `gfa_annotation_by_sample.tbl` | 每个 sample/haplotype 一个区块的 RepeatMasker 风格源基因组坐标汇总 |
| `gfa_annotation_summary.tsv` | 上述两种 scope 的机器可读分类、element 数、区间并集长度、百分比、分母、N/X 与 GC |
| `gfa_node_annotation.gff3` | node-local GFF3 |
| `gfa_node_annotation_conflicts.tsv` | TE–CDS overlap 审计 |
| `gfa_annotation_unresolved.tsv` | path/GAF/seqid 无法精确解释的记录及原因 |
| `pangenome.annotated.gfa` | S 行增加 `ai` 身份、`it` 是否有 TE、`tf` TE fraction、`tc` class、`fm` family、`gs` gene structure、`cf` confidence |
| `gfa_annotation_provenance.json` | 输入 SHA-256、阈值、图模式、计数、坐标与缺失值契约 |

#### 7.4.7 RepeatMasker 风格 `.tbl` 的统计算法

`.tbl` 是对已经投影到图上的既有注释作汇总，不重新运行 RepeatMasker，也不把图拓扑分支数当成重复元件数。程序同时输出两个不可互换的坐标 scope。

`gfa_annotation.tbl` 的分母是所有 S 行节点长度之和，称为 **unique GFA node sequence space**。对每个投影区间使用 `(category,node_id,node_start,node_end)` 排序流式合并；同一 node 被多个样本、多个路径或同一路径的多个 CIGAR block 覆盖时，其相同碱基只计一次。因此 graph `bases masked`、`Total interspersed repeats` 和各 class 长度都不会因样本数而线性膨胀。若两个组装把同一共享 node 分别注释成 Gypsy 和 DNA/MuDR，Gypsy 行和 DNA 行各自都可覆盖该 node；这些 class 行是独立并集，允许相互重叠，不能相加。`Total interspersed repeats` 单独对 Retroelements、DNA transposons 和 Unclassified 三个聚合类求并集，仍是非冗余值。若所有 S 行都有实际 sequence，sequence 数、总长度、N/X-excluded 长度和 `GC/(total-N-X)` 可精确计算；只给 `*` 与 `LN` 时总长度仍可计算，但 N/X 与 GC 明确写 `NA`。

`gfa_annotation_by_sample.tbl` 在每个 `(sample,haplotype)` 的原基因组坐标中统计。先将该样本所有 occurrence 按 `(seqid,genome_start,genome_end)` 合并为 graph-callable 区间；再将每条源 annotation feature 与同 seqid callable union 求交。求交通过已经合并区间的 end 数组和 `bisect_right` 定位首个可能重叠区间，随后只线性扫描实际相邻区间，而不是对每条 feature 扫描整条染色体。交集片段按 `(sample,haplotype,category,seqid)` 的起点顺序流式求并集。manifest 有 `fasta` 时，报告分母为完整 FASTA，总长度、N/X 与 GC 精确；没有 FASTA 时，分母退化为 callable source interval union，GC/N/X 不可由坐标推断而写 `NA`。无论哪种分母，重复长度只来自 annotation 与 graph-callable 区域的交集，因而不会把 minigraph GAF 未覆盖区误说成已投影。

`bases masked` 对 TE、Simple_repeat、Low_complexity、Satellite 和常见 small-RNA 类的全部区间作联合并集；`Total interspersed repeats` 对 Retroelements、DNA transposons 与 Unclassified 三个聚合类的区间作并集，Helitron/RC 按 RepeatMasker 习惯单列 Rolling-circles，不并入该行，所以两者通常不同。分类映射保留 RepeatMasker 常见层次：SINE、Penelope、LINE 及 CRE/L2/R1/R2/RTE/L1 子类，LTR 及 BEL/Copia/Gypsy/retroviral 子类，DNA 及 hAT/Tc1/EnSpm/MULE/PiggyBac/Tourist/other 子类，不能解析的 TE 归 Unclassified。聚合行和子类行都是各自区间并集，聚合行不是对子行简单求和；某些历史 RepeatMasker alias 只进入其兼容子类行，这是为了复现 `.tbl` 的既有层次而不是重新定义生物学分类。

`number of elements` 只统计与 graph-callable 区域至少有 1 bp 交集的来源元素。RepeatMasker `.out` 按 `(class/family, RepeatMasker ID)` 合并被 insertion/deletion 拆开的 fragment，与原 `.tbl` 的 element 计数一致；GFF3/BED 等没有 fragment ID 的输入按去重后的源 annotation record 计。一个元素即使跨越许多节点仍只计一次。图级 element 数是各样本来源元素数，不是 node 数或 graph component 数。上述定义、坐标 scope、分母种类和精确数值同时写入 `gfa_annotation_summary.tsv`；该 TSV 是比较、作图和自动 QC 的权威统计接口，`.tbl` 是便于人工阅读的兼容布局。

复杂度方面，occurrence 和 feature 各进行一次排序；callable 区间先压缩，feature–callable 求交为二分定位加局部扫描；node-space 与 source-space 长度均用每个键只保存最后一个 union 尾端的流式累加器。除了已生成的投影明细及排序引用，不建立逐碱基数组。provenance 明确记录两个 `.tbl` 路径，并声明 graph 使用 unique-node union、sample 使用 callable-source union，以及 RepeatMasker fragment-ID/其他来源 record 的 element 语义。

以 `/scratch/shuoc/TE/Arabidopsis_thaliana/V2/GCA_020911765.2_ASM2091176v2_genomic.fna.tbl` 为参照，将对应 FASTA 的 5 条染色体构造成 5 个 `LN` 节点和 5 条完整 P path，并用显式 `seqid_map` 连接 `CP087126.2–CP087130.2` 与 PanSN 染色体名。程序读取 102,569 条坐标合法的 `.out` feature 并全部投影。逐样本表与原表的 sequence 数 5、总长度 133,917,131 bp、N/X-excluded 长度 133,916,631 bp、GC 36.35%、`bases masked` 43,771,690 bp 和 mask rate 32.69% 完全一致。原 `.out` 另有 18 行 `end < start`，不形成正长度区间，按输入契约排除；其中仅由非法行代表的 fragment ID 不进入 element count。Pan_TE 的 category occupied bp 使用明确的坐标区间并集，RepeatMasker 原 `.tbl` 还包含其内部 overlap/fragment 归属规则，因此分类行的 bp 可有小幅差异；自动比较 mask rate 时以 `masked` 行为准，不把分类行之和替代非冗余 masked union。

#### 7.4.8 拟南芥真实 PGGB 验证

使用 `/scratch/shuoc/TE/Arabidopsis_thaliana/upload` 下 7 套逐基因组 RepeatMasker `.out`，对 PGGB Chr1 smoothxg GFA 完整运行。图含 2,270,366 个节点、3,144,225 条 `0M` links 和 7 条 P paths；程序恢复 17,011,115 个 occurrence，读取 800,587 条全基因组 repeat features，投影 12,150,687 个区间，unresolved 为 0。节点汇总为 TE 853,457、PartialTE 10,872、Mixed 33,237、NonTE 1,372,800。

逐样本将 occurrence 表的 `te_covered_bp` 求和，再独立对原 `.out` 的 Chr1 TE 区间求并集，7 个样本逐 bp 完全相等：TAIR10 6,725,937；T2T 9,089,723；Cvi-0 10,928,548；Ey15-2 10,696,603；IP-Met-6 8,281,890；Rabacal-1 10,936,411；Tsu-0 9,235,551 bp。这同时验证了 PanSN 归一化、P path 坐标、bidirected link、反向节点变换和重叠区间去重。PGGB 将五条染色体拆为五个 community GFA 时应逐文件运行同一 manifest；node id 只在各 community 文件内部解释，不应未经命名空间化直接拼表。

#### 7.4.9 拟南芥真实 minigraph 验证

使用同一套 7 个拟南芥 RepeatMasker 结果，运行 `/scratch/shuoc/TE/Arabidopsis_thaliana/Ara.gfa` 与 `rgfa_validation_20260724/Ara_to_graph.gaf`。rGFA 含 74,660 个节点、总 node sequence 166,090,619 bp；1,376 条 primary GAF 全部通过结构、path length、坐标 token/命名 path 和 graph-edge 验证。默认质量阈值保留 521 条记录，其中 295 条为来源坐标 walk、226 条为命名来源 path，855 条低质量记录被审计性过滤。

最终构建 868,506 个 CIGAR match run、1,179,408 个精确 occurrence block，投影 1,133,307 条 TE 区间，`unresolved=0`。节点身份是 TE 24,596、PartialTE 17,621、Mixed 5,383、NonTE 27,060；所有节点都有可调用证据。node-local TE 区间并集为 66,307,931 / 166,090,619 bp（39.9227%）。这是图中唯一 node 序列的 TE 并集比例，不等同于任一线性基因组的 RepeatMasker mask rate，因为 homologous sequence 在图中共享，而一个 node 也可在不同基因组 context 中出现。

另作了不调用投影器内部数据结构的独立 bp 审计：从 occurrence TSV 合并每个 `sample#hap#contig` 的 callable genome interval；直接解析 7 个原始 `.out`，排除非 TE class 和 `end < start` 的异常记录后合并源 TE interval；再计算 `源 TE 并集 ∩ callable 并集`，与 interval TSV 中实际投影 TE 的 genome-coordinate 并集逐行比较。两者完全一致（`cmp=0`）。可映射 TE bp 为 T2T 41,719,175、TAIR10 28,571,472、Cvi-0 25,520,542、Ey15-2 26,638,897、IP-Met-6 26,671,318、Rabacal-1 26,251,204、Tsu-0 27,629,146，总计 203,001,754 bp。这里使用 callable 交集而不是全基因组 mask bp，避免把 GAF 未覆盖区域错误算作投影缺失。

连续 `M/=/X` 合并前，坐标相同但产生 9,371,468 个 block 和约 4.4 GB 输出；无损合并后 TE node-local 覆盖逐 bp不变，block 降至 1,179,408，输出约 1.2 GB。随后支持率改为 traversal-context 计数，消除了 CIGAR 碎片数对 `Mixed`/TE 支持率的影响。最终结果位于 `/tmp/pan_te_minigraph_anno_context_20260728/anno`；临时目录只用于本机验证，不作为长期发布资产。

## 8. PAF：插入候选处理算法

PAF 模式以多装配 all-vs-all 比对和所有对应 FASTA 为输入。FASTA identifier 必须在装配间唯一；程序先验证 PAF 至少 12 列、整数/链方向合法、名称和长度与 FASTA 一致，并为 `samtools faidx` 建索引。无参考模式自动选一个 base assembly，再加上 PAF 提取的 scaffold；富集模式只返回 scaffold 供追加到指定参考。

算法假设：在两个共线基因组片段间，若一侧有大缺口而另一侧无相当大小缺口，则较大一侧可能多出一个插入（常为 TE）。它不将所有 SV 都断言为 TE，而是抽取候选序列供后续三个 TE 发现器独立判断。

### 8.1 记录过滤与分组

解析为 `PAFRecord(qname,qlen,qstart,qend,strand,tname,tlen,tstart,tend,matches,block_len,mapq,cigar)`，按 `(qname,tname)` 分组。保留足够长、质量可用的比对；核心阈值为：

| 参数 | 默认 | 理由 |
|---|---:|---|
| `min_alignment_length` | 5,000 bp | 小 HSP 很难作为稳健共线锚点 |
| `min_mapq` | 20（minimap2） | 降低多重定位链造成的伪缺口 |
| `min_chain_alignments` | 3 | 至少三段锚点才有连续共线背景 |
| `min_chain_aligned` | max(50 kb, 2×30 kb)=60 kb | 两侧要有明显长于最大候选的支持，避免单个局部 SV 误判 |
| `max_chain_gap` | 500 kb | 防止跨巨大断裂/重排硬连链 |
| `overlap_tolerance` | 200 bp | 容忍比对边缘的小重叠/坐标噪声 |
| 候选长度 | 100–30,000 bp | 覆盖 MITE/SINE 到典型 LTR/LINE，排除微小噪声和巨型重排 |
| `gap_ratio_threshold` | 0.30 | 两侧缺口至多约 1:3 才称显著不对称 |

### 8.2 共线链

同一 query-target 对内，按 target start 排序，贪婪延长链。每一步 target 必须单调前进，target/query gap 都不得超过 500 kb，也不得反向重叠超过 200 bp。负链时 query 随 target 前进应在原坐标中递减，所以 gap 用 `prev.qstart - curr.qend` 计算。链的最坏复杂度是 O(n²)，但按 target 坐标排序后，一旦 gap 超过上限便可提前终止；该复杂度适合通常每对数百到数千条比对的场景。

### 8.3 两类候选提取

**A. 链间不对称缺口。** 对链中相邻 `curr,nxt` 计算：

```text
t_gap = nxt.tstart - curr.tend
q_gap = nxt.qstart - curr.qend       (+ 链)
q_gap = curr.qstart - nxt.qend       (- 链)
```

令 `max_gap=max(|q_gap|,|t_gap|)`、`min_gap=min(...)`。只有 `100 ≤ max_gap ≤ 30,000` 且 `min_gap/max_gap ≤ 0.3` 时保留。较大 gap 的一侧是插入候选；要求坐标为正增长并排除参考序列名称。这个判断的价值是利用两侧基因组相对信息，较少把两边共有的 assembly gap 当成 TE。

**B. minimap2 CIGAR 内部 indel。** 当 PAF 被识别为 minimap2 且带 `cg:Z` 时，逐操作走 CIGAR。`M,=,X` 同时推进 query/target；`I` 在 query 上提出候选，`D` 在 target 上提出候选；`N` 仅推进 target，`P` 不推进。对于负链，query 光标从 `qend` 向 `qstart` 递减。程序重算 CIGAR 消耗长度；有 `S/H`、重构字符串不一致或消耗不等于 PAF 坐标跨度的记录直接跳过。此验证很重要：一条损坏 CIGAR 不应把后续所有插入坐标平移。

### 8.4 合并、提取、QC 与打包

同一序列上的候选按坐标排序；重叠（严格 `start < current_end`）且合并后不超过 30 kb 的区间合并，source 记为 `gap`、`cigar` 或 `cigar+gap`。随后坐标 clamp 到 FASTA 长度并重新过滤。每 500 个区间一次调用 `samtools faidx` 提取；即使命令部分失败也按返回 header 精确匹配成功条目，避免一个坏区间让整批丢失。

提取 DNA 后按 100–30,000 bp、N≤20%、单核苷酸熵、DUST≤0.7 过滤；以 `cd-hit-est -c 0.8 -aS 0.8 -aL 0.8 -G 0 -r 1 -mask NX` 做跨单倍型候选去冗余。80% 是候选级别而不是最终家族级别的阈值：它尽早压缩显然同源的插入，却把真正家族归并留给下游共识/最终合并。

通过 QC 的序列以 5,000 个 `N` 间隔拼成最长 50 Mb 的 `paf_scaffold_N`。N spacer 让发现器不会把相邻候选当作一个连续元件，同时减少超多短 FASTA record 的 I/O。中间文件为 `paf_raw_insertions.fa`、`paf_filtered_insertions.fa`、`paf_te_insertions.fa`，统计写入 `paf_processing_stats.json`；主流程完成时通常删除前三者之外的临时原始/过滤文件。

## 9. TE-looker：原理、当前实现与门控

### 9.1 目标架构（原理层面）

项目设计把 TE-looker 定位在约 45–75% 共识分化的“twilight zone”。目标 `dtr run` 架构为：

1. 低复杂度/串联重复/NUMT 掩蔽；
2. 三轨证据：已知/mdl 衍生 HMM 的 nhmmer，sourmash sketch 证据，FastGA 自比对；
3. 将基因组窗口建成 k-NN 图；
4. Leiden CPM（可双 γ）社区发现；
5. 每社区经 abPOA/MAFFT 与边界修正建立共识；
6. 家族级合并（sketch 预筛、blastn、系统发育门）；
7. LTR/TIR/Helitron 等结构特征和 G1–G6 硬门，输出 `families.fasta` 与规范化 provenance。

其思想是不要只依赖单一高 identity 比对：不同证据轨在弱信号区域相互补充，窗口图把分散的局部相似性转化为社区，社区共识再抵消单拷贝突变。mdl-refine 的任务是补齐“已有 mdl 种子家族”的较远拷贝；TE-looker 的任务是发现根本没有 mdl 种子的更深远家族，两者可重叠但不相同。

### 9.2 当前 Pan_TE 实际调用

当前 `Pan_TE` 对处理后的**未整体硬掩蔽**基因组调用：

```text
dtr run --genome genome/genome.fa --out te-looker/run --threads N
```

当所安装的 `dtr` 已包含本项目新增的内部完整阶段时，推荐在 Pan_TE 中使用正式开关：

```text
Pan_TE ... --te-looker-final-gated
```

它会先执行 `dtr run --help`，确认二进制明确声明 `--final-gated`；不支持时会在 discovery 前失败，而不会把一个同版本但能力不同的 dtr 误当作完整实现。为兼容旧调用，`--te-looker-extra-args "--final-gated"` 也会触发同一能力检查。该模式依次运行 discovery、refine、成员共线片段拼接、copy/length gate 与跨 seed CD-HIT 合并；原始精修库保存为 `families.discovery.fasta`，拼接中间库为 `families.stitched.fasta`，门控中间库为 `families.gated.fasta`，最终库提升为 `families.fasta`。它在根目录写入明确的 `te-looker-gated-v1` provenance，因此 Pan_TE 可直接将该 `families.fasta` 作为轨道最终库。`--final-gated` 是显式 opt-in：未具备该能力的旧 dtr 继续使用下述宿主侧 fallback，不能根据版本号推断该 flag 一定存在。

如果 `mdl-repeat/refiner_output/consensus_masking.fa` 存在，Pan_TE 为每条至少 30 bp 的共识生成单序列 Stockholm block，使用 `hmmbuild --dna` 生成 `mdl_repeat_seed.hmm`，并传入 `--dfam-hmm ... --enable-known-hmm-discovery`。当前 `dtr-hmm` 已真实执行 `nhmmer --tblout`：把 1-based inclusive hit 坐标转为 0-based half-open BED，同一模型/染色体的重叠域（包括双链同位点）只保留最高分命中，再要求每模型至少 5 个、每命中至少 80 bp，交给 `te-refine` 建 POA 共识。HMM 共识与 de novo 共识随后共同进入 stitch/gate/merge；选择 masking library（T1+T2）而非宽松 T3 的理由是 HMM 种子需要保守、边界较可靠的训练序列。Track 1 扩展已知家族召回，但不阻止其他轨发现新家族。

`dtr` 的能力按实际二进制探测，不能只相信命令行接受某个 flag。当前实现特别区分：

- `--w`：真正控制 `(w,k)` minimizer 采样的种子密度。默认 1，等同每位置精确发射；`w>1` 约保留 `2/(w+1)`，显著减内存但改变结果。
- `--min-count`：canonical 16-mer 的出现次数下限。默认会随基因组大小设为：<300 Mb 50；300 Mb–1 Gb 100；1–3 Gb 200；3–10 Gb 300；>10 Gb 500。大基因组背景 k-mer 碰撞/低拷贝噪声更多，`dtr` 内置 20 过宽；若 `w>1`，该值按密度重新缩放，但不低于 5。
- `--dfam-hmm` + `--enable-known-hmm-discovery`：同时给出才启用真实 Track 1；因此预检把 `nhmmer` 与 `spoa` 列为 TE-looker 的必要依赖。HMM 命中只是候选成员，不绕过最终的 copy/length/CD-HIT 提升边界。
- `--window-stride`：只在探测到支持它的未来版本（≥0.4.0）才转发。它是窗口图顶点间距，不是 `--w` 的别名；旧版即使接受该参数也可能忽略，因此不得把它误写成已生效的降采样。
- `RUST_MIN_STACK=64 MiB`：为 rayon 工作线程提供足够栈，降低大基因组并行收集中的 stack overflow 风险。

### 9.3 为什么 Pan_TE 还要做宿主侧后处理

`dtr` 的基础 `run` 仍可只做每个重复 16-mer seed 的粗共识（discovery + spoa refine）；若将这种输出直接合并，单个元件可能产生大量 offset/fragment family，导致数万条冗余共识主导库。`--final-gated` 已将 Pan_TE 需要的 fragment stitch、copy/length gate 和跨 seed merge 内置为可审计阶段，但 Pan_TE 保持防御性验证：仅当 `te-looker/run/provenance.json` 是有效 JSON，且显式声明 `contract=te-looker-gated-v1`、`mode=gated`、`final_library=families.fasta` 时，才原样信任 `families.fasta`；仅有同名文件、损坏 JSON 或旧格式均走宿主后处理。否则执行：

1. 若有 `disc.members.bed`，以成员实例的共线链检测同一元件的碎片家族；跨至少两个家族的重叠/近邻拷贝重新提取、MSA、建全长 stitched consensus。拼接失败是 non-fatal，保留原输出而非虚构序列。
2. 从 header 读 `members=N`，丢弃 `N<5` 或序列 <100 bp 的家族。没有 `members=` 时按最低 5 的保守默认解释。
3. `cd-hit-est -c 0.80 -aS 0.6 -n 5 -r 1 -g 1` 合并跨 seed 的部分重叠片段。`aS=0.6` 允许短片段与较长共识有足够覆盖而被合并；这比最终跨轨 0.9 有意更松。
4. 若最终家族数超过 `15 × genome_Mb`，硬失败并要求检查门控，而不是把膨胀库交给 combine。

因此使用结果时应区分：`run/families.fasta` 是 dtr 原始输出；`te-looker/consensus.fa` 是 Pan_TE 送入 combine 的最终轨道库；只有前者的 `provenance.json` 满足上述 gated 契约时二者才可等价，否则后者已被 stitch/gate/merge。每次运行还会写入 `run/pan_te_handoff.json`，记录输入/中间/输出 FASTA 的 SHA-256、家族数、总碱基数、门槛、CD-HIT 参数及拼接状态。 

### 9.3.1 嵌套实例分解与自适应边界重估（已实现，当前为显式后处理）

`dtr-nested`、`dtr-boundary` 与 `dtr-boundary-gate` 为已门控的 `families.fasta` 提供保守、可回退的额外 refinement 路径；目前不会由 `dtr run --final-gated` 或 Pan_TE 自动调用，必须显式运行并审阅审计文件后才替换库。这一边界避免把仍属 partial 的 Track 2 结果误混入最终库。

流程为：

1. `dtr-copy` 为原始 final library 建立全基因组 copy catalog（BED 为 0-based、half-open，保留方向）。
2. `dtr-nested --min-flank 80` 在同一染色体上寻找 inner 被不同 family 的 outer 完整包含、且 outer 两侧各至少保留 80 bp 的关系。它输出 `nested.events.tsv`，并将含内层元件的外层实例从 `clean.copies.bed` 排除。这里的分解是实例拓扑分解：内层家族仍由自己的实例代表，外层边界只从未嵌套实例取证；它不臆造或切割新的 consensus 序列。
3. `dtr-boundary` 对每个 family 的 clean instances 提取两侧 flank（默认 100 bp），按方向统一后以 POA 建候选共识。只有实例数达到 `--min-members`（默认 5）、候选长度处于原长度的 `1/max-length-ratio` 到该倍数之间（默认 ratio=1.5）时才提出替换；否则逐字保留原序列。`--max-poa-bases`（默认 4,000）和 `--max-refinement-families` 用于限制 POA 的二次复杂度，超出预算的 family 不是被丢弃，而是保留原边界，可在另一批次重估。
4. 对 proposal library 再运行 `dtr-copy`，`dtr-boundary-gate` 按 family/染色体对实例区间求并集；仅 proposal coverage 相对原库至少增加 `--min-relative-gain`（默认 0.01，即 1%）且序列实际改变时，才接受 proposal。其他 family 原样回退。这个 coverage proxy 只用于选择，不取代最终 RepeatMasker 验证。

拟南芥 `/scratch/shuoc/TE/Arabidopsis_thaliana/demo/genome/genome.fa` 的 69-family final library 测试中，原库有 2,307 个 copy instances；检测到 176 个跨 family 嵌套事件，保留 2,242 个 clean instances。受测试资源预算限制，本次只对最多 5 个、长度 ≤2,000 bp、每家族最多 5 个实例的 family 进行 POA；其中 2 个提案通过 coverage gate：`stitched_3_ref` 从 1,869 bp 到 2,099 bp（coverage proxy +13.18%），`stitched_21_ref` 从 1,196 bp 到 1,472 bp（+372.33%）。其余 67 条记录未改动。以与基线相同的 RepeatMasker 4.2.1、默认模式、`-pa 2` 重新掩蔽，masked bases 从 23,012,842（17.184454%）增加到 23,018,188（17.188446%）：净增 5,346 bp、+0.003992 个百分点（相对 +0.023231%）。完整证据位于 `validation/arabidopsis_dtr_final/boundary_refine/` 和 `validation/arabidopsis_dtr_final/mask_rate_boundary/`。

### 9.4 Track 2 的 H0 spaced-seed 成本门禁

在实现窗口图前，`dtr-spaced-seed` 对 `/scratch/shuoc/TE/Arabidopsis_thaliana/demo/genome/genome.fa` 进行了精确 canonical spaced-seed 计数；默认 pattern 为 `11101110111011101111`（span 20、weight 16，满足 v4 的 weight≥16 特异性下限），仅将出现次数 ≥3 的 seed 计为重复候选。结果为 96,732,256 个 distinct seed 中 3,754,268 个重复，合计 32,220,331 次发生，最高计数 27,943。即使每 seed 以 hash-bottom cap=200 限制，桶内配对上界仍为 346,260,134；不封顶为 56,151,837,481。故 Track 2 不能实行候选组全配对，必须采用 LSH 桶上限、每窗口固定邻居预算与精确复核。审计 JSON 位于 `validation/arabidopsis_track2/spaced_seed_h0.json`；该阶段只测量成本，不产生 TE 家族。

### 9.4.1 Track 2 的有界 LSH 候选图（当前部分实现）

`dtr-lsh` 接收候选窗口 FASTA，构建 canonical k-mer bottom-k MinHash（默认 k=15、128 hashes），按 LSH band 提议候选对；每个 bucket 以确定性 hash 截断到 128 个节点，再对候选做精确 Jaccard 复核（默认 ≥0.08），并以每节点最多 16 条边限制最终图。短候选若 distinct k-mer 少于完整 sketch，不会被丢弃，而是使用其可用的完整 band。这一阶段只输出边图和 `te-looker-lsh-graph-v1` provenance，**不将连通分量或社区直接提升为 TE 家族**，因此不会绕过 POA、结构证据或最终 gate。拟南芥的 5,000 条 sampled windows 在该设置下产生 428 条边、0 个超大桶；审计位于 `validation/arabidopsis_track2/window_lsh_provenance.json`。

### 9.4.2 候选窗口、确定性社区与组件 POA（当前部分实现）

`dtr-window` 从重复 spaced-seed occurrence 生成窗口。窗口由全局 hash-bottom 保留，而不是每个 seed 各取固定数，因此高频 seed 不会支配样本；窗口 header 使用 `(record_index, seed_start, seed_code)`，选择 hash 只用于排序，避免 hash 碰撞覆盖不同窗口。拟南芥的 32,220,331 个合格 occurrence 被全局预算限制为 5,000 个唯一窗口（`spaced_windows_5000.fa`）。

`dtr-community` 读取所有窗口节点及 LSH 边图，输出稳定排序的连通分量 membership；孤立节点被保留，分量只是后续细分容器，`family_call=false` 禁止直接作为 TE 家族。此试验得到 4,681 个分量、214 个多节点分量、最大 36 个节点（`window_components_provenance.json`）。

`dtr-community-refine` 在这些容器内运行确定性的 first-level Louvain local moving（resolution=1，固定节点排序和 tie-break），以模块度局部移动切断桥接；它**不是完整 Leiden**，provenance 明示 `deterministic_first_level_louvain_local_moving_not_leiden`。真实图由 4,681 个容器细分为 4,689 个社区，最大社区从 36 降至 13；满足 5–30 成员 POA 门槛的社区仍为 7 个。`dtr-component-poa` 现在消费细分后的 membership。它以共识消减单拷贝突变，虽然 9.3.1 的可选 final-library 后处理已实现边界重估，Track 2 本身尚未自动接入该阶段，也没有完成 Leiden aggregation、分类或最终 gate，故结果绝不进入 `families.fasta`。

### 9.4.3 全基因组 copy catalog（当前部分实现）

`dtr-copy` 把 provisional consensus 对全基因组运行 `blastn -task blastn -dust no -soft_masking false`，以 query/染色体/链分组；按 query 相邻间隙 ≤200 bp、genome 相邻间隙 ≤500 bp 合并 HSP，以 HSP 长度加权 identity，要求合并实例 query coverage ≥0.60、identity ≥75%。最后，同一 query/染色体/链上相距 ≤100 bp 的重复起点只保留一个实例。输出 BED 为 0-based half-open（第 4 列为 consensus ID，第 6 列为链），同时保留原始 HSP 和 `te-looker-copy-catalog-v1` provenance。

在 7 条拟南芥 provisional consensus 上，171,454 个 raw HSP 合并并筛为 253 个实例；与先前独立手工实现的 253 个实例在坐标归一化后完全一致（253/253）。文件为 `validation/arabidopsis_track2/copy_catalog/dtr_copy_instances.bed` 和 `dtr_copy_provenance.json`。copy catalog 是候选实例证据，不能替代结构检验或家族提升门控。

### 9.4.4 Stage 8 non-TE 守卫（rRNA 分支已验证）

`dtr-nonte-guard` 需要显式提供 non-TE reference（例如 rRNA、细胞器 FASTA）；它对 provisional consensus 运行 BLASTN，以 identity ≥80% 的 HSP 在 query 坐标上作区间并集，若高 identity query coverage ≥0.80，则将候选从 interspersed TE 候选中拒绝。该策略把“命中一个短保守片段”和“整体可由非 TE 序列解释”区分开；输出 retained/rejected FASTA、候选级 TSV 原因、原始 HSP 与 `te-looker-non-te-guard-v1` provenance。即使被拒序列可留作 rDNA/ncRNA-array 注释，它也不能进入 TE 家族库。

以已有 `nonte_rrna.fa` 对上述 7 条候选验证，6 条被拒（coverage 0.999593–1.0、identity 87.302–100%），只有 `track2_cc_001577_partial` 保留。审计文件位于 `validation/arabidopsis_track2/non_te_guard/`，但这只是 rRNA-only 历史中间结果，不能视为完整 E1。

随后下载官方 NCBI Arabidopsis chloroplast `NC_000932.1` 与 mitochondrion `NC_037304.1` FASTA，并与 supplied `nonte_rrna.fa` 合成 5-record panel；organelle FASTA SHA-256 为 `a4dc6bd89e88975f5be6b7c69224f234b64ca14f78bfb39ff0e72fbd5e823470`，完整 panel 为 `d45b5e19905f2343151ebe76803cba97d8f8467cb44618895fd6674e08e78875`。在既有拟南芥 combine library 的 3,787 条序列中，守卫拒绝 30 条：21 条最佳命中 `NC_037304.1`，9 条为 rRNA；被拒对象中包括原标为 LTR/Gypsy、LINE/L1 与 tRNA 的序列。这直接验证 E1 所针对的污染，而不是仅验证命令可运行。

更新后的 Track 2 运行使用 organelle+rRNA panel；由于社区细分改变了 POA 成员集，7 条 provisional consensus 全部由 rRNA reference 解释（0 retained）。因此 `validation/arabidopsis_track2/orchestrated_panel_refined/` 的真实结论是 eligible=0，不存在可提升候选。

随后从官方 `GCF_000001735.4_TAIR10.1` GFF 中提取 684 个 `tRNA` 与 82 个 `snRNA` feature：将 GFF 的 1-based inclusive 坐标转换为 BED 0-based half-open，使用同一官方 FASTA、`bedtools getfasta -s` 按链提取。v3 panel 含 771 records（organelle+rRNA+tRNA+snRNA），SHA-256=`cf852bb11a71d1e40dd655132b31107b878e5038beb57e7c6777fa98b6f4e8af`。在同一 combined library 上，v3 拒绝 34/3,787 条：11 条最佳为线粒体、9 条 rRNA、12 条 tRNA、2 条 snRNA；后两类中包含原标为 LTR、LINE、hAT、Copia 和 Unknown 的序列。审计在 `validation/arabidopsis_track2/non_te_panel_v3_curated/`。

v3 manifest 完整覆盖 `organelle/rrna/trna/snrna`，因此可令 **G6 的 non-TE 子门**为 true；在旧的单候选正例上它确实变为 true，但 `acceptance_complete=false` 仍保持，因为蛋白/宿主基因守卫尚未实现。故 provenance 仍为 `mode=partial, family_call=false`，禁止输出最终家族。

### 9.4.5 实例侧翼与末端结构审计（当前部分实现）

`dtr-structure` 读取 eligible consensus、0-based half-open copy BED 和基因组；对每个实例比较两侧 flanks，枚举 k=3–15 的 TSD（禁止 k<3，因为 mariner 的 TA 没有区分力），报告支持实例数、可用侧翼数和比例。另在共识两端最多 500 bp 中，报告 ≥12 bp 的 direct terminal seed、其 reverse-complement terminal seed，以及末端连续 poly-A/T 长度；实例 <20 时标注 `boundary_uncertain=true`。结构证据的条件是：≥5 个实例且 TSD 支持比例 ≥0.20，或 ≥12 bp 的 terminal seed/poly-A/T。它的 `te-looker-structure-audit-v1` provenance 永远是 `mode=partial, family_call=false`：这些特征可能来自内部重复或不精确边界，故只能帮助后续 G1–G6 审计，不能独立分类或提升家族。

对 rRNA 守卫留下的 1 条候选，86/86 个实例具可用侧翼；TSD 仅 2/86（3 bp），没有 TSD 支持；但其共识报告 80 bp direct、36 bp inverted terminal seed。结果位于 `validation/arabidopsis_track2/orchestrated_rrna_guard_with_structure/structure_audit.tsv`，是结构线索而非最终 TE 断言。

### 9.4.6 G1–G6 接纳审计骨架（当前部分实现）

`dtr-accept-audit` 不接受“一个总分”来掩盖缺失证据，而是输出每条候选的门级列：G1 为默认 ≥10 个全基因组实例的 copy support；G2 为结构审计线索；G3 为至少 20 个实例时才不标 `boundary_uncertain`；G4 为显式 non-TE report 的 `retain`；G5 标记“少于 20 拷贝且平均 identity ≥98%”的疑似 segmental duplication；G6 的 non-TE 子门由 `--non-te-manifest` 机器判定：manifest 必须覆盖 `organelle`、`rrna`、`trna` 与 `snrna` 四类；列名为 `g6_non_te_panel_complete`，避免把该子门误说成全部外部守卫。无 manifest 或缺少任何类别时该列为 false；即使该列为 true，蛋白/宿主基因守卫未完成时 provenance 仍固定 `acceptance_complete=false`，因此最强状态只能是 `eligible_for_full_guard`，从机制上阻止当前 Track 2 提升为最终 family。

统一拟南芥运行的剩余候选为 86 copies、mean identity 91.199%、G1–G5 分别为 true/true/true/retain/false，状态为 `eligible_for_full_guard`，但 G6=false。记录在 `validation/arabidopsis_track2/orchestrated_full_partial_audit/acceptance_audit.tsv` 与 `acceptance_provenance.json`；它同时验证了 gate 的积极路径和“不允许完成接纳”的安全边界。

### 9.4.6 可复现的 Track 2 partial 编排

`dtr-track2 --genome G --out O [--non-te-reference R --non-te-manifest M --max-windows N --threads T]` 顺序调用 window、LSH、community、component POA、可选 non-TE guard、copy catalog、structure audit 和 acceptance audit。若提供 R，guard 在编目之前运行；若同时提供 M，acceptance audit 将按 M 的类别完整性自动计算 G6 non-TE 子门，避免将已解释的 rRNA/细胞器候选计为 TE 拷贝；若未提供 R，流程仍能运行但总 provenance 保持 partial。它只产生 `provisional.fa`、`non_te_retained.fa`、`copy_instances.bed`、`structure_audit.tsv`、`acceptance_audit.tsv` 及逐阶段 JSON，不产生 `families.fasta`。

拟南芥的统一运行采用 `--max-windows 5000 --threads 4 --non-te-reference .../nonte_rrna.fa`，得到 5,000 windows、7 条 provisional、1 条 eligible consensus、86 个该 eligible consensus 的实例；总审计为 `validation/arabidopsis_track2/orchestrated_rrna_guard_with_structure/track2_provenance.json`，其中强制写入 `contract=te-looker-track2-v1`、`mode=partial`、`family_call=false`、`acceptance_complete=false`、`final_library=null`。这实现了当前各阶段的可重复闭环，但绝不是“完整 TE-looker 已完成”的声明。

### 9.4 已完成的拟南芥端到端验证

在 `/scratch/shuoc/TE/Arabidopsis_thaliana/demo/genome/genome.fa` 上，以 `--threads 4 --min-count 50 --w 1 --final-gated` 完成测试。发现/精修得到 337 个家族（673,452 bp）；共线成员链检测到 37 条可信链，重建 37 条 stitched consensus，并替代 85 个源碎片家族，拼接库为 289 个家族（673,719 bp）。随后 min-members=5、min-length=100 bp 门控保留 211 个家族（326,477 bp），`cd-hit-est -c 0.80 -aS 0.6 -n 5 -r 1 -g 1` 生成最终 69 个家族（167,509 bp）。最终 `provenance.json` 的所有完成阶段和 `families.fasta` 已由 Pan_TE 的契约检查接受。验证副本位于本仓库 `validation/arabidopsis_dtr_final/`；运行细节应以其中 JSON 审计记录为准。

在同一 69-family final library 上，9.3.1 的嵌套/边界后处理进一步产生一个 2-family 替换的 `families.boundary-gated.fasta`；它的独立 RepeatMasker mask rate 为 17.188446%，相对未边界重估的 17.184454% 增加 5,346 bp（+0.003992 pp）。这是小幅、正向的实测收益，不应外推为所有物种或所有家族均会提升；应始终同时保留原库、proposal、gate audit 及 matched RepeatMasker `.tbl`。

## 10. 三源合并、分类和最终输出

合并输入按优先文件路径为：

- LTR：`Look4LTRs/consensus_lib/ltr_library.fa`；
- mdl：`mdl-repeat/refiner_output/phase3_analysis_library.fa`；
- TE-looker：`te-looker/consensus.fa`（兼容旧名 `consensi.fa`）。

先直接串接为 `combine/raw_te_library.all.fa`，再运行：

```text
cd-hit-est -c 0.9 -aS 0.8 -aL 0.8 -G 0 -n 8 -r 1 -mask NX
```

输出为 `combine/raw_te_library.fa`，并保留兼容文件名 `raw_TEs_combined.fa`、`raw_TEs.fa`。90% 比轨道内部的 80% 更严格，原因是每个轨道已经自行合并；此处若仍以 80% 跨轨合并，可能错误吞并相关但独立的 TE 家族。`-r 1` 同时考虑反向互补，`-aS/-aL 0.8` 防止只凭很短局部命中吞并一条完整家族。

提供 `--model-dir` 时，`run_Classifier` 在隔离的 `ClassifyTE_env` 中依次/组合 RepeatClassifier 与 ClassifyTE，再以 `Combine_for_Two` 汇总，产出分类库 `combine/TEs.fa`。无模型时，`raw_te_library.fa` 是未分类的去冗余从头库，不应把它误称为已分类 TE 库。

## 11. 关键参数速查与调参原则

| 参数 | 默认/行为 | 何时调整 | 风险 |
|---|---|---|---|
| `--threads` | 4；建议 4 的倍数 | 依 CPU、I/O、MAFFT/rmblastn 并发调整 | 过高会让多个外部工具争抢内存 |
| `-M` | 0（不额外限制） | 与调度器资源一致 | 过小会在预检阻止 dtr；不可指望自动降采样 |
| `--fragment-size` | 40,000 | 兼容分片相关工具时 | 不是 GFA 长等位窗口参数 |
| `--te-looker-w` | 1 | 内存确实不足且接受结果变化时 | seed 集改变，跨运行不再完全可比 |
| `--te-looker-extra-args` | 原样传递 | 传递非正式 dtr 参数 | 仅将其作为高级兼容入口；若含 `--final-gated` 仍会做能力检查 |
| `--te-looker-final-gated` | 显式 opt-in | 需要 dtr 内部 stitch/gate/merge 审计时 | 只接受在 `dtr run --help` 中声明该能力的二进制；旧 dtr 应不加该开关并由宿主 fallback 处理 |
| `--gfa-gaf-min-mapq` | 20 | 比对质量已知高/低时 | 降低会放大错误 presence 和错误等位 |
| `--gfa-gaf-min-identity` | 0.90 | 更远样本可能需要降低 | 低于此值要特别审计错误路径/锚链 |
| `--gfa-window-overlap` | 1,000 bp | 很长 TE 或重复边界需更大重叠时 | 过大增加重复候选和计算 |
| `--gfa-spacer-length` | 0 | 只在下游不能处理多 record 时设定 | 人工 N 可能影响某些工具的边界行为 |
| `--gfa-mode` | `alleles` | 仅旧 rGFA 无 GAF 时用 `legacy` | legacy 不具样本/位点解释能力 |
| `--enable-copy-recruit` | 关闭 | 需要用真实基因组复制数作严格门控 | 增加 rmblastn 计算；阈值过严丢低拷贝 TE |
| `--hard-min-copies` | 0（关闭） | 目标明确为高置信高拷贝库时 | 把稀有/近期插入直接丢掉 |

调参总原则是先改变能说明其统计/生物学含义的显式阈值，并把原始统计、manifest 和版本一同保存；不要因为内存或输出数量不理想而悄然改变 `w`、GAF 阈值或 copy floor。前者改变发现空间，后者改变结果的可比性。

## 12. 交付物、溯源与建议 QC

主要文件如下（小写目录为当前主路径，旧大写名仅为兼容）：

| 文件 | 内容与用途 |
|---|---|
| `genome/genome.fa` | 清洗、索引后的发现输入 |
| `genome/paf_enrichment/paf_processing_stats.json` | PAF 各过滤阶段计数与比对类型 |
| `genome/gfa_enrichment/gfa_allele_manifest.tsv` | 每个 GFA 候选 occurrence 的位点、样本、hash、QC、GAF 溯源 |
| `genome/gfa_enrichment/gfa_processing_stats.json` | GFA 格式、候选、去重、QC、spacer 等统计 |
| `genome/gfa_enrichment/gfa_rgfa_locus_presence.tsv` | rGFA+GAF 的 0/1/NA locus × sample 矩阵 |
| `Look4LTRs/consensus_lib/ltr_library.fa` | LTR 轨道共识 |
| `mdl-repeat/tmp/mdl_repeat.raw.fa` | 原始 MDL 家族 |
| `mdl-repeat/refiner_output/refiner_mdl_stats.json` | refine 各阶段统计/门控审计 |
| `mdl-repeat/refiner_output/te_signal_annotation.tsv` | TE 结构信号标记（默认不删除） |
| `te-looker/run/families.fasta` | dtr 原始输出；仅在 `provenance.json` 满足 `te-looker-gated-v1` 契约时才是最终门控库 |
| `te-looker/run/pan_te_handoff.json` | Pan_TE 交接审计：模式、输入/输出 SHA-256、家族数、bp、门控与合并参数、拼接状态 |
| `te-looker/consensus.fa` | 进入最终 combine 的 TE-looker 库 |
| `combine/raw_te_library.all.fa` | 合并前所有轨道共识 |
| `combine/raw_te_library.fa` | 跨轨 0.9 去冗余库 |
| `combine/TEs.fa` | 有 `--model-dir` 时的分类最终库 |
| `anno/gfa_node_annotation_intervals.tsv` | GFA node-local 的 TE/基因结构投影区间 |
| `anno/gfa_node_occurrence_annotations.tsv` | occurrence 级坐标、callability 与身份（权威层） |
| `anno/gfa_node_annotations.tsv` | 节点汇总身份与样本状态 |
| `anno/gfa_annotation.tbl` | 唯一 GFA node sequence 空间的 RepeatMasker 风格汇总 |
| `anno/gfa_annotation_by_sample.tbl` | 各 sample/haplotype 的源坐标 RepeatMasker 风格汇总 |
| `anno/gfa_annotation_summary.tsv` | 两种 scope 的机器可读长度、数量、百分比与分母 |
| `anno/pangenome.annotated.gfa` | 带节点身份、TE fraction、class、gene structure 和 confidence 标签的 GFA |
| `anno/gfa_annotation_provenance.json` | GFA 注释输入 checksum、阈值与缺失值契约 |

建议最少进行以下 QC：

1. 检查每阶段检查点是否对应实际非空输出，尤其是从旧运行恢复时。
2. 对 PAF 检查 `gap_candidates`、`cigar_candidates`、`merged_candidates`、`qc_passed` 的级联比例；若候选骤增，优先检查 PAF aligner/CIGAR 和 reference 名称映射。
3. 对 GFA 审阅 manifest 中 `locus_key`、`support_sample_count`、`qc_status`、`sequence_source`；rGFA 矩阵中把 `NA` 与 `0` 分开统计。
4. 对 mdl 检查原始数、Phase 0 保留数、copy gate、Phase 1 输出、Phase 2 MDL quality gate 的数目及 `mdl_quality_dropped.tsv`；若结构资产缺失，确认门控是否按设计跳过。
5. 对 TE-looker 验证 `provenance.json` 是否满足 gated 契约，检查 `pan_te_handoff.json` 的输入/输出 SHA-256、原始与最终家族数/bp、stitch/gate/merge 审计及每 Mb 15 家族的 sanity bound；Track 2 的 partial 输出还须检查 window/LSH/component/POA/copy/non-TE provenance，且不得被误作为最终库。
6. 对最终库抽样用 RepeatMasker 或独立参考注释验证；GFA 注释则逐样本比较 occurrence `te_covered_bp` 总和与源注释区间并集，并要求 unresolved、Conflict 和 Unknown 比例有可解释原因。

## 13. 局限性与解释边界

- 从头发现的“重复”不自动等于 TE。Pan_TE 的结构标记、非 TE 分离、蛋白负过滤和可选分类用于降低该差距，但最终生物学解释仍应结合分类证据和人工审阅。
- PAF 不对称缺口代表候选插入，而不是 TE 的直接证明；复杂重排、assembly gap、重复比对仍可产生候选。
- PGGB 等位提取依赖路径命名/contig 匹配和可靠锚点；高度重复区域可能只能以不稳定锚链或未解析候选表示。
- rGFA 候选的 DNA 来源为完整图节点，不是对 query FASTA 作 CIGAR 重放；末端部分节点被保守排除，因此结果偏向高置信但可能少于理论可恢复的序列。
- GFA 注释的 PGGB absence 依赖完整 P/W 路径；rGFA 若无每个样本的高质量 GAF 只能报告来源节点，不能推断其他样本缺失。含 indel 的 GAF 已作精确 CIGAR block 投影，但 query-only insertion 没有对应 graph node，本身不会凭空生成可注释节点。
- >2 Gb 的 mdl 路径在发现/精修阶段使用样本；Phase 2b 可恢复全基因组 copy-count 证据，却不能凭空发现样本完全遗漏的家族。
- TE-looker 的完整多轨社区/门控架构应以所安装 `dtr` 的显式 `te-looker-gated-v1` provenance 契约和版本能力为证。当前 Pan_TE 的宿主侧 gate 是对 discovery-only 输出的必要补偿，并通过 `pan_te_handoff.json` 审计；它不等同于完整 dtr 的内部 G1–G6 审计。

## 14. 代码溯源索引

本文各行为可在以下文件追溯：

- 总调度、检查点、dtr 内存/门控/合并：`bin/Pan_TE`
- mdl 包装、Dust/TRF、>2 Gb 采样：`bin/build_mdl`
- Refiner 调度与参数：`bin/Refiner_mdl/main.py`、`config.py`
- Refiner 的 Phase 0/1/2：`phase0_triage.py`、`phase1_consensus.py`、`phase2_library_split.py`
- refine 的设计目标与实施差异：`bin/Refiner_mdl/REFINE_STRATEGY_DESIGN_v2.md`
- PAF 解析、链、CIGAR、QC、scaffold：`bin/paf_processor.py`
- GFA 编排：`bin/gfa_processor.py`
- PGGB P 路径等位算法：`bin/gfa_alleles.py`
- rGFA+GAF 验证、锚链、presence：`bin/rgfa_gaf.py`
- 图节点注释：`bin/gfa_annotator.py`
- 片段共线重组：`bin/fragment_stitcher.py`
