#!/usr/bin/env python3
coverage_file = "repeat_coverage_result.txt"
lib_fasta = "lib.fa"
output_fasta = "filtered_lib.fa"

# Step 1: 读取coverage信息
coverage_dict = {}
with open(coverage_file, "r") as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        # 格式：repeat_name    cov_str（cov_str为"1,2,3..."）
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        repeat_name = parts[0]
        cov_str = parts[1]
        cov_list = list(map(int, cov_str.split(",")))
        coverage_dict[repeat_name] = cov_list

# Step 2: 读取lib.fa中的序列
# 假设lib.fa是标准FASTA格式，形如：
# >repeat_name
# AGCTAGCTAGC...
# >another_repeat
# ...
# 我们将其以字典形式读入：{repeat_name: sequence}
sequences = {}
current_name = None
current_seq_lines = []
with open(lib_fasta, "r") as f:
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            # 遇到新的序列名，先将上一个序列存入
            if current_name is not None:
                sequences[current_name] = "".join(current_seq_lines)
            # 初始化新的记录
            current_name = line[1:].strip()  # 去掉">"
            current_seq_lines = []
        else:
            # 序列行
            current_seq_lines.append(line)

# 不要忘记最后一条序列
if current_name is not None:
    sequences[current_name] = "".join(current_seq_lines)

# Step 3: 根据coverage对序列进行过滤
# coverage < 3 的碱基切除，只保留 coverage≥3 的碱基
with open(output_fasta, "w") as out:
    for rname, seq in sequences.items():
        if rname not in coverage_dict:
            # 如果没有对应的coverage信息，可选择跳过或原样输出
            # 此处选择跳过或原样输出自行决定，这里假设跳过：
            continue
        
        cov_list = coverage_dict[rname]
        
        # 序列长度和coverage长度可能存在不一致的情况，需要注意
        # 假设coverage长度 >= 序列长度（通常是这样），如果不是则截断
        min_len = min(len(seq), len(cov_list))
        
        filtered_bases = []
        for i in range(min_len):
            if cov_list[i] >= 3:
                filtered_bases.append(seq[i])
        
        # 拼接过滤后的碱基
        filtered_seq = "".join(filtered_bases)

        # 如果过滤后序列不为空则输出
        if filtered_seq:
            out.write(f">{rname}\n")
            # 序列一般FASTA格式可每行60-80bp，如需要可做分行处理，这里简单直接一行输出
            # 如果需要分行，每60个碱基一行，可额外处理
            out.write(filtered_seq + "\n")

print("切割完成，结果已写入：", output_fasta)

