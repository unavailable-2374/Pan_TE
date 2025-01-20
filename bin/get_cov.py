#!/usr/bin/env python3
input_file = "tmp.fa.ori.out"
output_file = "repeat_coverage_result.txt"

coverage = {}  # 用于存储coverage信息, key是repeat名称，value是列表，索引为碱基位置(0-based)，值为覆盖次数
repeat_lengths = {}  # 用于记录每个repeat的consensus长度

with open(input_file, "r") as f:
    for line in f:
        line = line.strip()
        # 跳过空行和注释行
        if not line or line.startswith("#"):
            continue

        parts = line.split()
        # 检查列数，标准RepeatMasker输出行通常有>=14列数据
        if len(parts) < 14:
            continue
        
        # 典型列：
        # 0: SW score
        # 1: %div
        # 2: %del
        # 3: %ins
        # 4: query_sequence
        # 5: q_begin
        # 6: q_end
        # 7: q_(left)
        # 8: strand (+/-/C)
        # 9: repeat_name
        # 10: class/family
        # 11: rep_begin
        # 12: rep_end
        # 13: (rep_left)
        
        strand = parts[8]
        repeat_name = parts[9]
        rep_begin = int(parts[11])
        rep_end = int(parts[12])

        # 去掉rep_left字段的括号
        rep_left_str = parts[13].strip("()")
        rep_left = int(rep_left_str)

        # 计算consensus长度
        consensus_length = rep_end + rep_left

        # 初始化或扩展coverage数组
        if repeat_name not in coverage:
            coverage[repeat_name] = [0] * consensus_length
            repeat_lengths[repeat_name] = consensus_length
        else:
            if consensus_length > repeat_lengths[repeat_name]:
                old_len = repeat_lengths[repeat_name]
                coverage[repeat_name].extend([0] * (consensus_length - old_len))
                repeat_lengths[repeat_name] = consensus_length

        # 计数覆盖度
        for i in range(rep_begin - 1, rep_end):
            coverage[repeat_name][i] += 1

# 输出结果
# 格式：repeat_name\tcov1,cov2,cov3,... (一行一个repeat)
with open(output_file, "w") as out:
    for rname, cov_list in coverage.items():
        coverage_str = ",".join(map(str, cov_list))
        out.write(f"{rname}\t{coverage_str}\n")

print("统计完成，结果已写入：", output_file)

