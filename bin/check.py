#!/usr/bin/env python3
# 导入需要的库
import sys

# 用于存储读段ID和其比对次数
read_counts = {}

# 打开SAM文件进行读取
with open('LTR.sam', 'r') as file:
    for line in file:
        if line.startswith('@'):  # 跳过文件头部的注释行
            continue
        fields = line.split('\t')  # 按制表符分割行
        flag = int(fields[1])  # FLAG位是每行的第二个字段
        if flag & 4:  # 检查FLAG位，判断是否未比对上
            continue  # 如果未比对上，则跳过当前读段
        read_id = fields[0]  # 读段的ID是每行的第一个字段
        if read_id in read_counts:
            read_counts[read_id] += 1  # 如果ID已存在，比对次数加1
        else:
            read_counts[read_id] = 1  # 如果是新ID，则初始化比对次数为1

# 筛选出比对次数2次及以上的读段ID
filtered_ids = [read_id for read_id, count in read_counts.items() if count >= 1]

# 输出结果
for read_id in filtered_ids:
    print(read_id)

