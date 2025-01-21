#!/usr/bin/env python3
import sys

read_counts = {}

# 打开SAM文件进行读取
with open('LTR.sam', 'r') as file:
    for line in file:
        if line.startswith('@'): 
            continue
        fields = line.split('\t')  
        flag = int(fields[1])  
        if flag & 4:  
            continue  
        read_id = fields[0] 
        if read_id in read_counts:
            read_counts[read_id] += 1  
        else:
            read_counts[read_id] = 1  

filtered_ids = [read_id for read_id, count in read_counts.items() if count >= 1]

for read_id in filtered_ids:
    print(read_id)

