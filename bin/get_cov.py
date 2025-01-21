#!/usr/bin/env python3
input_file = "tmp.fa.ori.out"
output_file = "repeat_coverage_result.txt"

coverage = {} 
repeat_lengths = {}  

with open(input_file, "r") as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith("#"):
            continue

        parts = line.split()
        if len(parts) < 14:
            continue
        
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

        rep_left_str = parts[13].strip("()")
        rep_left = int(rep_left_str)

        consensus_length = rep_end + rep_left

        if repeat_name not in coverage:
            coverage[repeat_name] = [0] * consensus_length
            repeat_lengths[repeat_name] = consensus_length
        else:
            if consensus_length > repeat_lengths[repeat_name]:
                old_len = repeat_lengths[repeat_name]
                coverage[repeat_name].extend([0] * (consensus_length - old_len))
                repeat_lengths[repeat_name] = consensus_length

        for i in range(rep_begin - 1, rep_end):
            coverage[repeat_name][i] += 1

# Formate：repeat_name\tcov1,cov2,cov3,... 
with open(output_file, "w") as out:
    for rname, cov_list in coverage.items():
        coverage_str = ",".join(map(str, cov_list))
        out.write(f"{rname}\t{coverage_str}\n")

print("Complete, Write to：", output_file)

