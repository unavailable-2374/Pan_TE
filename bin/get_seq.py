#!/usr/bin/env python3
coverage_file = "repeat_coverage_result.txt"
lib_fasta = "lib.fa"
output_fasta = "filtered_lib.fa"

coverage_dict = {}
with open(coverage_file, "r") as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        repeat_name = parts[0]
        cov_str = parts[1]
        cov_list = list(map(int, cov_str.split(",")))
        coverage_dict[repeat_name] = cov_list

sequences = {}
current_name = None
current_seq_lines = []
with open(lib_fasta, "r") as f:
    for line in f:
        line = line.strip()
        if line.startswith(">"):
    
            if current_name is not None:
                sequences[current_name] = "".join(current_seq_lines)
            
            current_name = line[1:].strip()  
            current_seq_lines = []
        else:
            
            current_seq_lines.append(line)

if current_name is not None:
    sequences[current_name] = "".join(current_seq_lines)


with open(output_fasta, "w") as out:
    for rname, seq in sequences.items():
        if rname not in coverage_dict:
            continue
        
        cov_list = coverage_dict[rname]
        
        min_len = min(len(seq), len(cov_list))
        
        filtered_bases = []
        for i in range(min_len):
            if cov_list[i] >= 3:
                filtered_bases.append(seq[i])
        
        filtered_seq = "".join(filtered_bases)

        if filtered_seq:
            out.write(f">{rname}\n")
            out.write(filtered_seq + "\n")

print("Complete, Write to：", output_fasta)

