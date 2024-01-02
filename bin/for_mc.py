#!/usr/bin/env python3
import random
import argparse

def read_gfa_file(gfa_file):
    segments = {}
    links = {}
    segment_lengths = {}  # 用于存储每个片段的长度

    with open(gfa_file, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            if parts[0] == 'S':
                segment_id, sequence = parts[1], parts[2]
                segments[segment_id] = sequence
                segment_lengths[segment_id] = len(sequence)  # 计算并存储片段长度
                links[segment_id] = ([], [])  # 初始化每个片段的连接列表

            elif parts[0] == 'L':
                from_seg, from_dir, to_seg, to_dir = parts[1], parts[2], parts[3], parts[4]
                if from_seg not in links:
                    links[from_seg] = ([], [])
                if to_seg not in links:
                    links[to_seg] = ([], [])

                if from_dir == '+':
                    links[from_seg][1].append(to_seg)
                else:
                    links[from_seg][0].append(to_seg)

                if to_dir == '+':
                    links[to_seg][0].append(from_seg)
                else:
                    links[to_seg][1].append(from_seg)

    return segments, links, segment_lengths

def extend_segment(segment_id, links, segment_lengths, used_segments, min_length, max_length, core_uniq_length):
    path = [segment_id]
    current_length = segment_lengths.get(segment_id, 0)
    used_length = 0
    local_used_segments = set()

    while current_length < max_length:
        extended = False
        for direction in [0, 1]:
            connection_options = links.get(path[0] if direction == 0 else path[-1])
            if connection_options and len(connection_options) > direction:
                options = connection_options[direction]
                choice = next((seg for seg in options if seg not in used_segments and seg not in local_used_segments), None)

                if not choice:
                    choice = next((seg for seg in options if seg in used_segments or seg in local_used_segments), None)

                if choice:
                    path.insert(0, choice) if direction == 0 else path.append(choice)
                    segment_length = segment_lengths.get(choice, 0)
                    current_length += segment_length
                    if choice in used_segments or choice in local_used_segments:
                        used_length += segment_length
                    else:
                        local_used_segments.add(choice)
                    extended = True
                    break

        if not extended:
            break

        # 检查在刚达到 min_length 时的已使用片段占比
        if current_length >= min_length:
            if (used_length / current_length) > core_uniq_length:
                return None  # 如果占比超过40%，延展失败
            break  # 如果达到 min_length 且占比未超过40%，停止延展

    return path if current_length >= min_length else None

def process_ids_only(gfa_file, output_fasta, output_order, min_length, core_min_length, core_max_length, max_length, core_uniq_length):
    # 读取GFA文件
    segments, links, segment_lengths = read_gfa_file(gfa_file)

    # 初始化核心片段集合
    core_segments = {sid for sid, length in segment_lengths.items() if core_min_length <= length <= core_max_length}
    
    # 用于跟踪哪些片段已被使用
    used_segments = set()
    
    # 存储成功延展的路径
    successful_paths = []
    num = 0

    # 循环处理每个核心片段
    while core_segments:
        segment_id = random.choice(list(core_segments))
        core_segments.remove(segment_id)
        # 执行延展逻辑
        path = extend_segment(segment_id, links, segment_lengths, used_segments, min_length, max_length, core_uniq_length)
        if path is None:  # 放弃这次延展
            print ("Failed")
            used_segments.add(segment_id)
            continue  # 直接进入下一次循环
        
        if path:
            print ("Successed")
            # 如果延展成功，记录路径并标记所有片段为已使用
            successful_paths.append(path)
            used_segments.update(path)
            core_segments.difference_update(path)

    with open(output_order, 'w') as order_file, open(output_fasta, 'w') as fasta_file:
        for path in successful_paths:
            orders = ','.join(path)
            num=num+1
            order_file.write(f'Segment_{num}\t{orders}\n')
            # 通过ID与片段的对应关系组成序列
            sequence = ''.join(segments[sid] for sid in path if sid in segments)
            fasta_file.write(f'>Segment_{num}\n{sequence}\n')

        # 输出未参与延展且长度超过40kb的序列
        for segment_id, length in segment_lengths.items():
            if segment_id not in used_segments and length > min_length:
                num=num+1
                fasta_file.write(f'>Segment_{num}\n{segments[segment_id]}\n')
    
    return successful_paths

def main(gfa_file, output_fasta, output_order, min_length, core_min_length, core_max_length, max_length, core_uniq_length):
    process_ids_only(gfa_file, output_fasta, output_order, min_length, core_min_length, core_max_length, max_length, core_uniq_length)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process GFA file for segment extension.")
    parser.add_argument("gfa_file", type=str, help="Path to the GFA file.")
    parser.add_argument("output_fasta", type=str, help="Path to the output FASTA file.")
    parser.add_argument("output_order", type=str, help="Path to the output order file.")
    parser.add_argument("--min_length", type=int, default=40000, help="Minimum length of extended segments.")
    parser.add_argument("--max_length", type=int, default=40000, help="Maximum length of extended segments.")
    parser.add_argument("--core_min_length", type=int, default=500, help="Minimum length for a segment to be considered as core.")
    parser.add_argument("--core_max_length", type=int, default=40000, help="Maximum length for a segment to be considered as core.")
    parser.add_argument("--core_uniq_length", type=float, default=0.4, help="Minimum length for a uniq sequence on core.")
    args = parser.parse_args()

    main(args.gfa_file, args.output_fasta, args.output_order, args.min_length, args.core_min_length, args.core_max_length, args.max_length, args.core_uniq_length)

