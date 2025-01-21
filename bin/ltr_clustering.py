#!/usr/bin/env python3

import networkx as nx
from Bio import SeqIO
from collections import defaultdict
import subprocess
import os
import tempfile
from itertools import islice
import multiprocessing
import shutil
from Refiner import ConservativeMSAConsensusBuilder
import random
from multiprocessing import Pool

def batch_iterator(iterator, batch_size):
    while True:
        batch = list(islice(iterator, batch_size))
        if not batch:
            break
        yield batch

class MMseqsRunner:
    def __init__(self, threads=None):
        self.threads = threads or multiprocessing.cpu_count()
        self.FAMILY_SIZE_CUTOFF = 15
        self.MIN_SCORE = 250
        self.GAP_INIT = -25
        self.GAP_EXT = -5
        self.MIN_MATCH = 7
        self.MAX_ELEMENTS = 100

        try:
            devnull = open(os.devnull, 'w')
            subprocess.run(['mmseqs'], stdout=devnull, stderr=devnull)
            devnull.close()
        except:
            raise RuntimeError("MMseqs2 not found. Please install MMseqs2 first.")

    def run_search(self, input_fasta, output_dir, min_seq_id=0.8):
        tmp_dir = os.path.join(output_dir, 'tmp')
        os.makedirs(tmp_dir, exist_ok=True)

        db_path = os.path.join(output_dir, 'DB')
        result_path = os.path.join(output_dir, 'res')
        result_file = os.path.join(output_dir, 'search_results.txt')

        try:
            print("Creating MMseqs2 database...")
            retcode = subprocess.call([
                'mmseqs', 'createdb',
                input_fasta,
                db_path
            ])
            if retcode != 0:
                raise RuntimeError("Database creation failed")

            print("Running MMseqs2 search...")
            retcode = subprocess.call([
                'mmseqs', 'search',
                db_path, db_path,
                result_path,
                tmp_dir,
                '--threads', str(self.threads),
                '--min-seq-id', str(min_seq_id),
                '--search-type', '3',
                '--min-ungapped-score', '250',
                '--gap-open', '25',
                '--gap-extend', '5',
                '-s', '5.7',
                '--max-seqs', '1000',
                '--min-length', str(self.MIN_MATCH),
                '--cov-mode', '3',
                '--mask-lower-case', '1',
                '--e-profile', '0.001'
            ])
            if retcode != 0:
                raise RuntimeError("Search failed")

            print("Converting results to readable format...")
            retcode = subprocess.call([
                'mmseqs', 'convertalis',
                db_path, db_path,
                result_path,
                result_file,
                '--format-output', 'query,target,pident,alnlen,qlen,tlen'
            ])
            if retcode != 0:
                raise RuntimeError("Results conversion failed")

            return result_file

        except:
            print("Error in MMseqs2 processing")
            raise

        finally:
            if os.path.exists(tmp_dir):
                shutil.rmtree(tmp_dir)
            for ext in ['', '.index', '.dbtype', '.lookup']:
                db_file = db_path + ext
                if os.path.exists(db_file):
                    os.remove(db_file)
            for ext in ['', '.index', '.dbtype']:
                res_file = result_path + ext
                if os.path.exists(res_file):
                    os.remove(res_file)

def parse_mmseqs_results(results_file, min_identity=60, min_coverage=0.3, batch_size=10000):
    edges = []
    try:
        with open(results_file) as f:
            while True:
                batch = list(islice(f, batch_size))
                if not batch:
                    break
                for line in batch:
                    qseqid, sseqid, pident, length, qlen, slen = line.strip().split('\t')
                    pident = float(pident)
                    length = int(length)
                    qlen = int(qlen)
                    slen = int(slen)
                    coverage = length / min(qlen, slen)
                    if pident >= min_identity or coverage >= min_coverage:
                        if qseqid < sseqid:
                            edges.append((qseqid, sseqid, {'weight': pident}))
                if len(edges) >= batch_size:
                    yield edges
                    edges = []
        if edges:
            yield edges
    except:
        print("Error parsing MMseqs2 results")
        raise

def build_similarity_network(edge_batches):
    G = nx.Graph()
    for edges in edge_batches:
        G.add_edges_from(edges)
    return G

def cluster_sequences(G, min_family_size=15):
    communities = {}
    for i, comp in enumerate(nx.connected_components(G)):
        if len(comp) >= min_family_size:
            for node in comp:
                communities[node] = i
    return communities

def parallel_process_family(args):
    sequences, family_id, output_dir, threads_per_family = args
    return process_family_sequences(sequences, family_id, output_dir, threads_per_family)

def process_family_sequences(sequences, family_id, output_dir, threads=None):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    temp_input = os.path.join(output_dir, f"family_{family_id}_temp_input.fasta")
    temp_output = os.path.join(output_dir, f"family_{family_id}_consensus.fasta")

    with open(temp_input, 'w') as f:
        for seq_id, seq in sequences:
            f.write(f">{seq_id}\n{seq}\n")

    builder = ConservativeMSAConsensusBuilder(
        min_similarity=0.4,
        min_coverage=0.2,
        min_base_freq=0.3,
        threads=threads
    )

    try:
        builder.build_consensus(temp_input, temp_output)
        return temp_output
    except:
        print(f"Error processing family {family_id}")
        return None
    finally:
        if os.path.exists(temp_input):
            os.remove(temp_input)

def process_families_parallel(family_sequences, output_dir, total_threads, min_family_size=15):
    valid_families = []
    for family_id, sequences in family_sequences.items():
        if len(sequences) >= min_family_size:
            if len(sequences) > 100:
                sequences = sorted(sequences, key=lambda x: x[2], reverse=True)[:100]
            valid_families.append((sequences, family_id, output_dir, max(1, total_threads // 4)))

    if not valid_families:
        return []

    num_processes = min(len(valid_families), max(1, total_threads // 4))
    print(f"Processing {len(valid_families)} families with {num_processes} parallel processes...")

    with Pool(processes=num_processes) as pool:
        consensus_files = pool.map(parallel_process_family, valid_families)

    return [f for f in consensus_files if f is not None]

def main(fasta_file, output_dir, min_identity=80, min_coverage=0.5,
         batch_size=10000, threads=None, max_family_size=100, min_family_size=15):
    if threads is None:
        threads = multiprocessing.cpu_count()

    os.makedirs(output_dir, exist_ok=True)

    mmseqs_runner = MMseqsRunner(threads)
    search_results = mmseqs_runner.run_search(
        fasta_file,
        output_dir,
        min_seq_id=min_identity/100
    )

    print("Parsing MMseqs2 results and building similarity network...")
    edge_batches = parse_mmseqs_results(search_results, min_identity, min_coverage, batch_size)
    G = build_similarity_network(edge_batches)

    print("Clustering sequences...")
    communities = cluster_sequences(G, min_family_size)

    family_sequences = defaultdict(list)

    for batch in batch_iterator(SeqIO.parse(fasta_file, "fasta"), batch_size):
        for record in batch:
            if record.id in communities:
                family_id = communities[record.id]
                family_sequences[family_id].append((record.id, str(record.seq)))

    if max_family_size:
        for family_id, sequences in family_sequences.items():
            if len(sequences) > max_family_size:
                family_sequences[family_id] = random.sample(sequences, max_family_size)

    consensus_files = process_families_parallel(
        family_sequences,
        output_dir,
        threads,
        min_family_size
    )

    final_output = os.path.join(output_dir, "all_consensus_sequences.fasta")
    with open(final_output, 'w') as outfile:
        for consensus_file in consensus_files:
            if os.path.exists(consensus_file):
                with open(consensus_file) as infile:
                    outfile.write(infile.read())
                os.remove(consensus_file)

    if os.path.exists(search_results):
        os.remove(search_results)

    return final_output

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Cluster LTR sequences and generate consensus sequences')
    parser.add_argument('--fasta', required=True, help='Input FASTA file')
    parser.add_argument('--output', required=True, help='Output directory')
    parser.add_argument('--min-identity', type=float, default=80,
                      help='Minimum sequence identity (default: 80)')
    parser.add_argument('--min-coverage', type=float, default=0.5,
                      help='Minimum sequence coverage (default: 0.5)')
    parser.add_argument('--batch-size', type=int, default=10000,
                      help='Batch size for processing (default: 10000)')
    parser.add_argument('--threads', type=int, default=None,
                      help='Number of threads (default: all available)')
    parser.add_argument('--max-family-size', type=int, default=100,
                      help='Maximum number of sequences per family (default: 100)')
    parser.add_argument('--min-family-size', type=int, default=15,
                      help='Minimum number of sequences required to form a family (default: 15)')

    args = parser.parse_args()

    final_output = main(args.fasta, args.output,
                       args.min_identity, args.min_coverage,
                       args.batch_size, args.threads,
                       args.max_family_size, args.min_family_size)

    print(f"All consensus sequences have been written to: {final_output}")