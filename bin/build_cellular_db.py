#!/usr/bin/env python3
"""Build the cellular protein BLAST DB used by final_gene_filter.py.

Input is a Swiss-Prot FASTA already restricted by the UniProt query
'(reviewed:true) NOT (keyword:KW-0814) NOT (taxonomy_id:10239)'. Proteins are removed when
their description names a TE product, or when they are homologous to RepeatPeps (DIAMOND
e <= 1e-10 covering >= 30% of the protein).

Usage: build_cellular_db.py <query.fa> <RepeatPeps.lib> <output_dir> <threads>
"""
import os
import re
import subprocess
import sys

TE_WORDS = re.compile(
    r"transpos|retrotranspos|reverse transcriptase|integrase|\bgag\b|gag-pol|pol polyprotein|"
    r"copia|gypsy|\bTy[1-5]\b|LINE-1|ORF1p|ORF2p|piggyBac|mariner|harbinger|\bPIF\b|mutator|"
    r"\bhAT\b|CACTA|helitron|\bTc1\b|retrovir|endogenous retrovirus|\bERV|tigger|PGBD|jockey", re.I)


def names_te_product(header: str) -> bool:
    description = header.split(" ", 1)[1] if " " in header else ""
    return bool(TE_WORDS.search(description.split(" OS=")[0]))


def copy_fasta(source: str, target: str, keep) -> int:
    kept, write = 0, False
    with open(target, "w") as out:
        for line in open(source):
            if line.startswith(">"):
                write = keep(line)
                kept += write
            if write:
                out.write(line)
    return kept


def main(query: str, peps: str, out_dir: str, threads: str) -> None:
    os.makedirs(out_dir, exist_ok=True)
    no_keyword = os.path.join(out_dir, "sprot_kw.fa")
    final = os.path.join(out_dir, "sprot_cell.fa")
    hits_path = os.path.join(out_dir, "vs_peps.tsv")
    db = os.path.join(out_dir, "peps")
    copy_fasta(query, no_keyword, lambda h: not names_te_product(h))
    subprocess.run(["diamond", "makedb", "--in", peps, "-d", db, "-p", threads, "--quiet"], check=True)
    subprocess.run(["diamond", "blastp", "-q", no_keyword, "-d", db, "-o", hits_path, "-p", threads,
                    "--evalue", "1e-10", "--query-cover", "30", "--max-target-seqs", "1",
                    "--outfmt", "6", "qseqid", "--quiet"], check=True)
    homologous = {line.strip() for line in open(hits_path)}
    kept = copy_fasta(no_keyword, final, lambda h: h[1:].split()[0] not in homologous)
    if kept == 0:
        raise RuntimeError("no proteins left after filtering")
    subprocess.run(["makeblastdb", "-in", final, "-dbtype", "prot", "-out",
                    os.path.join(out_dir, "sprot_cell"), "-parse_seqids"],
                   check=True, stdout=subprocess.DEVNULL)
    for path in (no_keyword, hits_path, db + ".dmnd"):
        os.remove(path)
    print(f"sprot_cell: {kept} proteins in {out_dir}")


if __name__ == "__main__":
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    main(*sys.argv[1:])
