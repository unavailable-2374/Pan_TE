#!/usr/bin/env python3
"""Build the cellular protein BLAST DB used by final_gene_filter.py.

Input is a Swiss-Prot FASTA already restricted by the UniProt query
'(reviewed:true) NOT (keyword:KW-0814) NOT (taxonomy_id:10239)'.

- Proteins are removed when their description names a TE product, or when RepeatPeps
  homology (DIAMOND e <= 1e-10) covers at least 30% of the protein.
- Remaining proteins can still carry TE-derived segments: an Alu has no protein of its own,
  so RepeatPeps cannot flag a host protein that translates one, and such a segment would
  make the filter cut real Alu entries. Segments homologous to RepeatPeps, or to the
  translated Dfam curated TE consensi (identity >= 35%), are replaced by X. The rest of the
  protein stays usable. Dfam segments are located with DIAMOND and with BLASTX at the same
  settings the gene filter uses, so the mask is at least as sensitive as the filter.

Usage: build_cellular_db.py <query.fa> <RepeatPeps.lib> <Dfam-curated_only.embl.gz> <output_dir> <threads>
"""
import collections
import gzip
import os
import re
import subprocess
import sys

TE_WORDS = re.compile(
    r"transpos|retrotranspos|reverse transcriptase|integrase|\bgag\b|gag-pol|pol polyprotein|"
    r"copia|gypsy|\bTy[1-5]\b|LINE-1|ORF1p|ORF2p|piggyBac|mariner|harbinger|\bPIF\b|mutator|"
    r"\bhAT\b|CACTA|helitron|\bTc1\b|retrovir|endogenous retrovirus|\bERV|tigger|PGBD|jockey", re.I)
NON_TE_TYPES = {"Simple_repeat", "Low_complexity", "Satellite", "rRNA", "tRNA", "snRNA", "scRNA",
                "srpRNA", "RNA", "ARTEFACT", "Segmental", "Other"}
MIN_IDENTITY = "35"


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


def embl_te_consensi(lines, out) -> int:
    """Write the TE families of a Dfam EMBL stream as FASTA; returns the number written."""
    written = 0
    name, te_type, seq, in_seq = None, None, [], False
    for line in lines:
        if line.startswith("ID"):
            name, te_type, seq, in_seq = line.split()[1].rstrip(";"), None, [], False
        elif "Type:" in line and line.startswith("CC") and te_type is None:
            te_type = line.split("Type:")[1].split()[0] if line.split("Type:")[1].split() else None
        elif line.startswith("SQ"):
            in_seq = True
        elif line.startswith("//"):
            if name and seq and te_type and te_type not in NON_TE_TYPES:
                out.write(f">{name}\n{''.join(seq)}\n")
                written += 1
            name, in_seq = None, False
        elif in_seq:
            seq.extend(part for part in line.split() if part.isalpha())
    return written


def mask_sequence(seq: str, intervals) -> str:
    """Replace 1-based inclusive intervals with X."""
    chars = list(seq)
    for start, end in intervals:
        for i in range(max(start, 1) - 1, min(end, len(chars))):
            chars[i] = "X"
    return "".join(chars)


def write_masked_fasta(source: str, target: str, mask) -> tuple:
    masked_proteins = masked_aa = total = 0
    with open(target, "w") as out:
        header, seq = None, []

        def flush():
            nonlocal masked_proteins, masked_aa, total
            if header is None:
                return
            total += 1
            protein = header[1:].split()[0]
            text = "".join(seq)
            if protein in mask:
                new = mask_sequence(text, mask[protein])
                masked_proteins += 1
                masked_aa += new.count("X") - text.count("X")
                text = new
            out.write(f"{header}\n{text}\n")

        for line in open(source):
            if line.startswith(">"):
                flush()
                header, seq = line.rstrip("\n"), []
            else:
                seq.append(line.strip())
        flush()
    return total, masked_proteins, masked_aa


def blastx_segments(query: str, db: str, threads: int, work_dir: str):
    """Protein intervals hit by translated TE consensi, at the gene filter's BLASTX settings."""
    records, name = {}, None
    for line in open(query):
        if line.startswith(">"):
            name = line[1:].split()[0]; records[name] = []
        else:
            records[name].append(line.strip())
    chunks = [os.path.join(work_dir, f"dfam_q{j}.fa") for j in range(threads)]
    handles = [open(c, "w") for c in chunks]
    for i, (k, v) in enumerate(records.items()):
        handles[i % threads].write(f">{k}\n{''.join(v)}\n")
    for h in handles:
        h.close()

    def run(chunk):
        out = chunk + ".tsv"
        subprocess.run(["blastx", "-query", chunk, "-db", db, "-evalue", "1e-5", "-max_target_seqs", "100000",
                        "-num_threads", "1", "-outfmt", "6 sseqid sstart send pident", "-out", out], check=True)
        return out

    import concurrent.futures
    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as pool:
        outs = list(pool.map(run, chunks))
    for out in outs:
        for line in open(out):
            protein, start, end, pident = line.split()
            if float(pident) >= float(MIN_IDENTITY):
                yield min(int(start), int(end)), max(int(start), int(end)), protein
    for c in chunks + outs:
        os.remove(c)


def diamond(*args) -> None:
    subprocess.run(["diamond", *args, "--quiet"], check=True)


def main(query: str, peps: str, dfam_embl: str, out_dir: str, threads: str) -> None:
    os.makedirs(out_dir, exist_ok=True)
    path = lambda name: os.path.join(out_dir, name)
    copy_fasta(query, path("sprot_kw.fa"), lambda h: not names_te_product(h))

    diamond("makedb", "--in", peps, "-d", path("peps"), "-p", threads)
    diamond("blastp", "-q", path("sprot_kw.fa"), "-d", path("peps"), "-o", path("vs_peps.tsv"), "-p", threads,
            "--evalue", "1e-10", "--query-cover", "30", "--max-target-seqs", "1", "--outfmt", "6", "qseqid")
    homologous = {line.strip() for line in open(path("vs_peps.tsv"))}
    kept = copy_fasta(path("sprot_kw.fa"), path("sprot_unmasked.fa"),
                      lambda h: h[1:].split()[0] not in homologous)
    if kept == 0:
        raise RuntimeError("no proteins left after filtering")

    mask = collections.defaultdict(list)
    diamond("blastp", "-q", path("sprot_unmasked.fa"), "-d", path("peps"), "-o", path("seg_peps.tsv"), "-p", threads,
            "--evalue", "1e-10", "--max-target-seqs", "25", "--outfmt", "6", "qseqid", "qstart", "qend")
    for line in open(path("seg_peps.tsv")):
        protein, start, end = line.split()
        mask[protein].append((int(start), int(end)))

    with gzip.open(dfam_embl, "rt") as lines, open(path("dfam_te.fa"), "w") as out:
        families = embl_te_consensi(lines, out)
    if families == 0:
        raise RuntimeError(f"no TE families parsed from {dfam_embl}")
    diamond("makedb", "--in", path("sprot_unmasked.fa"), "-d", path("sprot"), "-p", threads)
    diamond("blastx", "-q", path("dfam_te.fa"), "-d", path("sprot"), "-o", path("seg_dfam.tsv"), "-p", threads,
            "--more-sensitive", "--evalue", "1e-5", "--id", MIN_IDENTITY, "--max-target-seqs", "0",
            "--outfmt", "6", "sseqid", "sstart", "send")
    for line in open(path("seg_dfam.tsv")):
        protein, start, end = line.split()
        mask[protein].append(tuple(sorted((int(start), int(end)))))

    subprocess.run(["makeblastdb", "-in", path("sprot_unmasked.fa"), "-dbtype", "prot", "-out", path("sprot_unmasked"),
                    "-parse_seqids"], check=True, stdout=subprocess.DEVNULL)
    for start, end, protein in blastx_segments(path("dfam_te.fa"), path("sprot_unmasked"), int(threads), out_dir):
        mask[protein].append((start, end))

    total, masked_proteins, masked_aa = write_masked_fasta(path("sprot_unmasked.fa"), path("sprot_cell.fa"), mask)
    subprocess.run(["makeblastdb", "-in", path("sprot_cell.fa"), "-dbtype", "prot", "-out", path("sprot_cell"),
                    "-parse_seqids"], check=True, stdout=subprocess.DEVNULL)
    for name in ("sprot_kw.fa", "sprot_unmasked.fa", "vs_peps.tsv", "seg_peps.tsv", "seg_dfam.tsv",
                 "dfam_te.fa", "peps.dmnd", "sprot.dmnd"):
        os.remove(path(name))
    for name in os.listdir(out_dir):
        if name.startswith("sprot_unmasked."):
            os.remove(path(name))
    print(f"sprot_cell: {total} proteins, {families} Dfam TE families, "
          f"{masked_proteins} proteins with {masked_aa} residues masked, in {out_dir}")


if __name__ == "__main__":
    if len(sys.argv) != 6:
        sys.exit(__doc__)
    main(*sys.argv[1:])
