#!/usr/bin/env python3
"""Drop library entries that are fragments contained in a longer entry of the same class.

The cross-detector merge in the Combine step (cd-hit-est -c 0.9 -aS 0.8 -aL 0.8)
only collapses near-equal-length entries, so a short consensus that is fully
covered by a longer one survives. These contained fragments (mostly mdl-repeat
pieces of LTR elements) make up most of the library's redundancy while adding
almost nothing to masking.

A second cd-hit-est pass clusters entries whose shorter member is covered
>=80% at >=80% identity by the longer representative. Within each cluster a
member is dropped when it has the representative's top-level class or is
Unknown. Members of a different class are kept: these are usually nested
insertions (for example a SINE or DNA element inside an LTR), and dropping them
would lose that class from the annotation.

Usage: drop_contained_fragments.py <in.fa> <out.fa> [--threads N] [--report report.tsv]
"""
import argparse
import os
import shutil
import subprocess
import sys
import tempfile

CDHIT_ARGS = ['-c', '0.8', '-aS', '0.8', '-G', '0', '-n', '5', '-r', '1',
              '-mask', 'NX', '-M', '0', '-d', '0']


def top_class(seq_id):
    """Top-level class from a RepeatMasker-style id ('name#LTR/Gypsy' -> 'LTR')."""
    if '#' not in seq_id:
        return 'Unknown'
    label = seq_id.split('#', 1)[1].split('/', 1)[0]
    return label or 'Unknown'


def parse_clusters(clstr_lines):
    """Yield (representative, [members]) from cd-hit .clstr lines."""
    rep, members = None, []
    for line in clstr_lines:
        if line.startswith('>Cluster'):
            if rep is not None:
                yield rep, members
            rep, members = None, []
            continue
        if '>' not in line:
            continue
        seq_id = line.split('>', 1)[1].split('...', 1)[0]
        if line.rstrip().endswith('*'):
            rep = seq_id
        else:
            members.append(seq_id)
    if rep is not None:
        yield rep, members


def select_drops(clusters):
    """Return {member: representative} for members to drop, plus kept cross-class pairs."""
    drops, kept = {}, []
    for rep, members in clusters:
        rep_class = top_class(rep)
        for member in members:
            member_class = top_class(member)
            if member_class == 'Unknown' or member_class == rep_class:
                drops[member] = rep
            else:
                kept.append((member, rep))
    return drops, kept


def filter_fasta(in_path, out_path, drops):
    kept = total = 0
    with open(in_path) as src, open(out_path, 'w') as dst:
        keep = False
        for line in src:
            if line.startswith('>'):
                total += 1
                keep = line[1:].split()[0] not in drops
                kept += keep
            if keep:
                dst.write(line)
    return total, kept


def run(in_path, out_path, threads=1, report=None):
    if shutil.which('cd-hit-est') is None:
        raise RuntimeError('cd-hit-est not found in PATH')
    with tempfile.TemporaryDirectory(dir=os.path.dirname(os.path.abspath(out_path))) as tmp:
        clustered = os.path.join(tmp, 'clustered.fa')
        subprocess.run(['cd-hit-est', '-i', in_path, '-o', clustered, '-T', str(threads)] + CDHIT_ARGS,
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
        with open(clustered + '.clstr') as fh:
            drops, kept_cross = select_drops(parse_clusters(fh))
    total, kept = filter_fasta(in_path, out_path, drops)
    if report:
        with open(report, 'w') as fh:
            fh.write('member\trepresentative\taction\n')
            for member, rep in sorted(drops.items()):
                fh.write(f'{member}\t{rep}\tdropped\n')
            for member, rep in kept_cross:
                fh.write(f'{member}\t{rep}\tkept_cross_class\n')
    return {'input': total, 'kept': kept, 'dropped': len(drops), 'kept_cross_class': len(kept_cross)}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n', 1)[0])
    parser.add_argument('input')
    parser.add_argument('output')
    parser.add_argument('--threads', type=int, default=1)
    parser.add_argument('--report')
    args = parser.parse_args(argv)
    stats = run(args.input, args.output, args.threads, args.report)
    print('\t'.join(f'{k}={v}' for k, v in stats.items()))
    return 0


if __name__ == '__main__':
    sys.exit(main())
