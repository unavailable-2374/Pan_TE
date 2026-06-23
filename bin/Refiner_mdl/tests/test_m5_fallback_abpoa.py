"""
M5 unit tests — bounded blastn fallback + abPOA large-family MSA path.

Covers the M5 deliverables in isolation (no genome / DB required for most):
  * recruit_by_blastn returns [] (never crashes / never fabricates) when no DB is set.
  * _blast_hits_to_instances coordinate + strand + pident-filter math.
  * _parse_fasta_with_gaps tolerant reader for abPOA FASTA-with-gaps output.
  * build_msa(regime='abpoa') equal-length validation + failure-to-None behavior,
    including a REAL abpoa run when the binary is present (skipped otherwise).
  * ensure_genome_blast_db idempotency when the DB already exists on disk.

Run:  cd bin/Refiner_mdl && python3 -m unittest tests.test_m5_fallback_abpoa
"""

import os
import sys
import tempfile
import unittest
from types import SimpleNamespace

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from config import RefinerMdlConfig                      # noqa: E402
from phase1_extract import Copy                          # noqa: E402
import phase1_align as align                             # noqa: E402
import phase1_fallback as fb                             # noqa: E402

ABPOA_EXE = "/home/shuoc/tool/abPOA/bin/abpoa"
_HAS_ABPOA = os.path.exists(ABPOA_EXE)


def _cfg(**over):
    c = RefinerMdlConfig(genome_file="/nonexistent/genome.fa")
    for k, v in over.items():
        setattr(c, k, v)
    return c


# ═══════════════════════════════════════════════════════════════════════
# recruit_by_blastn — safety: no DB -> [] (never crash, never fabricate)
# ═══════════════════════════════════════════════════════════════════════

class TestRecruitNoDB(unittest.TestCase):

    def test_empty_db_returns_empty_list(self):
        cfg = _cfg(genome_blast_db="")
        rec = {'id': 'mdl_R9', 'sequence': 'ACGT' * 50}
        self.assertEqual(fb.recruit_by_blastn(rec, cfg), [])

    def test_missing_db_attr_returns_empty_list(self):
        # config without genome_blast_db attribute at all -> getattr default "" -> [].
        cfg = SimpleNamespace(min_recruit_pident=70.0,
                              enable_divergent_blast_recruitment=False)
        rec = {'id': 'mdl_R9', 'sequence': 'ACGT' * 50}
        self.assertEqual(fb.recruit_by_blastn(rec, cfg), [])

    def test_empty_sequence_returns_empty_list(self):
        cfg = _cfg(genome_blast_db="/some/db")
        self.assertEqual(fb.recruit_by_blastn({'id': 'x', 'sequence': ''}, cfg), [])


# ═══════════════════════════════════════════════════════════════════════
# _blast_hits_to_instances — coordinate / strand / pident math
# ═══════════════════════════════════════════════════════════════════════

class TestBlastHitsToInstances(unittest.TestCase):

    def test_plus_strand_coords(self):
        # 1-based inclusive [101, 300] -> 0-based half-open [100, 300), '+'.
        out = fb._blast_hits_to_instances("chr1\t101\t300\t95.0\t200", 70.0)
        self.assertEqual(len(out), 1)
        i = out[0]
        self.assertEqual((i.chrom, i.start, i.end, i.strand), ('chr1', 100, 300, '+'))
        self.assertAlmostEqual(i.divergence, 0.05, places=6)

    def test_minus_strand_from_descending_coords(self):
        # sstart > send => minus strand; span normalized to [lo-1, hi).
        out = fb._blast_hits_to_instances("chr2\t800\t601\t88.0\t200", 70.0)
        self.assertEqual(len(out), 1)
        i = out[0]
        self.assertEqual((i.chrom, i.start, i.end, i.strand), ('chr2', 600, 800, '-'))

    def test_pident_filter_drops_low_identity(self):
        txt = "chrA\t1\t100\t95.0\t100\nchrB\t1\t100\t60.0\t100"
        keep = fb._blast_hits_to_instances(txt, 70.0)
        self.assertEqual([i.chrom for i in keep], ['chrA'])

    def test_divergent_floor_keeps_more(self):
        txt = "chrA\t1\t100\t95.0\t100\nchrB\t1\t100\t50.0\t100"
        keep = fb._blast_hits_to_instances(txt, 45.0)
        self.assertEqual({i.chrom for i in keep}, {'chrA', 'chrB'})

    def test_seqid_prefix_stripped(self):
        out = fb._blast_hits_to_instances("ref|NC_003070.9|\t1\t100\t99.0\t100", 70.0)
        self.assertEqual(out[0].chrom, 'NC_003070.9')

    def test_malformed_lines_skipped(self):
        txt = "garbage\nchrA\t1\t100\t99.0\t100\nchrB\tNaN\tx\t99\t1"
        out = fb._blast_hits_to_instances(txt, 70.0)
        self.assertEqual([i.chrom for i in out], ['chrA'])

    def test_multiple_hsps_same_subject_all_kept(self):
        # Many HSPs on ONE chromosome (the contiguous-genome case '-max_hsps 1' broke):
        # all are parsed into separate Instances, not collapsed to one.
        txt = "\n".join(f"chr1\t{i*1000+1}\t{i*1000+300}\t90.0\t300\t200"
                        for i in range(8))
        out = fb._blast_hits_to_instances(txt, 70.0)
        self.assertEqual(len(out), 8)
        self.assertEqual({i.chrom for i in out}, {'chr1'})

    def test_max_keep_caps_by_score(self):
        # 5 hits with ascending scores; max_keep=2 keeps the two highest-score ones.
        lines = [
            "chrA\t1\t100\t80.0\t100\t50",
            "chrB\t1\t100\t80.0\t100\t400",
            "chrC\t1\t100\t80.0\t100\t100",
            "chrD\t1\t100\t80.0\t100\t300",
            "chrE\t1\t100\t80.0\t100\t10",
        ]
        out = fb._blast_hits_to_instances("\n".join(lines), 70.0, max_keep=2)
        self.assertEqual({i.chrom for i in out}, {'chrB', 'chrD'})


# ═══════════════════════════════════════════════════════════════════════
# _parse_fasta_with_gaps — tolerant abPOA reader
# ═══════════════════════════════════════════════════════════════════════

class TestParseFastaWithGaps(unittest.TestCase):

    def _write(self, text):
        fd, path = tempfile.mkstemp(suffix='.msa')
        os.close(fd)
        with open(path, 'w') as fh:
            fh.write(text)
        return path

    def test_order_preserved_and_gaps_kept(self):
        path = self._write(">s0\nACGT----ACGT\n>s1\nACGTACGTACGT\n")
        try:
            out = align._parse_fasta_with_gaps(path)
        finally:
            os.remove(path)
        self.assertEqual([n for n, _ in out], ['s0', 's1'])
        self.assertEqual(out[0][1], 'ACGT----ACGT')
        self.assertEqual(len(out[0][1]), len(out[1][1]))

    def test_multiline_rows_concatenated(self):
        path = self._write(">s0\nACGT\n----\n>s1\nACGT\nACGT\n")
        try:
            out = align._parse_fasta_with_gaps(path)
        finally:
            os.remove(path)
        self.assertEqual(out[0][1], 'ACGT----')
        self.assertEqual(out[1][1], 'ACGTACGT')


# ═══════════════════════════════════════════════════════════════════════
# build_msa(regime='abpoa') — validation + failure-to-None
# ═══════════════════════════════════════════════════════════════════════

class TestBuildMsaAbpoa(unittest.TestCase):

    def _copies(self, seqs):
        return [Copy(id=f"c{i}", sequence=s, divergence=0.1, strand='+')
                for i, s in enumerate(seqs)]

    def test_missing_binary_returns_none(self):
        cfg = _cfg(abpoa_exe="/definitely/not/here/abpoa",
                   msa_subsample_cap_abpoa=300, subsample_seed=42)
        # os.path.sep in path + not exists -> falls back to 'abpoa' on PATH; if that
        # is also absent FileNotFoundError -> None. Either way: never a fabricated MSA.
        copies = self._copies(["ACGTACGTACGT"] * 4)
        # Force PATH lookup to fail too by pointing PATH away is overkill; instead use
        # a bogus bare name that is not on PATH.
        cfg.abpoa_exe = "abpoa_definitely_absent_xyz"
        self.assertIsNone(align.build_msa(copies, 'abpoa', cfg))

    def test_unequal_rows_rejected_via_parse(self):
        # Exercise the equal-length guard: monkeypatch the parser to return ALL the
        # expected names (s0..s3) but with unequal row lengths, simulating a corrupt
        # abPOA output. build_msa must return None — never silently truncate.
        if not _HAS_ABPOA:
            self.skipTest("abpoa binary absent; needs a real run to reach the parser")
        cfg = _cfg(abpoa_exe=ABPOA_EXE, msa_subsample_cap_abpoa=300, subsample_seed=42)
        orig = align._parse_fasta_with_gaps
        align._parse_fasta_with_gaps = lambda p: [
            ('s0', 'ACGTACGT'), ('s1', 'ACGTACGT'),
            ('s2', 'ACGT'), ('s3', 'ACGTACGT')]   # s2 shorter -> unequal lengths
        try:
            copies = self._copies(["ACGTACGTACGT", "ACGTACGTACGT",
                                   "ACGTACGTACGT", "ACGTACGTACGT"])
            self.assertIsNone(align.build_msa(copies, 'abpoa', cfg))
        finally:
            align._parse_fasta_with_gaps = orig

    @unittest.skipUnless(_HAS_ABPOA, "abpoa binary not installed")
    def test_real_abpoa_run_produces_equal_length_rows(self):
        cfg = _cfg(abpoa_exe=ABPOA_EXE, msa_subsample_cap_abpoa=300,
                   subsample_seed=42)
        seqs = [
            "ACGTACGTACGTACGTACGTACGT",
            "ACGTACGTAGGTACGTACGTACGT",
            "ACGTACGTACGTACGAAGGTACGTACGT",
            "ACGTACGTACGTACGTACGTACGT",
            "ACGTACGTACGTACGTACGTTCGT",
        ]
        res = align.build_msa(self._copies(seqs), 'abpoa', cfg)
        self.assertIsNotNone(res)
        rows, metas = res
        self.assertEqual(len(rows), len(seqs))
        self.assertEqual(len(metas), len(seqs))
        lens = {len(r) for r in rows}
        self.assertEqual(len(lens), 1, f"rows not equal length: {lens}")
        # Rows are uppercase, gaps normalized to '-'.
        for r in rows:
            self.assertTrue(set(r) <= set('ACGTN-'))


# ═══════════════════════════════════════════════════════════════════════
# ensure_genome_blast_db — idempotency
# ═══════════════════════════════════════════════════════════════════════

class TestEnsureGenomeDB(unittest.TestCase):

    def test_existing_db_path_returned_without_rebuild(self):
        tmp = tempfile.mkdtemp(prefix='m5_db_')
        try:
            db = os.path.join(tmp, 'genome')
            open(db + '.nsq', 'w').close()  # pretend a built DB exists
            cfg = _cfg(genome_blast_db=db)
            self.assertEqual(fb.ensure_genome_blast_db(cfg), db)
        finally:
            import shutil
            shutil.rmtree(tmp, ignore_errors=True)

    def test_db_on_disk_in_tempdir_is_discovered(self):
        tmp = tempfile.mkdtemp(prefix='m5_db_')
        try:
            db_dir = os.path.join(tmp, 'genome_blastdb')
            os.makedirs(db_dir)
            open(os.path.join(db_dir, 'genome.nhr'), 'w').close()
            cfg = _cfg(temp_dir=tmp, genome_blast_db="")
            got = fb.ensure_genome_blast_db(cfg)
            self.assertEqual(got, os.path.join(db_dir, 'genome'))
            self.assertEqual(cfg.genome_blast_db, got)
        finally:
            import shutil
            shutil.rmtree(tmp, ignore_errors=True)


if __name__ == '__main__':
    unittest.main(verbosity=2)
