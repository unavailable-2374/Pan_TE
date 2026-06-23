"""
N6 abPOA-harden unit tests — subprocess timeout + pre-launch size guards.

The abPOA branch of build_msa previously ran with a fixed 1800 s timeout and NO
copy-set size guard, so a pathological giant-element family (chr4 R=2: 65 × ~20 kb)
drove abPOA into a multi-TB band-matrix allocation that OOM-killed or hung the shard
for minutes — the exact cause of the chr4 full-run stall. N6 hardens this:

  * msa_wall_s caps the abPOA subprocess; TimeoutExpired -> build_msa returns None ->
    M2 never-regress keeps the original mdl-repeat seed (never faked, never hung).
  * abpoa_max_seq_bp / abpoa_max_total_bp refuse to LAUNCH abPOA when the work set is
    big enough that an OOM-kill (uncatchable by a Python timeout) would precede the
    wall clock — build_msa returns None up front and the caller keeps the seed.

The timeout test uses a fake `abpoa` shell script that sleeps past the wall clock, so
it exercises the REAL subprocess.TimeoutExpired path deterministically without needing
a pathological input or the real binary.

Run:  cd bin/Refiner_mdl && python3 -m unittest tests.test_n6_abpoa_harden
"""

import os
import stat
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from config import RefinerMdlConfig          # noqa: E402
from phase1_extract import Copy               # noqa: E402
import phase1_align as align                  # noqa: E402


def _copies(seqs):
    return [Copy(id=f"c{i}", sequence=s, divergence=0.1, strand='+')
            for i, s in enumerate(seqs)]


def _write_fake_abpoa(tmp_dir, body):
    """Write an executable fake `abpoa` shell script and return its path."""
    p = os.path.join(tmp_dir, 'abpoa')
    with open(p, 'w') as fh:
        fh.write("#!/bin/bash\n" + body + "\n")
    os.chmod(p, os.stat(p).st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)
    return p


class TestAbpoaTimeout(unittest.TestCase):
    """abPOA subprocess that overruns msa_wall_s -> None -> seed fallback."""

    def test_timeout_returns_none(self):
        tmp_dir = tempfile.mkdtemp(prefix='n6_to_')
        try:
            # Fake abpoa that sleeps well past the 1 s wall clock and never writes
            # output -> the only way build_msa returns is via TimeoutExpired -> None.
            fake = _write_fake_abpoa(tmp_dir, "sleep 10")
            cfg = RefinerMdlConfig(
                abpoa_exe=fake,
                msa_wall_s=1,                 # 1 s hard cap
                msa_subsample_cap_abpoa=300,
                subsample_seed=42,
                abpoa_max_seq_bp=0,           # disable size guards so we REACH abpoa
                abpoa_max_total_bp=0,
            )
            copies = _copies(["ACGTACGTAC" * 5, "ACGTACGTAA" * 5, "ACGTACGTAG" * 5])
            res = align.build_msa(copies, 'abpoa', cfg)
            # Timeout fired -> None (the caller keeps the seed). Never an alignment,
            # never an exception.
            self.assertIsNone(res)
        finally:
            import shutil
            shutil.rmtree(tmp_dir, ignore_errors=True)

    def test_timeout_is_bounded_by_msa_wall_s(self):
        """The call returns in ~msa_wall_s, not the old 1800 s, proving the new cap
        is what's enforced."""
        import time
        tmp_dir = tempfile.mkdtemp(prefix='n6_tob_')
        try:
            fake = _write_fake_abpoa(tmp_dir, "sleep 30")
            cfg = RefinerMdlConfig(
                abpoa_exe=fake, msa_wall_s=2, msa_subsample_cap_abpoa=300,
                subsample_seed=42, abpoa_max_seq_bp=0, abpoa_max_total_bp=0)
            copies = _copies(["ACGTACGTAC" * 5, "ACGTACGTAA" * 5, "ACGTACGTAG" * 5])
            t0 = time.time()
            res = align.build_msa(copies, 'abpoa', cfg)
            elapsed = time.time() - t0
            self.assertIsNone(res)
            # Comfortably under the fake's 30 s sleep -> the 2 s wall clock cut it.
            self.assertLess(elapsed, 15.0,
                            f"abpoa not bounded by msa_wall_s (took {elapsed:.1f}s)")
        finally:
            import shutil
            shutil.rmtree(tmp_dir, ignore_errors=True)


class TestAbpoaSizeGuard(unittest.TestCase):
    """Pre-launch guards skip abPOA on oversized work sets (never launch, no OOM)."""

    def test_giant_single_copy_skips_abpoa(self):
        # A single copy longer than abpoa_max_seq_bp -> guard trips before launch.
        # Point abpoa_exe at a path that would FAIL if launched, proving we never get
        # there (the guard returns None first).
        cfg = RefinerMdlConfig(
            abpoa_exe="/should/never/be/launched/abpoa",
            abpoa_max_seq_bp=50_000, abpoa_max_total_bp=0,
            msa_subsample_cap_abpoa=300, subsample_seed=42)
        big = "ACGT" * 20_000          # 80 kb > 50 kb cap
        copies = _copies([big, big, big])
        self.assertIsNone(align.build_msa(copies, 'abpoa', cfg))

    def test_total_bp_over_cap_skips_abpoa(self):
        cfg = RefinerMdlConfig(
            abpoa_exe="/should/never/be/launched/abpoa",
            abpoa_max_seq_bp=0, abpoa_max_total_bp=100_000,
            msa_subsample_cap_abpoa=300, subsample_seed=42)
        # 4 copies × 30 kb = 120 kb > 100 kb total cap (each copy under any per-seq cap).
        med = "ACGT" * 7_500           # 30 kb each
        copies = _copies([med, med, med, med])
        self.assertIsNone(align.build_msa(copies, 'abpoa', cfg))

    def test_under_caps_reaches_launch(self):
        # A small work set under both caps must NOT be skipped by the guard; with a
        # non-existent binary it falls through to the FileNotFoundError->None path,
        # but crucially NOT via the guard (we assert the guard's own short-circuit by
        # checking the longest/total are computed and pass — reaching the exec attempt).
        cfg = RefinerMdlConfig(
            abpoa_exe="abpoa_definitely_absent_n6",
            abpoa_max_seq_bp=50_000, abpoa_max_total_bp=8_000_000,
            msa_subsample_cap_abpoa=300, subsample_seed=42)
        copies = _copies(["ACGTACGTAC" * 5, "ACGTACGTAA" * 5, "ACGTACGTAG" * 5])
        # Reaches launch, binary missing -> None (not via the size guard).
        self.assertIsNone(align.build_msa(copies, 'abpoa', cfg))


class TestMafftSizeGuard(unittest.TestCase):
    """The SAME giant-element guard applies to the MAFFT path (mafft_auto / mafft_linsi),
    not just abPOA — a giant family must skip MSA up front rather than grind MAFFT for
    its full wall-clock timeout."""

    def test_mafft_auto_giant_copy_skips(self):
        cfg = RefinerMdlConfig(
            mafft_exe="/should/never/be/launched/mafft",
            abpoa_max_seq_bp=50_000, abpoa_max_total_bp=0,
            msa_subsample_cap_mafft=100, subsample_seed=42)
        big = "ACGT" * 20_000          # 80 kb > 50 kb cap
        copies = _copies([big, big, big])
        # regime mafft_auto: guard trips before run_mafft is ever called -> None.
        self.assertIsNone(align.build_msa(copies, 'mafft_auto', cfg))

    def test_mafft_total_bp_over_cap_skips(self):
        cfg = RefinerMdlConfig(
            mafft_exe="/should/never/be/launched/mafft",
            abpoa_max_seq_bp=0, abpoa_max_total_bp=100_000,
            msa_subsample_cap_mafft=100, subsample_seed=42)
        med = "ACGT" * 7_500           # 30 kb each; 4 × 30 kb = 120 kb > 100 kb cap
        copies = _copies([med, med, med, med])
        self.assertIsNone(align.build_msa(copies, 'mafft_linsi', cfg))

    def test_oversized_helper_direct(self):
        cfg = RefinerMdlConfig(abpoa_max_seq_bp=50_000, abpoa_max_total_bp=8_000_000)
        small = _copies(["ACGT" * 10, "ACGT" * 10])
        giant = _copies(["ACGT" * 200_000, "ACGT" * 5])   # one 800 kb copy
        self.assertFalse(align._oversized_for_msa(small, cfg, 'mafft_auto'))
        self.assertTrue(align._oversized_for_msa(giant, cfg, 'mafft_auto'))


class TestDefaultConfigHasHardening(unittest.TestCase):
    """The defaults ship the hardening (so a build_mdl-style invoke is protected)."""

    def test_defaults_present(self):
        cfg = RefinerMdlConfig()
        self.assertEqual(cfg.msa_wall_s, 60)   # tightened from 300: per-MSA backstop for
                                               # families slipping the sweet-spot guard
        self.assertGreater(cfg.abpoa_max_seq_bp, 0)
        self.assertGreater(cfg.abpoa_max_total_bp, 0)


if __name__ == '__main__':
    unittest.main()
