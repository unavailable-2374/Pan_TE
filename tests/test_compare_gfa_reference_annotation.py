#!/usr/bin/env python3

import sys
import unittest
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))

from compare_gfa_reference_annotation import (  # noqa: E402
    Hit,
    feature_comparison,
    intersect_bp,
    merge_intervals,
    subtract_intervals,
)


class ReferenceComparisonTests(unittest.TestCase):
    def test_interval_set_operations_use_half_open_coordinates(self):
        left = merge_intervals([(0, 10), (10, 15), (20, 30)])
        right = merge_intervals([(5, 22)])
        self.assertEqual(left, [(0, 15), (20, 30)])
        self.assertEqual(intersect_bp(left, right), 12)
        self.assertEqual(subtract_intervals(left, right), [(0, 5), (22, 30)])

    def test_feature_scan_handles_nested_and_fragmented_hits(self):
        gfa = [
            Hit("Chr1", 0, 100, "A", "LTR", 10, "g1", "n1"),
            Hit("Chr1", 20, 30, "B", "DNA", 20, "g2", "n2"),
            Hit("Chr1", 100, 150, "A", "LTR", 10, "g3", "n3"),
        ]
        direct = [Hit("Chr1", 10, 140, "A", "LTR", 30, "d1")]
        row = feature_comparison(gfa, direct, {"Chr1": [100]})[0]
        self.assertEqual(row["reference_nodes_crossed"], 2)
        self.assertEqual(row["any_overlap_bp"], 130)
        self.assertEqual(row["same_family_bp"], 130)
        self.assertEqual(row["same_family_gfa_fragments"], 2)


if __name__ == "__main__":
    unittest.main()
