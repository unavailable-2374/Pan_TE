import random
import sys
import unittest
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
UTILS = ROOT / "bin" / "Refiner" / "utils"
sys.path.insert(0, str(UTILS))

from complexity_utils import calculate_dust_score  # noqa: E402


def legacy_dust(sequence, window=64, threshold=7):
    if not sequence or len(sequence) < window:
        return 0.0
    sequence = sequence.upper()
    scores = []
    for start in range(len(sequence) - window + 1):
        counts = Counter(
            sequence[index:index + 3]
            for index in range(start, start + window - 2)
            if "N" not in sequence[index:index + 3]
        )
        score = sum(count * (count - 1) / 2
                    for count in counts.values() if count > 1)
        scores.append(score)
    maximum = (window - 2) * (window - 3) / 2
    return min((sum(scores) / len(scores)) / maximum * threshold, 1.0)


class RollingDustTests(unittest.TestCase):
    def test_matches_legacy_algorithm(self):
        randomizer = random.Random(20260723)
        for length in (1, 20, 63, 64, 65, 100, 1000):
            for _ in range(20):
                sequence = "".join(
                    randomizer.choice("ACGTN") for _ in range(length))
                self.assertAlmostEqual(
                    calculate_dust_score(sequence), legacy_dust(sequence),
                    places=12)

    def test_low_complexity_extremes(self):
        self.assertEqual(calculate_dust_score("A" * 200), 1.0)
        self.assertLess(calculate_dust_score("ACGTGCACTGATCGTA" * 20), 0.7)


if __name__ == "__main__":
    unittest.main()
