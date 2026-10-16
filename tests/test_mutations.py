import unittest
from pathlib import Path

import pandas as pd


class TestMutationRecovery(unittest.TestCase):

    def test_three_expected_mutations_recovered(self):
        results_file = Path("results/mismatches.csv")

        self.assertTrue(
            results_file.exists(),
            "Run python pairwise_alignment.py before running this test."
        )

        results = pd.read_csv(results_file)

        expected_mutations = ["L7K", "V29I", "N164I"]

        self.assertEqual(len(results), 3)
        self.assertEqual(
            results["Mutation"].tolist(),
            expected_mutations
        )


if __name__ == "__main__":
    unittest.main()
