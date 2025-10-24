"""
Tests the algorithms in 'multiple_sequence_analysis.py'.
"""

from multiple_sequence_analysis import *
import unittest


class TestSequenceAnalysis(unittest.TestCase):

    def test_01_global_alignment(self):
        """Test the basic global alignment (Needleman-Wunsch)."""
        print("\n--- Testing Global Alignment ---")
        seq1 = "ATGCGT"
        seq2 = "AGCT"

        # Expected optimal alignment for A-T-G-C-G-T and A-G-C-T
        # A T G C G T
        # A - G C - T
        # Score: (1*2) + (1*-2) + (1*2) + (1*2) + (1*-2) + (1*2) = 4

        score, align1, align2 = basic_global_alignment(seq1, seq2, match=2, mismatch=-2, gap=-2)

        print(f"  Seq1: {seq1} -> {align1}")
        print(f"  Seq2: {seq2} -> {align2}")
        print(f"  Score: {score}")

        self.assertEqual(score, 4, "Alignment score is incorrect.")
        self.assertEqual(align1, "ATGCGT", "Aligned sequence 1 is incorrect.")
        self.assertEqual(align2, "A-GC-T", "Aligned sequence 2 is incorrect.")

        print("  PASSED: Global Alignment")

    def test_02_multiple_sequence_lcs(self):
        """Test the iterative Longest Common Subsequence (LCS)."""
        print("\n--- Testing Multiple Sequence LCS ---")
        sequences = [
            "GATTACA",
            "GACTAC",
            "GATTCGA"
        ]
        # LCS(GATTACA, GACTAC) = GATAC
        # LCS(GATAC, GATTCGA) = GATC
        expected_lcs = "GATC"

        lcs = multiple_sequence_lcs(sequences)
        print(f"  Sequences: {sequences}")
        print(f"  LCS: {lcs}")

        self.assertEqual(lcs, expected_lcs, "Multiple Sequence LCS is incorrect.")
        print("  PASSED: Multiple Sequence LCS")

    def test_03_evolutionary_distance(self):
        """Test Hamming and P-distance models."""
        print("\n--- Testing Evolutionary Distance ---")
        seq_a = "ATGCCCTA"
        seq_b = "ATGCTACA"  # 3 mismatches (C/T, C/A, T/C)
        length = len(seq_a)  # 8

        # Hamming distance (count of differences)
        hamming_dist = simple_evolutionary_distance(seq_a, seq_b, model='hamming')
        self.assertEqual(hamming_dist, 3.0, "Hamming distance is incorrect.")
        print(f"  Hamming Distance: {hamming_dist}")

        # P-distance (proportion of differences)
        p_dist = simple_evolutionary_distance(seq_a, seq_b, model='p_distance')
        self.assertAlmostEqual(p_dist, 3.0 / 8.0, 4, "P-distance is incorrect.")
        print(f"  P-Distance: {p_dist:.4f}")

        print("  PASSED: Evolutionary Distance")

    def test_04_consensus_sequence(self):
        """Test consensus sequence and profile generation."""
        print("\n--- Testing Consensus Sequence ---")
        aligned_sequences = [
            "ATCGG-CTA",
            "ATGGG-GTA",
            "ATTGG-CTA",
            "ATGGG-TTA",
            "ATCGG-CTA"
        ]
        # Pos 1 2 3 4 5 6 7 8 9
        # A   5 0 0 0 0 0 2 0 5
        # C   0 0 1 0 0 0 3 0 0
        # G   0 0 0 5 0 0 0 1 0
        # T   0 0 4 0 0 0 0 4 0
        # Cns A T T G G - C T A (Gaps are ignored in majority vote)

        consensus, profile = generate_consensus_sequence(aligned_sequences)

        expected_consensus = "ATGGG-CTA"  # T wins at pos 3, G at pos 5
        self.assertEqual(consensus, expected_consensus, "Consensus sequence is incorrect.")

        # Check profile matrix accuracy at a key position (Pos 3: T=4, C=1)
        self.assertEqual(profile['T'][2], 4)
        self.assertEqual(profile['C'][2], 1)

        print(f"  Consensus Sequence: {consensus}")
        print("  Profile Matrix:")
        print(format_profile(profile, consensus))
        print("  PASSED: Consensus Sequence")


if __name__ == '__main__':
    unittest.main(argv=['first-arg-is-ignored'], exit=False)