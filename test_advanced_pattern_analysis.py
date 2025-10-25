"""
Tests the algorithms in 'advanced_pattern_analysis.py'.
"""

import unittest
import numpy as np
from advanced_pattern_analysis import (
    build_suffix_array,
    search_suffix_array,
    find_palindromes,
    find_simple_repeats,
    simulate_assembly,
    build_distance_matrix
)


class TestAdvancedPatternAnalysis(unittest.TestCase):

    def test_01_suffix_array_and_search(self):
        """Test suffix array construction and pattern search."""
        print("\n--- Testing Suffix Array ---")
        text = "GATTACCA"
        pattern = "CCA"

        # Suffixes:
        # 0: GATTACCA
        # 1: ATTACCA
        # 2: TTACCA
        # 3: TACCA
        # 4: ACCA
        # 5: CCA
        # 6: CA
        # 7: A

        # Sorted Suffixes:
        # 7: A
        # 4: ACCA
        # 1: ATTACCA
        # 6: CA
        # 5: CCA  <-- Index 5
        # 0: GATTACCA
        # 3: TACCA
        # 2: TTACCA

        expected_sa = [7, 4, 1, 6, 5, 0, 3, 2]
        sa = build_suffix_array(text)
        self.assertEqual(sa, expected_sa, "Suffix Array build is incorrect.")

        # Test pattern search
        results = search_suffix_array(text, sa, pattern)
        expected_results = [5]  # CCA starts at index 5
        self.assertEqual(sorted(results), sorted(expected_results), "Suffix Array search is incorrect.")

        results_multi = search_suffix_array("BANANAS", build_suffix_array("BANANAS"), "NA")
        self.assertEqual(sorted(results_multi), [2, 4], "Suffix Array search (multi-hit) is incorrect.")

        print(f"  Text: {text}, Pattern: {pattern}")
        print(f"  SA: {sa}")
        print(f"  Results: {results}")
        print("  PASSED: Suffix Array")

    def test_02_repeats_and_palindromes(self):
        """Test algorithms for finding repeats and palindromes."""
        print("\n--- Testing Repeats and Palindromes ---")
        seq = "ATGCATTATTATTACGC"

        # Palindromes (min_len=4):
        # 1. TATTAT (index 4, length 6)
        # 2. ATGCAT (index 0, length 6) - *NOT* a plain palindrome, must be exact match

        # Palindromes in ATGCATTATTATTACGC:
        # ATTA (3, 4)
        # TATTAT (4, 6)
        expected_pals = [(3, 4), (4, 6), (9, 4), (10, 3),
                         (11, 2)]  # The O(L^3) finds all, including len 2/3 but we filter on min_len=4

        # Filter expected results to min_len=4
        expected_pals = [(3, 4), (4, 6), (9, 4)]

        palindromes = find_palindromes(seq, min_len=4)
        self.assertEqual(sorted(palindromes), sorted(expected_pals), "Palindrome finding is incorrect.")
        print(f"  Palindromes: {palindromes}")

        # Repeats (unit_len=3, min_count=2): ATT repeats 3 times
        # Unit: ATT, Start: 4, Total Length: 9 (ATTATTATT)
        expected_repeats = [("ATT", 4, 9)]
        repeats = find_simple_repeats(seq, repeat_len=3, min_count=2)
        self.assertEqual(repeats, expected_repeats, "Simple repeats finding is incorrect.")
        print(f"  Simple Repeats: {repeats}")

        print("  PASSED: Repeats and Palindromes")

    def test_03_assembly_simulation(self):
        """Test sequence assembly simulation (Greedy OLC)."""
        print("\n--- Testing Assembly Simulation ---")
        reads = [
            "GATTACCA",
            "ACCAGCTA",
            "ATCGATT",
            "CCTAG"
        ]
        # Merges (min_overlap=4):
        # 1. GATTACCA + ACCAGCTA (overlap 4: ACCA). New: GATTACCAGCTA
        # 2. GATTACCAGCTA + ATCGATT (overlap 3: ATT). New: GATTACCAGCTA (no, wait, 3 is not >= 4)
        # 3. GATTACCA + ACCAGCTA -> GATTACCAGCTA
        # 4. GATTACCAGCTA + CCTAG (no overlap)

        # Using min_overlap=3 to ensure a chain:
        reads = [
            "ATGC",
            "GCAT",
            "CATG"
        ]
        # ATGC + GCAT (overlap 2: GC). New: ATGCAT
        # ATGCAT + CATG (overlap 3: CAT). New: ATGCATG
        # Expected Contig: ATGCATG

        contig = simulate_assembly(reads, min_overlap=2)  # Using 2 for clearer test
        self.assertEqual(contig, "ATGCATG", "Assembly simulation contig is incorrect.")

        print(f"  Reads: {reads}")
        print(f"  Contig: {contig}")
        print("  PASSED: Assembly Simulation")

    def test_04_phylogenetic_distance_matrix(self):
        """Test building the distance matrix."""
        print("\n--- Testing Phylogenetic Distance Matrix ---")
        aligned_sequences = [
            "ATGCATGC",  # Seq1
            "ATGCATTT",  # Seq2 (1 mismatch at pos 7, Dist=1/8=0.125)
            "TTGCATGC"  # Seq3 (2 mismatches at pos 1, 2. Dist=2/8=0.25)
        ]

        names, matrix = build_distance_matrix(aligned_sequences)

        self.assertEqual(names, ['Seq1', 'Seq2', 'Seq3'], "Sequence names/headers are incorrect.")
        self.assertEqual(matrix.shape, (3, 3), "Matrix shape is incorrect.")

        # Check Seq1 vs Seq2
        self.assertAlmostEqual(matrix[0, 1], 0.125)
        # Check Seq1 vs Seq3
        self.assertAlmostEqual(matrix[0, 2], 0.25)
        # Check Seq2 vs Seq3
        # Seq2 (ATGCATTT) vs Seq3 (TTGCATGC) -> 4 mismatches (A/T, T/T, G/G, C/C, A/A, T/T, T/G, T/C)
        # Mismatches: 1, 2, 7, 8 (A/T, T/G, T/C). Dist = 3/8 = 0.375
        self.assertAlmostEqual(matrix[1, 2], 0.375)

        print("  Distance Matrix:")
        print(matrix)
        print("  PASSED: Phylogenetic Distance Matrix")


if __name__ == '__main__':
    unittest.main(argv=['first-arg-is-ignored'], exit=False)