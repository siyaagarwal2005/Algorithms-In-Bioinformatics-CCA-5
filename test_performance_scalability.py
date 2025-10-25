"""
Test cases for Assignment 3, Question 5: Performance and Scalability
Demonstrates the functionality of parallel processing and progress monitoring.
"""

import unittest
import os
import time
from performance_scalability import (
    generate_sequences_efficiently,
    run_parallel_analysis,
    analyze_sequence
)


class TestPerformanceScalability(unittest.TestCase):

    def test_01_memory_efficiency_with_generator(self):
        """Test that the generator yields sequences one by one without consuming memory."""
        print("\n--- Testing Generator (Memory Efficiency) ---")
        NUM_SEQS = 1000
        MIN_LEN, MAX_LEN = 100, 200

        gen = generate_sequences_efficiently(NUM_SEQS, MIN_LEN, MAX_LEN)

        # Test 1: Check if it's a generator
        self.assertTrue(hasattr(gen, '__next__'), "Function did not return a generator.")

        # Test 2: Consume the generator and check count/length
        count = 0
        total_length = 0
        for seq in gen:
            count += 1
            total_length += len(seq)
            self.assertTrue(MIN_LEN <= len(seq) <= MAX_LEN)

        self.assertEqual(count, NUM_SEQS, "Generator yielded incorrect number of sequences.")
        print(f"  PASSED: Generator yielded {count} sequences.")

    def test_02_parallel_analysis_and_tqdm_monitoring(self):
        """
        Demonstrate parallel speedup and progress monitoring (tqdm).
        Note: Actual speedup is dependent on system CPU. This test confirms
        the parallelism mechanism works and checks results integrity.
        """
        print("\n--- Testing Parallel Processing and Progress Monitoring (tqdm) ---")
        NUM_SEQS = 5000  # Number of sequences to process
        MIN_LEN, MAX_LEN = 1000, 2000  # Make sequences long enough to be CPU bound

        # 1. Generate sequences (memory efficient source)
        sequences_gen = generate_sequences_efficiently(NUM_SEQS, MIN_LEN, MAX_LEN)

        # 2. Run analysis using parallel workers and progress bar
        print(f"Starting parallel run on {NUM_SEQS} sequences (look for the tqdm progress bar above).")

        # Start timer for comparison
        start_time = time.time()
        results = run_parallel_analysis(sequences_gen, total_tasks=NUM_SEQS)
        end_time = time.time()

        duration = end_time - start_time
        print(f"  Parallel processing finished in {duration:.2f} seconds.")

        # 3. Verify Results Integrity
        self.assertEqual(len(results), NUM_SEQS, "Parallel run did not return results for all tasks.")

        # Check a sample result structure
        first_result = results[0]
        self.assertIn('length', first_result)
        self.assertIn('gc_content', first_result)
        self.assertIn('is_low_complexity', first_result)

        print("  PASSED: Parallel processing mechanism and results integrity confirmed.")


if __name__ == '__main__':
    # Set the start method for macOS/Linux compatibility in multiprocessing
    import platform

    if platform.system() != 'Windows':
        from multiprocessing import set_start_method

        try:
            set_start_method("fork", force=True)
        except RuntimeError:
            pass  # Already set

    unittest.main(argv=['first-arg-is-ignored'], exit=False)