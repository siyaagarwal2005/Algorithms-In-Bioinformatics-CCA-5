"""
This script tests the functionality of the 'GenomicDatabase' class in genomic_db.py.
"""

import os
import unittest
import sqlite3
from genomic_db import GenomicDatabase, mock_fasta_parser, DB_NAME, TABLE_NAME


class TestGenomicDB(unittest.TestCase):
    """
    Test suite for the GenomicDatabase class.
    We use a temporary in-memory database or clean up the file database.
    """

    # Use the file-based DB for easier inspection and cleanup
    TEST_DB_NAME = "test_genomic_sequences.db"

    def setUp(self):
        """Setup: Initialize a fresh database before each test."""
        if os.path.exists(self.TEST_DB_NAME):
            os.remove(self.TEST_DB_NAME)
        self.db = GenomicDatabase(db_name=self.TEST_DB_NAME)
        print(f"\nSetting up database: {self.TEST_DB_NAME}")

    def tearDown(self):
        """Teardown: Close connection and remove the test database file."""
        del self.db  # Ensures __del__ is called to close connection
        if os.path.exists(self.TEST_DB_NAME):
            os.remove(self.TEST_DB_NAME)
        print(f"Tearing down and removing database: {self.TEST_DB_NAME}")

    def test_01_db_initialization_and_table_creation(self):
        """Test the database connection and table creation."""
        conn = sqlite3.connect(self.TEST_DB_NAME)
        cursor = conn.cursor()

        # Check if the table exists
        cursor.execute(f"SELECT name FROM sqlite_master WHERE type='table' AND name='{TABLE_NAME}';")
        self.assertIsNotNone(cursor.fetchone(), "Table was not created successfully.")

        conn.close()

    def test_02_qc_and_fasta_loading(self):
        """Test loading from FASTA (via generator) and QC calculations."""
        # Simulate loading from a FASTA file using the mock parser
        fasta_gen = mock_fasta_parser("dummy_path.fasta")
        self.db.load_from_fasta_parser(fasta_gen)

        # Verify total count
        count = self.db.cursor.execute(f"SELECT COUNT(*) FROM {TABLE_NAME}").fetchone()[0]
        self.assertEqual(count, 4, "Incorrect number of sequences loaded.")

        # Retrieve the High GC sequence for QC check
        high_gc_seq = self.db.search_by_header("High_GC")[0]
        self.assertEqual(high_gc_seq['length'], 10)
        # GC content calculation: 8 Gs, 2 Cs = 10 GC / 10 total = 100% (Error in mock data)
        # Recalculate based on mock data: CGGCGGCCGG -> 8 G/C out of 10. GC%=80.0
        self.assertAlmostEqual(high_gc_seq['gc_content'], 80.0)
        self.assertEqual(high_gc_seq['is_valid'], 1)

    def test_03_fastq_loading_and_source_format(self):
        """Test loading from FASTQ basics and checking the source tag."""
        fastq_data = [
            ("@Read1_A", "ATGCATGC"),
            ("@Read2_B", "GATTACA"),
        ]
        self.db.load_from_fastq_basics(fastq_data)

        # Verify total count
        count = self.db.cursor.execute(f"SELECT COUNT(*) FROM {TABLE_NAME}").fetchone()[0]
        self.assertEqual(count, 2, "Incorrect number of sequences loaded from FASTQ.")

        # Verify source format and header cleanup
        read1 = self.db.search_by_header("Read1_A")[0]
        self.assertEqual(read1['source_format'], "FASTQ")
        self.assertEqual(read1['header'], "Read1_A")  # Header should be clean

    def test_04_search_and_retrieval(self):
        """Test search functions (header and QC range)."""
        self.db.load_from_fasta_parser(mock_fasta_parser("dummy"))  # Load test data

        # Test search by header (partial match)
        mouse_seqs = self.db.search_by_header("Mouse")
        self.assertEqual(len(mouse_seqs), 1)
        self.assertTrue("Low_GC" in mouse_seqs[0]['header'])

        # Test search by GC range (0.0 to 1.0 - should return only SeqB (0%))
        low_gc_seqs = self.db.search_by_gc_range(0.0, 1.0)
        self.assertEqual(len(low_gc_seqs), 1)
        self.assertTrue("Low_GC" in low_gc_seqs[0]['header'])

        # Test search by GC range (49.0 to 51.0 - should return SeqD (50%))
        medium_gc_seqs = self.db.search_by_gc_range(49.0, 51.0)
        self.assertEqual(len(medium_gc_seqs), 1)
        self.assertTrue("Medium_GC" in medium_gc_seqs[0]['header'])

    def test_05_data_validation_retrieval(self):
        """Test the ability to find invalid sequences."""
        self.db.load_from_fasta_parser(mock_fasta_parser("dummy"))  # Load test data

        invalid_seqs = self.db.get_invalid_sequences()
        self.assertEqual(len(invalid_seqs), 1)
        self.assertTrue("Invalid_Base" in invalid_seqs[0]['header'])
        self.assertEqual(invalid_seqs[0]['is_valid'], 0)
        self.assertEqual(invalid_seqs[0]['sequence'], "ATGCXA")


if __name__ == '__main__':
    unittest.main(argv=['first-arg-is-ignored'], exit=False)