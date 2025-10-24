"""
This module manages a SQLite database for storing, searching, and validating
genomic sequence information from various file formats (FASTA, FASTQ basics).
"""

import sqlite3
from typing import List, Tuple, Generator, Dict, Any

# Define the structure of the database
DB_NAME = "genomic_sequences.db"
TABLE_NAME = "sequences"


class GenomicDatabase:
    """
    Manages the SQLite database for genomic sequences, providing methods for
    initialization, insertion, search, and data validation.
    """

    def __init__(self, db_name: str = DB_NAME):
        """Initializes the database connection and ensures the table exists."""
        self.conn = sqlite3.connect(db_name)
        self.cursor = self.conn.cursor()
        self._create_table()

    def __del__(self):
        """Closes the database connection when the object is deleted."""
        self.conn.close()

    def _create_table(self):
        """Creates the sequences table if it does not already exist."""
        self.cursor.execute(f"""
            CREATE TABLE IF NOT EXISTS {TABLE_NAME} (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                header TEXT NOT NULL,
                sequence TEXT NOT NULL,
                source_format TEXT NOT NULL,
                length INTEGER,
                gc_content REAL,
                is_valid INTEGER -- 1 for True, 0 for False
            )
        """)
        self.conn.commit()
        print(f"Database table '{TABLE_NAME}' ensured.")

    # --- Data Validation and QC Measures ---
    def _calculate_qc(self, sequence: str) -> Tuple[int, float, int]:
        """Calculates length, GC content, and checks sequence validity."""
        length = len(sequence)

        # Check for non-standard DNA bases (basic validation)
        valid_bases = {'A', 'T', 'G', 'C', 'N'}
        if not set(sequence.upper()).issubset(valid_bases):
            is_valid = 0
            gc_content = 0.0
        else:
            is_valid = 1
            gc_count = sequence.upper().count('G') + sequence.upper().count('C')
            # Avoid division by zero for empty sequences
            gc_content = (gc_count / length) * 100 if length > 0 else 0.0

        return length, gc_content, is_valid

    # --- Handle Different File Formats (FASTA/FASTQ Basics) ---
    def insert_sequence(self, header: str, sequence: str, source_format: str):
        """
        Inserts a single sequence record, calculating QC metrics automatically.
        """
        length, gc_content, is_valid = self._calculate_qc(sequence)

        try:
            self.cursor.execute(f"""
                INSERT INTO {TABLE_NAME} 
                (header, sequence, source_format, length, gc_content, is_valid)
                VALUES (?, ?, ?, ?, ?, ?)
            """, (header, sequence, source_format, length, gc_content, is_valid))
            self.conn.commit()
        except sqlite3.Error as e:
            print(f"Database error during insertion: {e}")

    def load_from_fasta_parser(self,
                               fasta_generator: Generator[Tuple[str, str], None, None]):
        """
        Loads sequences from a FASTA parser generator (e.g., from Q1 solution).
        """
        print("Starting batch insertion from FASTA...")
        count = 0
        for header, sequence in fasta_generator:
            self.insert_sequence(header, sequence, "FASTA")
            count += 1
        print(f"Finished loading {count} sequences from FASTA.")

    def load_from_fastq_basics(self, fastq_data: List[Tuple[str, str]]):
        """
        Loads sequence and header data derived from a FASTQ file.
        (NOTE: Quality scores are ignored for this basic implementation).
        Args: fastq_data is a list of (header, sequence) tuples.
        """
        print("Starting batch insertion from FASTQ basics...")
        count = 0
        for header, sequence in fastq_data:
            # FASTQ headers start with '@', FASTA with '>'
            self.insert_sequence(header.lstrip('@'), sequence, "FASTQ")
            count += 1
        print(f"Finished loading {count} sequences from FASTQ basics.")

    # --- Implement Search and Retrieval Functions ---
    def search_by_header(self, query: str) -> List[Dict[str, Any]]:
        """Searches for records where the header contains the query string."""
        self.cursor.execute(f"""
            SELECT header, sequence, length, gc_content, is_valid, source_format
            FROM {TABLE_NAME}
            WHERE header LIKE ?
        """, (f'%{query}%',))

        return self._fetch_results()

    def search_by_gc_range(self, min_gc: float, max_gc: float) -> List[Dict[str, Any]]:
        """Searches for records with GC content within the specified range."""
        self.cursor.execute(f"""
            SELECT header, sequence, length, gc_content, is_valid, source_format
            FROM {TABLE_NAME}
            WHERE gc_content BETWEEN ? AND ?
            ORDER BY gc_content DESC
        """, (min_gc, max_gc))

        return self._fetch_results()

    def get_invalid_sequences(self) -> List[Dict[str, Any]]:
        """Retrieves all sequences flagged as invalid during QC."""
        self.cursor.execute(f"""
            SELECT header, sequence, length, gc_content, is_valid, source_format
            FROM {TABLE_NAME}
            WHERE is_valid = 0
        """)
        return self._fetch_results()

    def _fetch_results(self) -> List[Dict[str, Any]]:
        """Internal helper to fetch all results as a list of dictionaries."""
        columns = [desc[0] for desc in self.cursor.description]
        results = [dict(zip(columns, row)) for row in self.cursor.fetchall()]
        return results


# Re-using the FASTA parser from Question 1 to load data
# NOTE: This assumes fasta_handlers.py is available or its functions are defined here.
# For simplicity, we define a basic mock-up parser here for testing purposes.
def mock_fasta_parser(file_path: str) -> Generator[Tuple[str, str], None, None]:
    """MOCK: Reads a dummy file to simulate the FASTA parser from Q1."""
    # In a real scenario, you would import parse_fasta from fasta_handlers
    data = [
        ("SeqA|Human|High_GC", "CGGCGGCCGG"),  # High GC (80%)
        ("SeqB|Mouse|Low_GC", "ATTATATTAA"),  # Low GC (0%)
        ("SeqC|Rat|Invalid_Base", "ATGCXA"),  # Invalid base 'X'
        ("SeqD|Yeast|Medium_GC", "ATGCATGCATGC"),  # 50% GC
    ]
    for header, sequence in data:
        yield header, sequence