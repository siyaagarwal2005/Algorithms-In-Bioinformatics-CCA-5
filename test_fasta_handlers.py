"""
This script tests the functions in 'fasta_handlers.py'.
"""

import os
from fasta_handlers import parse_fasta, write_fasta


def run_tests():
    """
    Provides test cases demonstrating the functionality of the
    FASTA handlers.
    """
    print("--- Running FASTA Handler Test Cases ---")

    # Define test data and filenames
    test_fasta_content = """
>seq1 | Homo sapiens | gene=ABC | metadata
ATGCGTACGTACGATGCGTACGTACGATGCGTACGTACGATGCGTACGTACG
ATGCGTACGTACGATGCGTACGTACGATGCGTACGTACGATGCGTACGTACG
GATTACA
>seq2 | Mus musculus | gene=XYZ
CCGTACGTACGATGCGTACGTACGATGCGTACGTACGATGCGTACGTACG
CCGTACGTACGATGCGTACGTACGATGCGTACGTACGATGCGTACGTACG
GATTACAGATTACAGATTACA
>seq3 | Empty Sequence

>seq4 | Single line sequence
AAAAAAAAAAAAAAAAAAAA
"""
    test_filename = "test_input.fasta"
    output_filename = "test_output.fasta"

    # 1. Create a dummy test FASTA file
    try:
        with open(test_filename, 'w') as f:
            f.write(test_fasta_content.strip())
        print(f"\n[Test 1] Created dummy file: {test_filename}")

        # 2. Test: Parse FASTA files and Extract Headers
        print("\n[Test 2] Parsing FASTA file...")

        parsed_data = list(parse_fasta(test_filename))

        assert len(parsed_data) == 4
        print(f"  -> PASSED: Correctly parsed {len(parsed_data)} sequences.")

        assert parsed_data[0][0] == "seq1 | Homo sapiens | gene=ABC | metadata"
        assert parsed_data[0][1].endswith("GATTACA")
        assert len(parsed_data[0][1]) == 107
        print("  -> PASSED: Seq1 header and sequence data are correct.")

        assert parsed_data[2][0] == "seq3 | Empty Sequence"
        assert parsed_data[2][1] == ""
        print("  -> PASSED: Correctly handled empty sequence (seq3).")

        # 3. Test: Implement FASTA writing
        print(f"\n[Test 3] Writing data to {output_filename}...")

        write_fasta(output_filename, parsed_data, line_width=50)
        print(f"  -> PASSED: Data written to {output_filename}.")

        # 4. Verification: Read the file we just wrote
        print("\n[Test 4] Verifying written file content...")
        written_data = list(parse_fasta(output_filename))

        assert len(written_data) == len(parsed_data)
        assert written_data[0][0] == parsed_data[0][0]
        assert written_data[0][1] == parsed_data[0][1]
        assert written_data[1][1] == parsed_data[1][1]
        print("  -> PASSED: Written file content matches original data.")

        # Check line wrapping
        with open(output_filename, 'r') as f:
            lines = f.readlines()

        assert len(lines[1].strip()) == 50
        assert len(lines[2].strip()) == 50
        assert len(lines[3].strip()) == 7
        print("  -> PASSED: Correct line wrapping (50 chars) verified.")

        print("\n--- All Tests Passed Successfully ---")

    except AssertionError as e:
        print(f"\n--- TEST FAILED ---")
        print(f"Assertion Error: {e}")
    except Exception as e:
        print(f"\n--- TEST FAILED ---")
        print(f"An unexpected error occurred: {e}")

    finally:
        # Cleanup
        if os.path.exists(test_filename):
            os.remove(test_filename)
        if os.path.exists(output_filename):
            os.remove(output_filename)
        print(f"\n[Cleanup] Removed {test_filename} and {output_filename}.")


if __name__ == "__main__":
    run_tests()