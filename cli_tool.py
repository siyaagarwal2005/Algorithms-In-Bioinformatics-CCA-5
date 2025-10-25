"""
Creates a professional-quality CLI tool that uses argparse, logging,
and comprehensive error handling.
"""

import argparse
import logging
import sys
import os
from typing import List, Tuple

# --- 1. Package Your Tools: Import functionality from a reusable module ---
try:
    from fasta_handlers import parse_fasta
except ImportError:
    # Fallback/error handling if the module is missing
    print("FATAL: Could not import 'fasta_handlers'. Ensure fasta_handlers.py is present.")
    sys.exit(1)

# --- 2. Implement Comprehensive Error Handling and Logging ---

# Configure logging to write to a file and the console
LOG_FILE = 'genomic_tool.log'
logging.basicConfig(
    level=logging.INFO,  # Default log level
    format='%(asctime)s - %(levelname)s - %(module)s - %(message)s',
    handlers=[
        logging.FileHandler(LOG_FILE),  # Log output to a file
        logging.StreamHandler(sys.stdout)  # Log output to the console
    ]
)
logger = logging.getLogger(__name__)


def process_fasta_cli(fasta_file: str, output_file: str, min_len: int):
    """
    The main processing logic function for the CLI.
    """

    logger.info(f"CLI tool started. Input: {fasta_file}, Output: {output_file}")
    logger.info(f"Applying filter: Minimum sequence length = {min_len} bp")

    processed_count = 0
    filtered_count = 0
    results: List[Tuple[str, int]] = []

    # Critical Try/Except Block for robust file operations
    try:
        # Use the reusable generator function
        sequence_generator = parse_fasta(fasta_file)

        # Process sequences one by one (memory efficient)
        for header, sequence in sequence_generator:
            processed_count += 1
            seq_len = len(sequence)

            if seq_len >= min_len:
                results.append((header, seq_len))
            else:
                filtered_count += 1

        # Log successful read operation
        logger.info(f"Read operation complete. Total sequences found: {processed_count}")

        # Write output
        with open(output_file, 'w') as out_f:
            out_f.write(f"--- Genomic Analysis Summary ---\n")
            out_f.write(f"Source File: {os.path.basename(fasta_file)}\n")
            out_f.write(f"Total Sequences Read: {processed_count}\n")
            out_f.write(f"Sequences Kept (>= {min_len} bp): {len(results)}\n")
            out_f.write(f"Sequences Filtered Out: {filtered_count}\n\n")

            for header, length in results:
                out_f.write(f"{length}\t{header}\n")

        logger.info(f"Results successfully written to: {output_file}")

    except (FileNotFoundError, IOError) as e:
        # Specific error handling for file I/O issues
        logger.error(f"FATAL I/O ERROR: Program terminated. Details: {e}")
        sys.exit(1)
    except Exception as e:
        # General catch-all for unexpected errors
        logger.critical(f"A CRITICAL, UNEXPECTED ERROR occurred. Program terminated. Details: {e}")
        sys.exit(1)

    logger.info("CLI tool finished successfully.")


# --- 3. Create Command-Line Interface (CLI) ---

def main():
    """
    Defines and runs the Command-Line Interface using argparse.
    """
    parser = argparse.ArgumentParser(
        description="A professional-grade CLI tool for genomic sequence filtering and summary reporting.",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Example Usage:\n  python cli_tool.py -i sequences.fasta -o filtered_summary.txt -l 500"
    )

    # Required argument: Input FASTA file
    parser.add_argument(
        '-i', '--input',
        required=True,
        help='Path to the required input FASTA file.'
    )

    # Optional arguments
    parser.add_argument(
        '-o', '--output',
        default='summary_output.txt',
        help='Path for the output summary file (default: summary_output.txt).'
    )
    parser.add_argument(
        '-l', '--min-len',
        type=int,
        default=100,
        help='Minimum sequence length (in bp) to keep (default: 100).'
    )

    args = parser.parse_args()

    # Call the core processing function
    process_fasta_cli(args.input, args.output, args.min_len)


if __name__ == "__main__":
    main()