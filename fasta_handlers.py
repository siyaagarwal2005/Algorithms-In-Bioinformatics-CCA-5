
""""
This module provides robust handlers for reading and writing FASTA files,
designed to be efficient for large genomic datasets.

Functions:
    parse_fasta(file_path): Parses a FASTA file, yielding sequences one by one.
    write_fasta(file_path, sequences, line_width, mode): Writes sequences
                                                         to a FASTA file.
"""

from typing import Iterator, Tuple, Iterable


def parse_fasta(file_path: str) -> Iterator[Tuple[str, str]]:
    """
   Parses a FASTA file efficiently, yielding one sequence at a time.

    This function is a generator, which allows it to handle large genomic
    files  without loading the entire file into memory. It correctly
    parses files with multiple sequences [cite: 24] and extracts the full
    header line [cite: 24] and the corresponding sequence.

    Args:
        file_path (str): The path to the FASTA file.

    Yields:
        Iterator[Tuple[str, str]]: An iterator where each item is a
        tuple of (header, sequence).
            - header (str): The full header line (without the '>').
            - sequence (str): The complete, concatenated sequence.

    Raises:
        FileNotFoundError: If the specified file_path does not exist.
        IOError: If an error occurs during file reading.

    Performance:
        Time Complexity: O(N), where N is the total number of characters
        in the file. Each line is read exactly once.
        Space Complexity: O(L), where L is the length of the *longest
        sequence* in the file. It does not store all sequences in memory.
    """
    header = None
    sequence_parts = []

    try:
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue  # Skip empty lines

                if line.startswith('>'):
                    if header is not None:
                        yield (header, "".join(sequence_parts))

                    header = line[1:]
                    sequence_parts = []
                else:
                    if header is not None:
                        sequence_parts.append(line)

            if header is not None:
                yield (header, "".join(sequence_parts))

    except FileNotFoundError:
        print(f"Error: File not found at {file_path}")
        raise
    except IOError as e:
        print(f"An error occurred while reading the file: {e}")
        raise


def write_fasta(file_path: str,
                sequences: Iterable[Tuple[str, str]],
                line_width: int = 70,
                mode: str = 'w') -> None:
    """
    Writes sequences to a FASTA file with proper formatting.

    This function takes an iterable of (header, sequence) tuples and writes
    them to the specified file. The sequence data is wrapped to the
    specified `line_width`.

    Args:
        file_path (str): The path to the output FASTA file.
        sequences (Iterable[Tuple[str, str]]): An iterable (e.g., a list or
            generator) of (header, sequence) tuples to write.
        line_width (int, optional): The maximum number of characters per
            sequence line. Defaults to 70.
        mode (str, optional): The file opening mode. 'w' to overwrite
            (default), 'a' to append.

    Raises:
        ValueError: If `line_width` is not a positive integer.
        IOError: If an error occurs during file writing.

    Performance:
        Time Complexity: O(N), where N is the total number of characters
        being written (sum of all header and sequence lengths).
        Space Complexity: O(W), where W is `line_width`. The function
        writes sequence chunks line by line without holding large
        amounts of data in memory.
    """
    if line_width <= 0:
        raise ValueError("Line width must be a positive integer.")

    try:
        with open(file_path, mode) as f:
            for header, sequence in sequences:
                f.write(f">{header}\n")

                for i in range(0, len(sequence), line_width):
                    f.write(f"{sequence[i:i + line_width]}\n")

    except IOError as e:
        print(f"Error writing to file {file_path}: {e}")
        raise

