"""
Implements advanced algorithms for pattern matching (Suffix Arrays),
repeat/palindrome finding, sequence assembly simulation, and basic phylogeny.
"""

from typing import List, Tuple, Dict
import numpy as np


# --- 1. Implement suffix array concepts for pattern matching ---

def build_suffix_array(text: str) -> List[int]:
    """
    Builds a Suffix Array for a given text.

    A Suffix Array is a sorted array of all suffixes of a string.

    Args:
        text (str): The input sequence (genome).

    Returns:
        List[int]: The suffix array (indices of the start of sorted suffixes).

    Performance:
        Time Complexity: O(L^2 log L) for simple Python sort,
                         or O(L log L) with specialized algorithms (like DC3).
                         We use the simple O(L^2 log L) for clarity/simplicity.
        Space Complexity: O(L) for storing the array.
    """
    suffixes = [(text[i:], i) for i in range(len(text))]
    # Sort by the suffix string (Python's default sort is efficient)
    suffixes.sort(key=lambda x: x[0])

    # Return only the indices
    suffix_array = [index for suffix, index in suffixes]
    return suffix_array


def search_suffix_array(text: str, suffix_array: List[int], pattern: str) -> List[int]:
    """
    Searches a pattern in the text using the pre-built Suffix Array.

    Uses binary search to efficiently find the range of suffixes that start
    with the pattern.

    Args:
        text (str): The original text string.
        suffix_array (List[int]): The suffix array of the text.
        pattern (str): The pattern to search for.

    Returns:
        List[int]: A list of start indices where the pattern is found.

    Performance:
        Time Complexity: O(P log L), where P is pattern length and L is text length.
        Space Complexity: O(1) beyond the input arrays.
    """
    L = len(text)
    P = len(pattern)
    start_index = -1
    end_index = -1

    # Binary search for the first occurrence (leftmost bound)
    low, high = 0, L - 1
    while low <= high:
        mid = (low + high) // 2
        suffix = text[suffix_array[mid]: suffix_array[mid] + P]
        if suffix >= pattern:
            start_index = mid
            high = mid - 1
        else:
            low = mid + 1

    # Binary search for the last occurrence (rightmost bound)
    low, high = 0, L - 1
    while low <= high:
        mid = (low + high) // 2
        suffix = text[suffix_array[mid]: suffix_array[mid] + P]
        if suffix <= pattern:
            end_index = mid
            low = mid + 1
        else:
            high = mid - 1

    # Collect results
    if start_index != -1 and end_index != -1 and start_index <= end_index:
        return [suffix_array[i] for i in range(start_index, end_index + 1)]
    return []


# --- 2. Develop algorithms for finding repeats and palindromes ---

def is_palindrome(s: str) -> bool:
    """Checks if a string is a palindrome (reads same forward and backward)."""
    return s == s[::-1]


def find_palindromes(seq: str, min_len: int = 4) -> List[Tuple[int, int]]:
    """
    Finds all plain palindromic sequences (e.g., ATTA) in a sequence.

    Args:
        seq (str): The input sequence.
        min_len (int): Minimum length of palindrome to report.

    Returns:
        List[Tuple[int, int]]: List of (start_index, length) of palindromes.

    Performance:
        Time Complexity: O(L^3), where L is sequence length (L^2 substrings,
                         L for checking each palindrome). Can be optimized.
        Space Complexity: O(R), where R is the number of palindromes found.
    """
    palindromes = []
    L = len(seq)
    for i in range(L):
        for j in range(i + min_len, L + 1):
            substring = seq[i:j]
            if is_palindrome(substring):
                palindromes.append((i, len(substring)))
    return palindromes


def find_simple_repeats(seq: str, repeat_len: int = 4, min_count: int = 2) -> List[Tuple[str, int, int]]:
    """
    Finds contiguous tandem repeats (e.g., ACACAC) with a fixed unit length.

    Args:
        seq (str): The input sequence.
        repeat_len (int): The length of the repeating unit (e.g., 2 for ACAC).
        min_count (int): Minimum number of times the unit must repeat.

    Returns:
        List[Tuple[str, int, int]]: List of (unit, start_index, total_length).

    Performance:
        Time Complexity: O(L) - linear scan after setting up the initial window.
        Space Complexity: O(1).
    """
    repeats = []
    L = len(seq)
    min_total_len = repeat_len * min_count

    i = 0
    while i <= L - min_total_len:
        unit = seq[i: i + repeat_len]
        current_len = repeat_len
        count = 1

        # Extend the repeat as long as the unit pattern is matched
        while i + current_len + repeat_len <= L and \
                seq[i + current_len: i + current_len + repeat_len] == unit:
            current_len += repeat_len
            count += 1

        if count >= min_count:
            repeats.append((unit, i, current_len))
            # Jump past the found repeat to search for the next one
            i += current_len
        else:
            i += 1

    return repeats


# --- 3. Create tools for sequence assembly simulation (Overlap-Layout-Consensus) ---

def find_overlap(read1: str, read2: str, min_overlap: int = 12) -> int:
    """
    Finds the maximum suffix of read1 that is a prefix of read2.
    Returns the overlap length, or 0 if overlap is too small.
    """
    max_len = min(len(read1), len(read2))
    for overlap in range(max_len, min_overlap - 1, -1):
        if read1.endswith(read2[:overlap]):
            return overlap
    return 0


def simulate_assembly(reads: List[str], min_overlap: int = 12) -> str:
    """
    Simulates a greedy Overlap-Layout-Consensus (OLC) assembly process.

    It iteratively finds the best overlap and merges the reads. This is a
    simplification of a real OLC assembler.

    Args:
        reads (List[str]): List of short sequences (reads).
        min_overlap (int): Minimum required overlap length for a merge.

    Returns:
        str: The assembled contig (simplification: only returns the single,
             longest assembled contig).

    Performance:
        Time Complexity: O(N^2 * L) in the worst case, where N is the number of
                         reads and L is the read length, due to checking all N^2
                         pairs in a loop.
        Space Complexity: O(C), where C is the length of the final contig.
    """
    if not reads:
        return ""

    current_contigs = list(reads)

    while True:
        best_overlap = 0
        best_pair = None

        # O(N^2) pairwise comparison loop
        for i in range(len(current_contigs)):
            for j in range(len(current_contigs)):
                if i == j:
                    continue

                # Check A -> B overlap
                overlap = find_overlap(current_contigs[i], current_contigs[j], min_overlap)
                if overlap > best_overlap:
                    best_overlap = overlap
                    best_pair = (i, j, overlap)

        if best_pair is None:
            # No sufficient overlaps found, stop assembly
            break

        # Merge the best pair
        i, j, overlap = best_pair
        merged_read = current_contigs[i] + current_contigs[j][overlap:]

        # Update contig list: remove the two merged reads and add the new contig
        contig_i = current_contigs.pop(i)

        # J index needs adjustment if i < j (since i was removed first)
        j_adjusted = j - 1 if i < j else j
        contig_j = current_contigs.pop(j_adjusted)

        current_contigs.append(merged_read)

    return max(current_contigs, key=len) if current_contigs else ""


# --- 4. Build phylogenetic relationship analysis tools (Distance Matrix and UPGMA) ---

# We reuse the simple_evolutionary_distance function from Question 3 here
def simple_evolutionary_distance(seq1: str, seq2: str, model: str = 'p_distance') -> float:
    """
    (Reused from Q3) Calculates P-distance between two sequences.
    """
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must be of equal length for distance calculation.")
    length = len(seq1)
    if length == 0: return 0.0
    mismatches = 0
    for base1, base2 in zip(seq1, seq2):
        if base1 != base2:
            mismatches += 1
    return mismatches / length


def build_distance_matrix(sequences: List[str]) -> Tuple[List[str], np.ndarray]:
    """
    Calculates a pairwise distance matrix (using P-distance) for a set of aligned sequences.

    Args:
        sequences (List[str]): Aligned sequences.

    Returns:
        Tuple[List[str], np.ndarray]: (Sequence headers/names, Distance matrix).

    Performance:
        Time Complexity: O(N^2 * L), where N is number of sequences and L is length.
        Space Complexity: O(N^2) for the matrix.
    """
    N = len(sequences)
    if N == 0:
        return [], np.array([])

    # Check for equal length (required for distance models)
    L = len(sequences[0])
    if not all(len(s) == L for s in sequences):
        raise ValueError("All sequences must be of equal length (aligned).")

    # Generate names/headers (simple index for this tool)
    names = [f"Seq{i + 1}" for i in range(N)]
    matrix = np.zeros((N, N))

    for i in range(N):
        for j in range(i + 1, N):
            dist = simple_evolutionary_distance(sequences[i], sequences[j], model='p_distance')
            matrix[i, j] = dist
            matrix[j, i] = dist  # Matrix is symmetric

    return names, matrix
