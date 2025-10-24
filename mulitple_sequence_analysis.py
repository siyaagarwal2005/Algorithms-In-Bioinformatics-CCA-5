"""
Implements basic algorithms for sequence alignment, common subsequence,
evolutionary distance, and consensus sequence generation.
"""

from typing import List, Tuple, Dict
import numpy as np


# --- 1. Implement algorithms for sequence alignment basics (Needleman-Wunsch) ---

def basic_global_alignment(seq1: str, seq2: str, match: int = 2, mismatch: int = -1, gap: int = -2) -> Tuple[
    int, str, str]:
    """
    Performs a basic Global Alignment (Needleman-Wunsch) between two sequences.

    This function calculates the optimal global alignment score and the
    aligned sequences using a simplified scoring matrix.

    Args:
        seq1 (str): The first sequence.
        seq2 (str): The second sequence.
        match (int): Score for a matching character.
        mismatch (int): Score for a mismatching character.
        gap (int): Score for opening or extending a gap.

    Returns:
        Tuple[int, str, str]: (Optimal alignment score, Aligned Seq1, Aligned Seq2)

    Performance:
        Time Complexity: O(M*N), where M and N are the lengths of seq1 and seq2.
        Space Complexity: O(M*N) for the score matrix and traceback matrix.
    """
    M, N = len(seq1), len(seq2)

    # 1. Initialization
    # Score matrix
    score = np.zeros((M + 1, N + 1), dtype=int)
    # Traceback matrix: 0=Stop, 1=Diagonal (Match/Mismatch), 2=Up (Gap in seq1), 3=Left (Gap in seq2)
    traceback = np.zeros((M + 1, N + 1), dtype=int)

    # Initialize first row and column
    for i in range(M + 1):
        score[i, 0] = i * gap
        traceback[i, 0] = 2  # Came from UP
    for j in range(N + 1):
        score[0, j] = j * gap
        traceback[0, j] = 3  # Came from LEFT

    # 2. Fill Matrix
    for i in range(1, M + 1):
        for j in range(1, N + 1):
            s = match if seq1[i - 1] == seq2[j - 1] else mismatch

            # Possible scores:
            diag_score = score[i - 1, j - 1] + s  # Diagonal (Match/Mismatch)
            up_score = score[i - 1, j] + gap  # Gap in seq2
            left_score = score[i, j - 1] + gap  # Gap in seq1

            score[i, j] = max(diag_score, up_score, left_score)

            # Set traceback pointer (prefer diagonal for ties)
            if score[i, j] == diag_score:
                traceback[i, j] = 1
            elif score[i, j] == up_score:
                traceback[i, j] = 2
            else:  # Must be left_score (or tie with left)
                traceback[i, j] = 3

    # 3. Traceback
    aligned_seq1, aligned_seq2 = "", ""
    i, j = M, N

    while i > 0 or j > 0:
        if traceback[i, j] == 1:  # Diagonal
            aligned_seq1 = seq1[i - 1] + aligned_seq1
            aligned_seq2 = seq2[j - 1] + aligned_seq2
            i -= 1
            j -= 1
        elif traceback[i, j] == 2:  # Up (Gap in seq2)
            aligned_seq1 = seq1[i - 1] + aligned_seq1
            aligned_seq2 = '-' + aligned_seq2
            i -= 1
        elif traceback[i, j] == 3:  # Left (Gap in seq1)
            aligned_seq1 = '-' + aligned_seq1
            aligned_seq2 = seq2[j - 1] + aligned_seq2
            j -= 1
        else:  # Should not happen, matrix is complete
            break

    return score[M, N], aligned_seq1, aligned_seq2


# --- 2. Find common subsequences between multiple sequences (Longest Common Subsequence) ---

def longest_common_subsequence_2d(seq1: str, seq2: str) -> str:
    """
    Finds the Longest Common Subsequence (LCS) for two sequences using
    Dynamic Programming.

    Args:
        seq1 (str): The first sequence.
        seq2 (str): The second sequence.

    Returns:
        str: The longest common subsequence string.

    Performance:
        Time Complexity: O(M*N)
        Space Complexity: O(M*N) for the DP matrix.
    """
    M, N = len(seq1), len(seq2)
    # DP matrix to store lengths of LCS
    L = np.zeros((M + 1, N + 1), dtype=int)

    # Fill DP table
    for i in range(1, M + 1):
        for j in range(1, N + 1):
            if seq1[i - 1] == seq2[j - 1]:
                L[i, j] = L[i - 1, j - 1] + 1
            else:
                L[i, j] = max(L[i - 1, j], L[i, j - 1])

    # Traceback to construct the LCS string
    lcs = ""
    i, j = M, N
    while i > 0 and j > 0:
        if seq1[i - 1] == seq2[j - 1]:
            lcs = seq1[i - 1] + lcs
            i -= 1
            j -= 1
        elif L[i - 1, j] > L[i, j - 1]:
            i -= 1
        else:
            j -= 1

    return lcs


# Note: For multiple sequences (N > 2), this problem is NP-hard.
# A common simplification is to find the LCS iteratively: LCS(S1, S2, S3) = LCS(LCS(S1, S2), S3).
def multiple_sequence_lcs(sequences: List[str]) -> str:
    """
    Finds the Longest Common Subsequence for multiple sequences iteratively.
    """
    if not sequences:
        return ""
    if len(sequences) == 1:
        return sequences[0]

    current_lcs = sequences[0]
    for i in range(1, len(sequences)):
        current_lcs = longest_common_subsequence_2d(current_lcs, sequences[i])

    return current_lcs


# --- 3. Calculate evolutionary distances using simple models (P-distance/Hamming) ---

def simple_evolutionary_distance(seq1: str, seq2: str, model: str = 'p_distance') -> float:
    """
    Calculates the evolutionary distance between two sequences based on
    simple models (Hamming distance or P-distance).

    Args:
        seq1 (str): The first sequence.
        seq2 (str): The second sequence.
        model (str): The distance model to use ('hamming' or 'p_distance').

    Returns:
        float: The calculated evolutionary distance.

    Raises:
        ValueError: If sequences are not of equal length or the model is unknown.

    Performance:
        Time Complexity: O(L), where L is the sequence length (linear scan).
        Space Complexity: O(1).
    """
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must be of equal length for distance calculation.")

    length = len(seq1)
    if length == 0:
        return 0.0

    mismatches = 0
    for base1, base2 in zip(seq1, seq2):
        if base1 != base2:
            mismatches += 1

    if model == 'hamming':
        # Hamming distance is the count of mismatches
        return float(mismatches)
    elif model == 'p_distance':
        # P-distance (p-hat) is the proportion of different sites
        return mismatches / length
    else:
        raise ValueError("Unknown distance model. Use 'hamming' or 'p_distance'.")


# --- 4. Develop consensus sequence generation algorithms ---

def generate_consensus_sequence(sequences: List[str]) -> Tuple[str, Dict[str, List[int]]]:
    """
    Generates a consensus sequence and a profile matrix from a list of
    aligned sequences.

    The consensus base at each position is the one that appears most frequently.
    Ties are broken alphabetically (A, C, G, T).

    Args:
        sequences (List[str]): A list of aligned sequences (must be same length).

    Returns:
        Tuple[str, Dict[str, List[int]]]:
        (Consensus sequence string, Profile matrix as a dict of lists)

    Performance:
        Time Complexity: O(N*L), where N is the number of sequences and L is
        the length of the sequences.
        Space Complexity: O(L) for the profile matrix.
    """
    if not sequences:
        return "", {}

    L = len(sequences[0])
    N = len(sequences)

    # Simple check for alignment
    if not all(len(s) == L for s in sequences):
        raise ValueError("All sequences must be of equal length (i.e., aligned).")

    bases = ['A', 'C', 'G', 'T']
    profile: Dict[str, List[int]] = {base: [0] * L for base in bases}
    consensus = []

    for i in range(L):  # Iterate over each column (position)
        column_counts = {base: 0 for base in bases}

        for seq in sequences:
            base = seq[i].upper()
            if base in column_counts:
                column_counts[base] += 1

        # Update profile matrix
        for base in bases:
            profile[base][i] = column_counts[base]

        # Determine consensus base for the column
        max_count = 0
        consensus_base = 'N'  # Default for positions with no valid bases

        # Iterate through bases in priority order (A, C, G, T) for tie-breaking
        for base in bases:
            count = column_counts.get(base, 0)
            if count > max_count:
                max_count = count
                consensus_base = base
            # Implicit tie-breaking: if count == max_count, the current base
            # (which is earlier alphabetically) is preferred.

        consensus.append(consensus_base)

    return "".join(consensus), profile


def format_profile(profile: Dict[str, List[int]], consensus_seq: str) -> str:
    """Formats the profile matrix and consensus sequence into a readable string."""
    header = "Pos: " + " ".join(f"{i + 1:2}" for i in range(len(consensus_seq)))
    profile_str = "\n" + header + "\n" + "=" * (len(header) + 5) + "\n"
    profile_str += f"Cns: " + "  ".join(consensus_seq) + "\n"

    for base in ['A', 'C', 'G', 'T']:
        counts = " ".join(f"{count:2}" for count in profile[base])
        profile_str += f"{base}:   {counts}\n"
    return profile_str