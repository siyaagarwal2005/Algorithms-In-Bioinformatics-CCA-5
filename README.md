# Assignment 3: Genomic Databases and Advanced Applications (CCA5)

## Project Overview

This repository contains the complete implementation for **Assignment 3: Genomic Databases and Advanced Applications (CCA5)**.  
The project progressively builds a suite of bioinformatics tools in Python — covering everything from efficient FASTA file handling and algorithmic problem-solving (sequence alignment, pattern matching, etc.) to advanced performance optimization and professional-grade command-line interface (CLI) development.

The final deliverable is a set of **well-documented, tested, and high-performance Python modules**.

---

## Setup and Installation

### Prerequisites
- Python **3.8+**
- `pip` package manager

### Dependencies

This project uses standard scientific computing libraries along with the **`tqdm`** library for progress monitoring.

Repository Structure
The repository is modular — each Python file corresponds to a specific assignment question for maximum clarity and reusability.


| File / Folder                   | Assignment Question        | Description                                                                                            |
| :------------------------------ | :------------------------- | :----------------------------------------------------------------------------------------------------- |
| `cli_tool.py`                   | **Q6 (Integration)**       | Main executable file implementing the CLI, logging, and error handling.                                |
| `fasta_handlers.py`             | **Q1 (FASTA Processing)**  | Efficient, generator-based parser for large FASTA files.                                               |
| `multiple_sequence_analysis.py` | **Q3 (Multiple Sequence)** | Implements global alignment, LCS, P-distance, and consensus sequence generation.                       |
| `advanced_pattern_analysis.py`  | **Q4 (Pattern Analysis)**  | Algorithms for suffix array construction, palindrome/repeat detection, and greedy assembly simulation. |
| `performance_scalability.py`    | **Q5 (Performance)**       | Demonstrates memory-efficient generators and multiprocessing with `tqdm`.                              |
| `test_*.py` files               | **Q3–Q5 (Testing)**        | Unit test files validating functionality for all major algorithms.                                     |
| `test_sequences.fasta`          | **Data**                   | Sample FASTA input file for testing the CLI.                                                           |
| `genomic_tool.log`              | **Q6 (Logging)**           | Automatically generated runtime log file.                                                              |

Key Algorithms Implemented


| Question | Requirement         | Implementation                                      | Time Complexity                       |
| :------- | :------------------ | :-------------------------------------------------- | :------------------------------------ |
| **Q3**   | Sequence Alignment  | **Needleman–Wunsch Global Alignment**               | O(M × N)                              |
| **Q3**   | Common Subsequences | **Longest Common Subsequence (LCS)** via 2D DP      | O(N × L²) for N sequences of length L |
| **Q4**   | Pattern Matching    | **Suffix Array** + Binary Search for pattern lookup | O(L log L)                            |
| **Q4**   | Assembly Simulation | **Greedy OLC (Overlap–Layout–Consensus)**           | O(N² × L) (worst case)                |
| **Q5**   | Memory Optimization | **Python Generators** for streaming large data      | O(1) memory                           |
| **Q5**   | Parallel Execution  | **`multiprocessing.Pool`** + **`tqdm`**             | ≈ O(L / P) where P = # of cores       |

