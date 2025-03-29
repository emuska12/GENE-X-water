# GENE-X-water
Lokálne zarovnávanie sekvencií (Smith-Waterman) pomocou DP a rekurzie v Streamlit app.

# Local Sequence Alignment Tool

This repository contains an implementation of the **Smith-Waterman algorithm** for local sequence alignment of biological sequences, specifically DNA (nucleotide) and protein (amino acid) sequences. The tool is designed as an interactive web application built with **Streamlit**, allowing users to align sequences efficiently using two methods: **dynamic programming** and **recursion with memoization**. It was developed as part of the bioinformatics course.

## Project Motivation
In this project, we imagine discovering biological samples from a meteorite, containing unknown DNA and protein sequences. After initial analysis, these sequences appear similar to those found on Earth. To investigate potential similarities and explore whether this could indicate a new form of life, we need a robust tool to compare these sequences against known Earth sequences. This application implements local sequence alignment to identify regions of maximum similarity, leveraging the Smith-Waterman algorithm for its precision in finding local matches.

## Features
- **Sequence Types**: Supports alignment of DNA sequences (using an EDNA-like scoring matrix) and protein sequences (using BLOSUM62 or PAM250 substitution matrices).
- **Alignment Methods**: Offers two approaches:
  - **Dynamic Programming**: Iterative computation for efficient alignment.
  - **Recursion with Memoization**: Recursive computation with caching for optimized performance.
- **Customizable Scoring**: Users can adjust scoring parameters:
  - Match score (default: +1)
  - Mismatch penalty (default: -1)
  - Gap opening penalty (default: -10)
  - Gap extension penalty (default: -0.5)
- **Input Options**: Upload FASTA files or manually enter sequences via text input.
- **Output**: Detailed alignment visualization, including:
  - Aligned sequences with annotations (e.g., `|` for matches, `:` for similar amino acids, `.` for mismatches).
  - Statistics such as identity, similarity, gap percentage, and score.
  - Execution time for performance comparison between methods.

## Technical Details
The Smith-Waterman algorithm is a local alignment method that differs from global alignment by setting negative scores to zero, focusing only on regions of high similarity. The implementation includes:
1. **Initialization**: Creates scoring and traceback matrices, with special handling for affine gap penalties.
2. **Matrix Filling**: Computes optimal scores using either dynamic programming or recursive memoization.
3. **Traceback**: Reconstructs the aligned sequences from the maximum score position back to a zero score.

The application uses:
- **NumPy** for efficient matrix operations.
- **BioPython** for parsing FASTA files and loading substitution matrices (BLOSUM62, PAM250).
- **Streamlit** for an intuitive user interface.
