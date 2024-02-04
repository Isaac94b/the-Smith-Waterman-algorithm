#  the Smith-Waterman algorithm

## Overview

The program demonstrates a local sequence alignment using the Smith-Waterman algorithm, allowing for the identification of similar regions within two sequences. The use of the termcolor library provides visual cues for better understanding the output.

## Features

- **Scoring Matrix:** The program generates a scoring matrix based on the input sequences and scoring parameters.
- **Alignment Paths:** It traces back the optimal alignment paths and displays the scores.
- **Aligned Sequences:** The program provides the aligned sequences for both input sequences.

## Usage

1. **Input Values:**
    - Run the program and enter the following values:
        - Match value
        - Mismatch value
        - Gap value

2. **Input Sequences:**
    - The program prompts for two input sequences. Alternatively, default sequences are provided for testing purposes.

3. **Output:**
    - The program displays the scoring matrix, the best alignment score, and the aligned sequences.
