# Global Pairwise Sequence Alignment Algorithm

The Global Pairwise Sequence Alignment algorithm, typically realized through the Needleman-Wunsch algorithm, is a cornerstone method in bioinformatics for aligning two biological sequences—like nucleotide or protein sequences. Its purpose is to discover regions of similarity that could be due to functional, structural, or evolutionary relationships between the sequences.

## How It Works
- **Initialization**: Begin by creating a scoring matrix with dimensions `(length(sequence1) + 1) x (length(sequence2) + 1)`. The first row and column are populated with gap penalties that often increase linearly.
- **Matrix Filling**: Fill the matrix using a scoring function that accounts for matches, mismatches, and gaps. Each cell's score is the maximum of either coming from the cell above (signifying an insertion), from the left (a deletion), or from the diagonal cell (a match or mismatch).
- **Trace-back**: The final step involves backtracking from the bottom-right cell to the top-left cell, piecing together the optimal alignment. This trace follows the path of maximum scores, reflecting the highest-scoring alignment under the given scoring scheme.

## Scoring Parameters
- **Match Score**: A positive value when two elements from the sequences align.
- **Mismatch Penalty**: A negative value when two elements do not match.
- **Gap Penalty**: A negative value assigned for every gap introduced in the alignment.

## Applications
The algorithm is key in several areas, including:
- **Comparative Genomics**: Identifying homologous genes or regulatory sequences.
- **Phylogenetics**: Inferring evolutionary relationships.
- **Protein Modeling**: Predicting the three-dimensional structure of proteins based on their sequence alignment.

Global alignment is ideal when the sequences in question are of similar length and composition, implying a comprehensive alignment across their entire span.

## Wavefront Parallelism
In order to parallelize the sequential algorithm, the dependencies between matrix cells need to be eliminated. This can be done via wavefront parallelism.

Wavefront parallelism is a parallel computing approach that processes matrix data diagonally, allowing for the execution of dependent computations concurrently. Each matrix point is calculated after its immediate predecessors—above, left, and diagonally left—ensuring data integrity. This method is efficient for algorithms where data points are interdependent, such as dynamic programming and numerical simulations.

## Implementation
The algorithm has been coded in C++, parallelization has been done via OpenMP. More specifically, each matrix cell of a wave (as specified in the wavefront parallelism description) is packaged in a task to be finished by consumer threads.

## Files

To run with a specific data set use: ./gpsa --x <sequence1_filename> --y <sequence2_filename>
By default, the program will look for X.txt and Y.txt. 

Here are the available sequences: 

1. X.txt, Y.txt, size: [18481x18961] 
Random, big sequences.

2. X2.txt, Y2.txt, size: [16383x16383] 
Same, but this has the same size dimensions that nicely divide.

3. simple1.txt, simple2.txt, size: [5x6] 
Small sequences that you can use for debugging. 

## Source

Cool article for the Needleman-Wunsch algorithm: https://medium.com/analytics-vidhya/sequence-alignment-and-the-needleman-wunsch-algorithm-710c7b1a23a4

Great for understanding wavefront parallelism: https://link.springer.com/article/10.1007/s10766-020-00658-y

Parts of the code and dataset used here are provided courtesy of the *University of Vienna*.
