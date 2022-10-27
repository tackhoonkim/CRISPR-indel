# CRISPR-indel
This is a C++ script used for aligning NGS reads against sgRNA library and detecting indel mutations.

The NGS data should be provided as paired end sequencing data with read 1 (forward read) containing reverse complement of the sgRNA target sequence (with potential indel mutations by the sgRNA paired with it), and read 2 (reverse read) containing the sgRNA. 

![NGS scheme](https://user-images.githubusercontent.com/62285278/198166170-58c29a94-639f-47ee-b89e-37e976770da9.png)

The script first identifies the sgRNA from the reverse read, then aligns the target sequence with the sgRNA sequence. If there is a mismatch, it was counted as an indel.


Input data format:
T1-1.fastq: read 1 (forward read) containing the sgRNA target sequence.
T1-2.fastq: read 2 (reverse read) containing the sgRNA sequence.

Output data:
statistics.txt:
Number of coupling_matches: Number of reads whose sgRNA and target sequences are paired.
Number of decoupling_matches: Number of reads whose sgRNA and target sequences are not paired. These reads were excluded from analysis.
Number of indel: Number of reads with indel mutation detected in read1.
Number of no_indel: Number of reads with no detected indel mutation.
Number of non perfect matches: Number of reads with sgRNA sequence not matching with any of the sgRNAs in the library.
Coupling ratios: The fraction of reads with sgRNA and target sequences paired: coupling_matches/(coupling_matches+decoupling_matches)

library_count.txt:
sequence: the sgRNA sequence
target count: The number of reads with the corresponding sgRNA sequence.
indel: The number of reads with corresponding sgRNA sequence and detected indel in target sequence.
no indel: The number of reads with corresponding sgRNA sequence and and no detected indel.
