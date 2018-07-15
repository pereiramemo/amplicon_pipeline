# Documentation

## Pre-processing
The pre-processing script is 'preprocess_workflow.bash'.
The workflow consists of the following tasks.
1) Merge pair-end reads with [pear](https://sco.h-its.org/exelixis/web/software/pear/doc.html) (default parameters).
2) Quality trim with [bbduk](https://sourceforge.net/projects/bbmap/). We trim bases of merged and unmerged reads with a quality lower than 25 (trimq=25) from both ends (qtrim=rl), and discard reads shorter than 100bp (minlength=100).
3) Concatenate all files in single multifastq
4) Convert fastq to fasta
5) Dereplicate with [vsearch](https://github.com/torognes/vsearch): we remove all unique reads (--minuniquesize 1).
6) Check for chimeras with [vsearch](https://github.com/torognes/vsearch). We use a min abundance ratio of parent vs chimera of 1.5 (--abskew  1.5) and the fasta is one line per sequence (--fasta_width 0)
7) Count sequences and clean intermediary files.




