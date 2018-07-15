# Documentation

## Pre-processing
The pre-processing script is `preprocess_workflow.bash`.
It consists of the following tasks:
1. Merge pair-end reads with [pear](https://sco.h-its.org/exelixis/web/software/pear/doc.html) (default parameters).
2. Quality trim with [bbduk](https://sourceforge.net/projects/bbmap/). We trim bases from both ends (qtrim=rl) in merged and unmerged reads with a quality lower than 25 (trimq=25), and discard reads shorter than 100bp (minlength=100).
3. Concatenate all files in single multifastq.
4. Convert fastq to fasta.
5. Dereplicate with [vsearch](https://github.com/torognes/vsearch): we remove all unique reads (--minuniquesize 1).
6. Check for chimeras with [vsearch](https://github.com/torognes/vsearch). We use a minimum abundance ratio of parent vs. chimera of 1.5 (--abskew  1.5) and the fasta is output as one line per sequence (--fasta_width 0).
7. Count sequences and clean intermediate files.


## Clustering
We applied two approaches [cd-hit(http://weizhongli-lab.org/cd-hit/) and [swarm](https://github.com/torognes/swarm).
The scripts for these are cd-hit_workflow.bash and swarm_workflow.bash.





