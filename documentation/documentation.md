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
We applied two different approaches, one using [cd-hit](http://weizhongli-lab.org/cd-hit/) and the other [swarm](https://github.com/torognes/swarm).  
The scripts are `cd-hit_workflow.bash` and `swarm_workflow.bash`. Both scripts consist of applying the clustering tool on the pre-processed fasta.
`swarm_workflow.bash` includes a dereprlication step with vsearch, where the abundance of each read is added as _[INT] at the end of the sequence id.




