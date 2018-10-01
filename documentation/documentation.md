# Documentation

## Pre-processing (original)
The pre-processing script is `preprocess_workflow.bash`.
It consists of the following tasks:
1. Merge pair-end reads with [pear](https://sco.h-its.org/exelixis/web/software/pear/doc.html) (default parameters).
2. Quality trim with [bbduk](https://sourceforge.net/projects/bbmap/). We trim bases with a quality lower than 25 (trimq=25) from both ends (qtrim=rl) in merged and unmerged reads, and discard reads shorter than 100bp (minlength=100).
3. Concatenate all files in single multifastq.
4. Convert fastq to fasta.
5. Dereplicate with [vsearch](https://github.com/torognes/vsearch) keeping singletons (--minuniquesize 1).
6. Check for chimeras with [vsearch](https://github.com/torognes/vsearch). We use a minimum abundance ratio of parent vs. chimera of 2 (--abskew  1.5) and the fasta is output as one line per sequence (--fasta_width 0).
7. Count sequence number and length.
8. Compute R1 and R2 stats with [vsearch](https://github.com/torognes/vsearch) --fastq_stats.
9. Clean intermediate files.


## Pre-processing (Mallorca pipeline)
The pre-processing script is `preprocess_mallorca_workflow.bash`.
It consists of the following tasks:
1. Quality trim with [SolexaQA](http://solexaqa.sourceforge.net/). We trim bases with a quality lower than 20 (--phredcutoff 20).
2. Rename files (SolexaQA output names are messy).
3. Remove short sequences with [SolexaQA](http://solexaqa.sourceforge.net/). Sequences shorter than 50bp are discarded (--length 50).
4. Rename files.
5. Trim adapters in pair-end reads with [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic), using the TruSeq3-PE.fa adapters (HiSeq and MiSeq machines).
6. Trim adapters in single-end reads with [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic), using the TruSeq3-SE.fa adapters (HiSeq and MiSeq machines).
7. Merge with Flash [Flash](https://ccb.jhu.edu/software/FLASH/).
8. Count sequence number and length.
9. Concatenate all preprocessed reads (single-end + merged reads).
10. Convert to fasta with [fq2fh.sh](https://github.com/pereiramemo/16S_analysis_pipelines/blob/master/scripts/fq2fa.sh).
11. Rename sequences: add sample name.
12. Clean intermediate files.

## Clustering
We applied two different approaches, one using [cd-hit](http://weizhongli-lab.org/cd-hit/) and the other [swarm](https://github.com/torognes/swarm).  
The scripts are `cd-hit_workflow.bash` and `swarm_workflow.bash`. Both scripts consist of applying the clustering tool on the pre-processed fasta.
`swarm_workflow.bash` includes a dereprlication step with vsearch, where the abundance of each read is added as "_[INT]" at the end of the sequence id.




