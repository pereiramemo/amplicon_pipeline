# Documentation

## Pre-processing (original)
The pre-processing script is `preprocess_workflow.bash`.
It consists of the following tasks:
1. Merge pair-end reads with [pear](https://sco.h-its.org/exelixis/web/software/pear/doc.html) (default parameters).
2. Quality trim with [bbduk](https://sourceforge.net/projects/bbmap/). We trim bases with a quality lower than 25 (trimq=25) from both ends (qtrim=rl) in merged and unmerged reads, and discard reads shorter than 100bp (minlength=100).
3. Concatenate all files in single multifastq.
4. Convert fastq to fasta.
5. Dereplicate with [vsearch](https://github.com/torognes/vsearch) keeping singletons (--minuniquesize 1).
6. Check for chimeras with [vsearch](https://github.com/torognes/vsearch). We use a minimum abundance ratio of parent vs. chimera of 2 (--abskew 2) and the fasta is output as one line per sequence (--fasta_width 0).
7. Count sequence number and length.
8. Compute R1 and R2 stats with [vsearch](https://github.com/torognes/vsearch) --fastq_stats.
9. Clean intermediate files.


## Pre-processing (workflow 2)
The pre-processing script is `preprocess_workflow2.bash`.
It consists of the following tasks:
1. Quality check with [bbduk](https://sourceforge.net/projects/bbmap/). We trim bases with a quality lower than 20 (trimq=20) from both ends (qtrim=rl), and discard reads shorter than 50bp (minlength=50).
2. Trim adapters in pair-end reads with [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic), using the TruSeq3-PE.fa adapters (HiSeq and MiSeq machines). Here again the minimum length (after trimming) is set to 50 (MINLEN:50).
3. Trim adapters in single-end reads with [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic), using the TruSeq3-SE.fa adapters (HiSeq and MiSeq machines). Minimum length set to 50 (MINLEN:50).
4. Merge with [Flash](https://ccb.jhu.edu/software/FLASH/).
5. Concatenate all preprocessed reads (single-end + merged reads).
6. Convert fastq to fasta with [fq2fh.sh](https://github.com/pereiramemo/16S_analysis_pipelines/blob/master/scripts/fq2fa.sh).
7. Check for chimeras with [vsearch](https://github.com/torognes/vsearch). We use a minimum abundance ratio of parent vs. chimera of 2 (--abskew 2) and the fasta is output as one line per sequence (--fasta_width 0).
8. Count sequence number and length.
9. Rename sequences: add sample name.  
10. Clean intermediate files.


## Pre-processing (workflow 3)
The pre-processing script is `preprocess_workflow2.bash`.
It consists of the following tasks:
1. Trim adapters in pair-end reads with [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic), using the TruSeq3-PE.fa adapters (HiSeq and MiSeq machines). Here again the minimum length (after trimming) is set to 50 (MINLEN:50).
2. Merge pair-end reads with [pear](https://sco.h-its.org/exelixis/web/software/pear/doc.html) (default parameters).
3. Quality check with [bbduk](https://sourceforge.net/projects/bbmap/). We trim bases with a quality lower than 20 (trimq=20) from both ends (qtrim=rl), and discard reads shorter than 50bp (minlength=50).
4. Convert fastq to fasta with [fq2fh.sh](https://github.com/pereiramemo/16S_analysis_pipelines/blob/master/scripts/fq2fa.sh).
5. Check for chimeras with [vsearch](https://github.com/torognes/vsearch). We use a minimum abundance ratio of parent vs. chimera of 2 (--abskew 2) and the fasta is output as one line per sequence (--fasta_width 0).
6. Rename sequences: add sample name and sequence id.
7. Count sequence number and length.
8. Clean intermediate files.

<p align="center">
<img src="https://github.com/pereiramemo/16S_analysis_pipelines/blob/master/figures/preprocess_workflow3.jpg">
</p>

## Clustering
We applied two different approaches, one using [cd-hit](http://weizhongli-lab.org/cd-hit/) and the other [swarm](https://github.com/torognes/swarm).  
The scripts are `cd-hit_workflow.bash` and `swarm_workflow.bash`. Both scripts consist of applying the clustering tool on the pre-processed fasta.
`swarm_workflow.bash` includes a dereprlication step with vsearch, where the abundance of each read is added as "_[INT]" at the end of the sequence id.

## Cross tables (abund2taxa_commands.bash)
1. Remove first three rows from TAXA_ANNOT; Select the sequence header and taxonomic classification, simplify header, convert tab to space in the taxonomic classification field, remove empty spaces at the end of the line.
2. Prepare .clster file for taxa annot propagation. Here we create a table with all the sequences matched to their representative.
3. Propagate the taxonomic annot to the table created in 2, with left_joiner2.perl.
4. Simplify header from SWARM_ABUND and remove empty spaces at the end of the line.
5. Check numbers: unique headers in TAXA_ANNOT and SWARM_ABUND, and the number of shared headers.
6. Cross tables with left_joiner2.perl




 



