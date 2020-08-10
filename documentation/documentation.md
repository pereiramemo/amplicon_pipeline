# Comments
All the pipelines are based on the Divisive Amplicon Denoising Algorithm implementation [DADA2](https://benjjneb.github.io/dada2/index.html) tool. 

This tool is able to resolve Amplicon Sequence Variants (ASVs), differing by as little as one nucleotide. 
ASVs provide biologically meaningful and consistent sequence labels, independently of reference databases. Also, ASVs provide a higher resolution of the taxonomic composition, compared to standard Operational Taxonomic Unit Clustering methods (Callahan 2017).

This repository contains three pipelines:
1. dada2_pipeline.R
2. primers_removal.R
3. taxa_annot.R

together with their corresponding wrap code:

1. dada2_pipeline_runner.R
2. primers_removal_runner.R
3. taxa_annot_runner.R

This way, all the pipelines can be executed from the command line. 
The pipelines customized the original code from the [DADA2 tutorial](https://benjjneb.github.io/dada2/tutorial.html).

# Bibliography
Callahan Benjamin J et al. (2017). "Exact sequence variants should replace operational taxonomic units in marker-gene data analysis". In: The ISME Journal 11, pp. 2639–2643. DOI:[10.1038/ismej.2017.119](https://www.nature.com/articles/ismej2017119)
