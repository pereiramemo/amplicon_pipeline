library(tidyverse)

###############################################################################
### 1. Load table
###############################################################################

ABUND_TBL_long <- read_tsv(file = "/bioinf/home/epereira/workspace/16S_analyses/\
lagunas_16S_analysis2/taxa_annot/amp2compare/clust2abund2taxa.tsv.gz",
                           col_names = T)


###############################################################################
### 2. Remove singletons
###############################################################################

# note: the column names are cluster, sample, abund, seq_rep and taxonomy

ABUND_TBL_long_redu <- ABUND_TBL_long %>%
                       group_by(cluster) %>%
                       mutate(total_abund = sum(abund)) %>%
                       filter(total_abund > 1) %>%
                       select(cluster, sample, abund, taxonomy)
 
                       
###############################################################################
### 3. Convert to wide format
###############################################################################

ABUND_TBL_wide <- spread(data = ABUND_TBL_long_redu, key = cluster, value = abund, fill = 0)



