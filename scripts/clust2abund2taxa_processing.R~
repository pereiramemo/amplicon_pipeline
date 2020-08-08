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
                       group_by(cluster) %>% # group table by cluster, this allows the application of functions within each cluster
                       mutate(total_abund = sum(abund)) %>% # create a new variable with the total abund per cluster
                       filter(total_abund > 1) %>% # select cluster with total_abund greater than x. In this case x = 1.
                       select(cluster, sample, abund , taxonomy) # select relevant columns

###############################################################################
### 3. Convert to wide format: otu x sample table with taxa annot
###############################################################################

ABUND_TBL_wide <- spread(data = ABUND_TBL_long_redu, key = sample, value = abund, fill = 0) 

###############################################################################
### 4. Format to abudnance matrix: sample x otu
###############################################################################

ABUND_TBL_wide_trans <- ABUND_TBL_wide %>% 
                        select(-taxonomy) %>% 
                        remove_rownames() %>%
                        column_to_rownames("cluster") %>%
                        t() %>%
                        as.data.frame()



