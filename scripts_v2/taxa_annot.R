#!/usr/bin/env Rscript

# This code is a modified version of the DADA2 tutorial: 
# https://benjjneb.github.io/dada2/tutorial.html

###############################################################################
### 1. Def. env
###############################################################################

suppressMessages(suppressWarnings(library(dada2)))
suppressMessages(suppressWarnings(library(tidyverse)))

REPO <- "/home/epereira/workspace/repositories/amplicon_pipelines/"
toolbox <- file.path(REPO, "scripts_v2/toolbox.R")
source(toolbox)

###############################################################################
### 2. Get parameters
###############################################################################

args = commandArgs(trailingOnly=TRUE)

INPUT_ASV_TABLE <- args[1]
INPUT_FASTA <- args[2]
OUTPUT_ASV_TABLE <- args[3]
METHOD <- args[4]
TRAIN_DB <- args[5]
REF_DB <- args[6]
BLAST_DB <- args[7]
TAXA_MAP <- args[8]
SEQ_MAP <- args[9]
EVALUE <- args[10]
MIN_IDENTITY <- args[11]
BLOUT <- args[12]
NSLOTS <- args[13] %>% as.numeric
SAVE_WORKSPACE <- args[14] %>% as.logical()

# INPUT_FASTA <- "~/workspace/indicadores_cuencas_2018/output/asv_run4/asvs_redu.fasta"
# INPUT_ASV_TABLE <- "~/workspace/indicadores_cuencas_2018/output/asv_run4/asv_table_redu.csv"
# OUTPUT_FASTA <- "~/workspace/indicadores_cuencas_2018/output/asv_run4/asvs_taxa_annot.fasta"
# OUTPUT_ASV_TABLE <- "~/workspace/indicadores_cuencas_2018/output/asv_run4/asv_table_taxa_annot.csv"
# REF_DB <- "/home/bioinf/resources/silva/silva_v138/silva_species_assignment_v138.fa.gz"
# TRAIN_DB <- "/home/bioinf/resources/silva/silva_v138/silva_nr_v138_train_seq.fa.gz"
# BLAST_DB <- "/home/bioinf/resources/silva/silva_v138/blastdb/SILVA_138_SSURef_NR99_tax_silva.fasta"
# TAXA_MAP <- "/home/bioinf/resources/silva/silva_v138/blastdb/taxmap_slv_ssu_ref_nr_138.txt"
# EVALUE <- "1e-15" %>% as.numeric
# MIN_IDENTITY <- "97" %>% as.numeric
# BLOUT <- "/home/epereira/workspace/indicadores_cuencas_2018/output/test_runs.blout"
# METHOD <- "NBC"
# NSLOTS <- "4" %>% as.numeric
# SAVE_WORKSPACE <- "TRUE" %>% as.logical()
# SEQ_MAP <- "~/workspace/indicadores_cuencas_2018/output/asv_run4/asvs_redu.tsv"

###############################################################################
### 3. Load and format data
###############################################################################

INPUT_ASV_TABLE <- read_csv(file = INPUT_ASV_TABLE, col_names = T) %>%
                   column_to_rownames("X1") %>%
                   as.matrix %>%
                   t

INPUT_ASV_TABLE_tdf <- INPUT_ASV_TABLE %>% 
                       t %>% 
                       as.data.frame %>%
                       rownames_to_column("asv")

###############################################################################
### 4. Run taxa annot: Naive Bayes Classifier (NBC)
###############################################################################

if (METHOD == "NBC") {
  
  print("Running NBC ...")
  TAXA <- assignTaxonomy(seqs = INPUT_ASV_TABLE, 
                         refFasta = TRAIN_DB, 
                         multithread = NSLOTS)  
  
  TAXA <- TAXA %>%
          as.data.frame %>%
          rownames_to_column("asv")
  
  INPUT_ASV_TABLE_ANNOT <- right_join(x = TAXA, y = INPUT_ASV_TABLE_tdf, by = "asv")
  
}

###############################################################################
### 5. Run taxa annot: NBC + Exact Matching (EM)
###############################################################################

if (METHOD == "NBCandEM") {
  
  print("Running NBC ...")
  TAXA <- assignTaxonomy(seqs =  INPUT_ASV_TABLE, 
                         refFasta = TRAIN_DB, 
                         multithread = NSLOTS)

  print("Running EM ...")
  TAXA <- addSpecies(taxtab = taxa, 
                     refFasta = REF_DB)
  
  TAXA <- TAXA %>%
          as.data.frame %>%
          rownames_to_column("asv")
  
  INPUT_ASV_TABLE_ANNOT <- right_join(x = TAXA, y = INPUT_ASV_TABLE_tdf, by = "asv")
  
}

###############################################################################
### 6. Run taxa annot: blastn
###############################################################################

if (METHOD == "BLAST") {
 
  print("Running BLAST ...")
  blastn_runner(db = BLAST_DB, 
                input_seqs = INPUT_FASTA,
                blout = BLOUT,
                evalue = EVALUE, 
                min_identity = MIN_IDENTITY, 
                nslots = NSLOTS)
  
  print("Mapping taxonomy to acc ...")
  
  # load blout as df
  BLOUT <- read_tsv(BLOUT, col_names = F)
  colnames(BLOUT) <- c("qseqid", "sseqid", "pident", "length", 
                       "mismatch", "gapopen", "qstart", "qend", 
                       "sstart", "send", "evalue", "bitscore")
  
  # load seq map as df
  SEQ_MAP <- read_tsv(SEQ_MAP, col_names = T)
  
  # load tax map as df
  TAXA_MAP <- read_tsv(TAXA_MAP, col_names = T)
  TAXA_MAP$sseqid <- paste(TAXA_MAP$primaryAccession, TAXA_MAP$start, TAXA_MAP$stop, sep = ".")
  
  # Check if all BLOUT sseqids are in TAXA_MAP
  if (sum(BLOUT$sseqid %in% TAXA_MAP$sseqid) != length(BLOUT$sseqid)) {
    stop("Not all sseqid in TAXA_MAP file")
  }
  
  # Check if all BLOUT qseqids are in SEQ_MAP
  if (sum(BLOUT$qseqid %in% SEQ_MAP$qseqid) != length(BLOUT$sseqid)) {
    stop("Not all sseqid in SEQ_MAP file")
  }
  
  # cross tables
  BLOUT_TAXA_MAPPED <- left_join(x = BLOUT, y = TAXA_MAP, by = "sseqid") 
  SEQ_TAXA_MAPPED <- left_join(x = SEQ_MAP, y = BLOUT_TAXA_MAPPED, by = "qseqid") %>%
                     select(asv, qseqid, sseqid, path, organism_name, pident) 
  
  # add taxonomy to asv table
  INPUT_ASV_TABLE_ANNOT <- right_join(x = SEQ_TAXA_MAPPED, y = INPUT_ASV_TABLE_tdf, by = "asv")
}

###############################################################################
### 7. Save asv annot table
###############################################################################

write.csv(x = INPUT_ASV_TABLE_ANNOT, file = OUTPUT_ASV_TABLE)
print("Output ASV table saved")
