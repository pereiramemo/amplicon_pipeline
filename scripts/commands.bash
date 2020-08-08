###############################################################################
## 1 - define variables
###############################################################################

WORKSPACE="/bioinf/home/epereira/workspace/pipeline_comparison/"

CONFIG="${WORKSPACE}/16S_analysis_pipelines/scripts/config"
source "${CONFIG}"

BIN="${WORKSPACE}/16S_analysis_pipelines/scripts/"
export BIN
export WORKSPACE

###############################################################################
## 2 - run preprocess
###############################################################################

SAMPLES=$( ls "${WORKSPACE}"/Cecilia/3363Raw/*.fastq | \
           sed "s/_R[1,2].fastq//" | sort | uniq )

parallel \
--xapply \
-j5 \
-S bigmem-1,bigmem-2,bigmem-3,bigmem-4 \
"bash ${BIN}/preprocess_workflow.bash {1}_R1.fastq {1}_R2.fastq \
${WORKSPACE}/16S_analysis_pipelines/preprocess_data/{2}" \
::: "${SAMPLES}" ::: $( echo "${SAMPLES}" | sed "s/.*\///")

###############################################################################
## 3 - concat all preprocessed data
###############################################################################

ALL_FASTA="${WORKSPACE}/16S_analysis_pipelines/preprocess_data/\
all_samples_prepro.fna"

ls "${WORKSPACE}/16S_analysis_pipelines/preprocess_data/"*/*.fasta | \
while read LINE; do

  SAMPLE=$( basename "${LINE}" _all_qc_derep_cc.fasta | \
  sed "s/sample_3363-//");
  sed "s/^>/>${SAMPLE}-/" "${LINE}"

done > "${ALL_FASTA}"

###############################################################################
## 4 - cluster
###############################################################################

OUTPUT_DIR="${WORKSPACE}/16S_analysis_pipelines/results/"
bash "${BIN}/vsearch_workflow.bash" \
--input "${ALL_FASTA}" \
--outdir "${OUTPUT_DIR}" \
--nslots 24 \
--sizein \
--identity 0.98

###############################################################################
## 5 - annotated centroids
###############################################################################

BLOUT="${WORKSPACE}/16S_analysis_pipelines/results/taxa_annot.blout"
CENTROIDS_FASTA="${WORKSPACE}/16S_analysis_pipelines/results/\
derep98_vsearch_out.fasta"

bash "${BIN}/taxa_annotation_wblast.bash" \
"${CENTROIDS_FASTA}" \
"${BLOUT}" \
64

###############################################################################
## 6 - cross tables
###############################################################################

OTU_TABLE="${WORKSPACE}/16S_analysis_pipelines/results/\
derep97_vsearch_out.clust.tbl"

OTU_TABLE_ANNOT="${OTU_TABLE/.tbl/_annot.tbl}"

awk 'BEGIN {FS="\t"; OFS="\t" }; {
  if ( NR == FNR ) {
    array_otu[$1]=$14
    next; }

  if (FNR ==  1) {
    print $0,"taxa";
  }

  if ( array_otu[$1] ) {
    print $0,array_otu[$1]
  }
}' "${BLOUT}" "${OTU_TABLE}" > "${OTU_TABLE_ANNOT}"

