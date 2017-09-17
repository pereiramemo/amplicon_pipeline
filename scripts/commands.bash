###############################################################################
## 1 - define variables
###############################################################################

WORKSPACE="/bioinf/home/epereira/workspace/pipeline_comparison/"

CONFIG="${WORKSPACE}/16S_analysis_pipelines/scripts/config"
source "${CONFIG}"

BIN="${WORKSPACE}/16S_analysis_pipelines/scripts/"
export BIN
export WORKSPACE

SAMPLES=$( ls "${WORKSPACE}"/Cecilia/3363Raw/* | \
           sed "s/_R[1,2].fastq//" | sort | uniq  | head )

###############################################################################
## 2 - run preprocess
###############################################################################

parallel \
--xapply \
-j10 \
-S bigmem-2,bigmem-3,bigmem-4 \
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

  SAMPLE=$( basename "${LINE}" _all_qc_derep_cc.fasta | sed "s/sample_//");
  sed "s/^>/${SAMPLE}/" "${LINE}"

done > "${ALL_FASTA}"

###############################################################################
## 4 - cluster
###############################################################################

OUTPUT_DIR="${WORKSPACE}/16S_analysis_pipelines/results/"

"${BIN}/vsearch_workflow.bash"
"${ALL_FASTA}" \
"${OUPUT_DIR}" \
24 \
seizein

###############################################################################
## 5 - annotated centroids
###############################################################################

BLOUT="${WORKSPACE}/16S_analysis_pipelines/results/taxa_annot.blout"
CENTROIDS_FASTA="${WORKSPACE}/16S_analysis_pipelines/results/"

"${BIN}/taxa_annotation_wblast.bash"
"${CENTROIDS_FASTA}" \
"${BLOUT}" \
24

###############################################################################
## 6 - cross tables
###############################################################################
