###############################################################################
## 1 - define variables
###############################################################################

WORKSPACE="/bioinf/home/epereira/workspace/pipeline_comparison/"

CONFIG="${WORKSPACE}/16S_analysis_pipelines/scripts/config"
source "${CONFIG}"

BIN="${WORKSPACE}/16S_analysis_pipelines/scripts/"
export BIN

SAMPLES=$( ls "${WORKSPACE}"/Cecilia/3363Raw/* | \
           sed "s/_R[1,2].fastq//" | sort | uniq )

###############################################################################
## 2 - run preprocess
###############################################################################

parallel \
--xapply \
-j10 \
-S bigmem-2,bigmem-3,bigmem-4 \
"bash ${BIN}/preprocess_workflow.bash {1}_R1.fastq {1}_R2.fastq {2}" \
::: "${SAMPLES}" ::: $( echo "${SAMPLES}" | sed "s/.*\///")
