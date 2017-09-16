###############################################################################
## define variables
###############################################################################

WORKSPACE="/bioinf/home/epereira/workspace/pipeline_comparison/"

CONFIG="${WORKSPACE}"/scripts/config
source "${CONFIG}"

INPUT=${1}
OUTPUT_DIR=${2}
NSLOTS=${3}
SIZEIN=${4}

# INPUT="${LAGUNAS}/tests/test.fastq"

ID99="0.99"
ID97="0.97"
DEREP_ID99_READS="${OUTPUT_DIR}/derep${ID99/0./}_vsearch_out.fasta"
DEREP_ID99_CLUST="${OUTPUT_DIR}/derep${ID99/0./}_vsearch_out.clust"
DEREP_ID97_READS="${OUTPUT_DIR}/derep${ID97/0./}_vsearch_out.fasta"
DEREP_ID97_CLUST="${OUTPUT_DIR}/derep${ID97/0./}_vsearch_out.clust"

##############################################################################
## OUTs clustering: id 99 
##############################################################################

if [[ -z "${SIZEIN}" ]]; then

  "${vsearch}" --cluster_fast \
  "${INPUT}" \
  --id "${ID99}" \
  --centroids "${DEREP_ID99_READS}" \
  --uc "${DEREP_ID99_CLUST}" \
  --otutabout "${DEREP_ID99_CLUST}.tbl" \
  --threads "${NSLOTS}" \
  --relabel OTU_

  if [[ $? -ne "0" ]]; then
    echo  "vsearch failed id 99"
    exit
  fi
fi

if [[ -n "${SIZEIN}" ]]; then

  "${vsearch}" --cluster_fast \
  "${INPUT}" \
  --id "${ID99}" \
  --centroids "${DEREP_ID99_READS}" \
  --uc "${DEREP_ID99_CLUST}" \
  --otutabout "${DEREP_ID99_CLUST}.tbl" \
  --threads "${NSLOTS}" \
  --relabel OTU_ \
  --sizein

  if [[ $? -ne "0" ]]; then
    echo  "vsearch failed id 99"
    exit
  fi

fi

##############################################################################
## OUTs clustering: id 97
##############################################################################

if [[ -z "${SIZEIN}" ]]; then

  "${vsearch}" --cluster_fast \
  "${INPUT}" \
  --id "${ID97}" \
  --centroids "${DEREP_ID97_READS}" \
  --uc "${DEREP_ID97_CLUST}" \
  --otutabout "${DEREP_ID97_CLUST}.tbl" \
  --threads "${NSLOTS}" \
  --relabel OTU_ 

  if [[ $? -ne "0" ]]; then
    echo  "vsearch failed id 97"
    exit 1
  fi

fi

if [[ -n "${SIZEIN}" ]]; then

  "${vsearch}" --cluster_fast \
  "${INPUT}" \
  --id "${ID97}" \
  --centroids "${DEREP_ID97_READS}" \
  --uc "${DEREP_ID97_CLUST}" \
  --otutabout "${DEREP_ID97_CLUST}.tbl" \
  --threads "${NSLOTS}" \
  --relabel OTU_ \
  --sizein

  if [[ $? -ne "0" ]]; then
    echo  "vsearch failed id 97"
    exit 1
  fi

fi

