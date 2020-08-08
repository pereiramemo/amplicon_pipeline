###############################################################################
## 0 - define variables
###############################################################################

WORKSPACE="/bioinf/home/epereira/workspace/16S_analyses/lagunas_16S_analysis/"

CONFIG="${WORKSPACE}"/scripts/config
source "${CONFIG}"

INPUT="${WORKSPACE}/filtered_data/all_samples.fasta"
# INPUT="${LAGUNAS}/tests/test.fastq"

OUTPUT_DIR=${RESULTS}/cd-hit_based
ID99="0.99"
ID97="0.97"
DEREP_ID99_READS="${OUTPUT_DIR}"/derep"${ID99/0./}"_cdhit_out.fasta
DEREP_ID97_READS="${OUTPUT_DIR}"/derep"${ID97/0./}"_cdhit_out.fasta

##############################################################################
## 1 - OUTs clustering: id 99 
##############################################################################

"${cd_hit_est}" \
  -i "${INPUT}" \
  -o "${DEREP_ID99_READS}" \
  -c "${ID99}" \
  -n 10 \
  -M 0 \
  -d 0 \
  -T "${NSLOTS}"

if [[ $? -ne "0" ]]; then
  echo  "cd-hit-failed id 99"
# exit
fi

##############################################################################
## 2 - OUTs clustering: id 97
##############################################################################

"${cd_hit_est}" \
  -i "${INPUT}" \
  -o "${DEREP_ID97_READS}" \
  -c "${ID97}" \
  -n 10 \
  -M 0 \
  -d 0 \
  -T "${NSLOTS}"

if [[ $? -ne "0" ]]; then
  echo  "cd-hit-failed id 97"
#  exit
fi

##############################################################################
## 3 - make abundance tables
##############################################################################

"${MODULES}"/clstr_parser.awk "${DEREP_ID99_READS}".clstr \
> "${OUTPUT_DIR}"/id"${ID99/0./}"_abund.tbl


"${MODULES}"/clstr_parser.awk "${DEREP_ID97_READS}".clstr \
> "${OUTPUT_DIR}"/id"${ID97/0./}"_abund.tbl
