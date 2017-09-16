###############################################################################
## 1- define variables
###############################################################################

WORKSPACE="/bioinf/home/epereira/workspace/pipeline_comparison/"

CONFIG="${WORKSPACE}"/scripts/config
source "${CONFIG}"

INPUT=${1}
OUTPUT=${2}
NSLOTS=${3}

###############################################################################
## 2- run blast
###############################################################################

DB="${RESOURCES}/SILVA_128_SSURef_Nr99_tax_silva.fasta"

if [[ ! -f "${DB}".nhr ]]; then

  "${makeblastdb}" \
  -in "${DB}" \
  -input_type fasta \
  -parse_seqids \
  -dbtype nucl

fi

"${blastn}" \
-db "${DB}" \
-query "${INPUT}" \
-outfmt '6 std qlen' \
-perc_identity 75 \
-max_target_seqs 1 \
-evalue 0.0001 \
-num_threads "${NSLOTS}" \
-out "${OUTPUT}"
