###############################################################################
## 1- define variables
###############################################################################

RUN_DIR="$(dirname "$(readlink -f "$0")")"
source "${RUN_DIR}"/config

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
-outfmt '6 std qlen stitle' \
-perc_identity 75 \
-max_target_seqs 1 \
-evalue 0.0001 \
-num_threads "${NSLOTS}" \
-out "${OUTPUT}"
