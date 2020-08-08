###############################################################################
## 1- define variables
###############################################################################

RUN_DIR="$(dirname "$(readlink -f "$0")")"
source "${RUN_DIR}"/config

R1="${1}"
R2="${2}"
OUTDIR="${3}"

mkdir "${OUTDIR}"

if [[ $? != 0 ]]; then
  echo "mkdir ${OUTDIR} failed"
  exit 1
fi

SAMPLE_NAME=$( basename "${OUTDIR}" )

###############################################################################
## 2 - merge
###############################################################################

"${vsearch}" \
--fastq_mergepairs "${R1}" \
--reverse "${R2}" \
--fastqout  "${OUTDIR}/sample_${SAMPLE_NAME}_assembled.fastq" \
--fastqout_notmerged_fwd \
"${OUTDIR}/sample_${SAMPLE_NAME}_unassembled_forward.fastq" \
--fastqout_notmerged_rev  \
"${OUTDIR}/sample_${SAMPLE_NAME}_unassembled_reverse.fastq" \
--threads "${NSLOTS}" \
 --fastq_minovlen 1 \
--fastq_allowmergestagger \
--fastq_maxdiffs 10

if [[ $? != 0 ]]; then
  echo "merge with ${vsearc} failed"
  exit 1
fi

MERGED="${OUTDIR}/sample_${SAMPLE_NAME}_assembled.fastq"
UNMERGED_FORWARD="${OUTDIR}/sample_${SAMPLE_NAME}_unassembled_forward.fastq"
UNMERGED_REVERSE="${OUTDIR}/sample_${SAMPLE_NAME}_unassembled_reverse.fastq"

###############################################################################
## 3 - concat all
###############################################################################

CONCAT="${OUTDIR}/sample_${SAMPLE_NAME}_all.fastq"

cat \
"${MERGED}" \
"${UNMERGED_FORWARD}" \
"${UNMERGED_REVERSE}" > "${CONCAT}"

if [[ $? != 0 ]]; then
  echo "concatenation failed"
  exit 1
fi

###############################################################################
## 3 - quality trimming 
###############################################################################

CONCAT_QC="${CONCAT/.fastq/_qc.fastq}"

"${vsearch}" \
-fastq_filter "${CONCAT}" \
-fastq_maxee 0.5 \
-fastq_minlen 100 \
-fastaout "${CONCAT_QC}" 

if [[ $? != 0 ]]; then
  echo "quality check 1 with ${vsearch} failed"
  exit 1
fi

###############################################################################
## 4 - convert to fasta
###############################################################################

CONCAT_QC_FASTA="${CONCAT_QC/.fastq/.fasta}"

awk 'NR % 4 == 1 {
       sub("@",">",$0);
       print $0};
     NR % 4 == 2 {
       print $0}' "${CONCAT_QC}" > "${CONCAT_QC_FASTA}"

if [[ $? != 0 ]]; then
  echo "fasta conversion failed"
  exit 1
fi

###############################################################################
## 6 - dereplication
###############################################################################

CONCAT_QC_DEREP="${CONCAT_QC_FASTA/_qc.fasta/_qc_derep.fasta}"

"${vsearch}" \
--derep_prefix "${CONCAT_QC_FASTA}" \
--output "${CONCAT_QC_DEREP}" \
--minuniquesize 1 \
--sizeout

if [[ $? != 0 ]]; then
  echo "dereplication ${vsearch} failed"
  exit 1
fi

###############################################################################
## 7 - chimera check
###############################################################################

CONCAT_QC_DEREP_CC="${CONCAT_QC_DEREP/.fasta/_cc.fasta}"

"${vsearch}" \
--uchime_denovo "${CONCAT_QC_DEREP}" \
--nonchimeras "${CONCAT_QC_DEREP_CC}" \
--fasta_width 0 \
--abskew 1.5

if [[ $? != 0 ]]; then
  echo "chimera check ${vsearch} failed"
  exit 1
fi

###############################################################################
## 8 - count sequences and clean
###############################################################################

FILES_FASTQ="${MERGED},${UNMERGED_FORWARD},${UNMERGED_REVERSE},${CONCAT},\
${CONCAT_QC}"

COUNTS="${OUTDIR}"/seq_counts.tbl

IFS=","
for F in $( echo "${FILES_FASTQ}" ); do

  NAME=$(basename "${F}")
  N=$( count_fastq "${F}" )

  L=$( infoseq "${F}" | awk ' NR > 1 { sum = sum + $6 } END { print sum }' )
  A=$( echo  "${L}" / "${N}" | bc -l )

  echo -e "${NAME}\t${N}\t${A}" >> "${COUNTS}"
  #rm "${F}"

done

if [[ $? != 0 ]]; then
  echo "seq counts fastq failed"
  exit 1
fi

FILES_FASTA="${CONCAT_QC_FASTA},${CONCAT_QC_DEREP},${CONCAT_QC_DEREP_CC}"

for F in $( echo "${FILES_FASTA}" ); do

  NAME=$(basename "${F}")
  N=$( count_fasta "${F}" )

  L=$( infoseq "${F}" | awk ' NR > 1 { sum = sum + $6 } END { print sum }' )
  A=$( echo "${L}" / "${N}" | bc -l )
  echo -e "${NAME}\t${N}\t${A}" >> "${COUNTS}"

done

if [[ $? != 0 ]]; then
  echo "seq counts fasta failed"
  exit 1
fi
#rm "${CONCAT_QC_FASTA}" \
#   "${CONCAT_QC_DEREP}" \
#   "${DISCARDED}"


