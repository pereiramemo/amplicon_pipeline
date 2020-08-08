###############################################################################
## 0 - define variables
###############################################################################

RUN_DIR="$(dirname "$(readlink -f "$0")")"
source "${RUN_DIR}"/config

R1="${1}"
R2="${2}"
OUTDIR="${3}"
NSLOTS="${4}"

# R1="${DATA}"/amp_size_248/3_S3_L001_R1_001.fastq
# R2="${DATA}"/amp_size_248/3_S3_L001_R2_001.fastq
# OUTDIR="${RESULTS}"/test


mkdir "${OUTDIR}"

if [[ $? != 0 ]]; then
  echo "mkdir ${OUTDIR} failed"
  exit 1
fi

SAMPLE_NAME=$( basename "${OUTDIR}" )

###############################################################################
## 1 - merge with pear
###############################################################################

"${pear}" \
-f "${R1}" \
-r "${R2}" \
-o "${OUTDIR}/sample_${SAMPLE_NAME}" \
-j "${NSLOTS}"

if [[ $? != 0 ]]; then
  echo "merge with ${pear} failed"
  exit 1
fi

MERGED="${OUTDIR}/sample_${SAMPLE_NAME}.assembled.fastq"
UNMERGED_FORWARD="${OUTDIR}/sample_${SAMPLE_NAME}.unassembled.forward.fastq"
UNMERGED_REVERSE="${OUTDIR}/sample_${SAMPLE_NAME}.unassembled.reverse.fastq"
DISCARDED="${OUTDIR}/sample_${SAMPLE_NAME}.discarded.fastq"

###############################################################################
## 2 - quality trimming
###############################################################################

# merged
MERGED_QC="${MERGED/.fastq/_qc.fastq}"

"${bbduk}" -Xmx1g \
in="${MERGED}" \
out="${MERGED_QC}" \
qtrim=rl \
minlength=100 \
overwrite=true \
trimq=20 \
threads="${NSLOTS}"


if [[ $? != 0 ]]; then
  echo "quality check 1 with ${bbduk} failed"
  exit 1
fi

# unmerged
UNMERGED_FORWARD_QC="${UNMERGED_FORWARD/.fastq/_qc.fastq}"
UNMERGED_REVERSE_QC="${UNMERGED_REVERSE/.fastq/_qc.fastq}"

"${bbduk}" -Xmx1g \
in="${UNMERGED_FORWARD}" \
in2="${UNMERGED_REVERSE}" \
out="${UNMERGED_FORWARD_QC}" \
out2="${UNMERGED_REVERSE_QC}" \
qtrim=rl \
minlength=100 \
overwrite=true \
trimq=20 \
threads="${NSLOTS}"

if [[ $? != 0 ]]; then
  echo "quality check 2 with ${bbduk} failed"
  exit 1
fi

###############################################################################
## 3 - concat all
###############################################################################

CONCAT_QC="${OUTDIR}/sample_${SAMPLE_NAME}_all_qc.fastq"

cat \
"${MERGED_QC}" \
"${UNMERGED_FORWARD_QC}" \
"${UNMERGED_REVERSE_QC}" > "${CONCAT_QC}"

if [[ $? != 0 ]]; then
  echo "concatenation failed"
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
## 5 - dereplication
###############################################################################

CONCAT_QC_DEREP_FASTA="${CONCAT_QC_FASTA/_qc.fasta/_qc_derep.fasta}"

"${vsearch}" \
--derep_prefix "${CONCAT_QC_FASTA}" \
--output "${CONCAT_QC_DEREP_FASTA}" \
--minuniquesize 1 \
--sizeout \
--threads "${NSLOTS}"

if [[ $? != 0 ]]; then
  echo "dereplication ${vsearch} failed"
  exit 1
fi

###############################################################################
## 6 - chimera check
###############################################################################

CONCAT_QC_DEREP_CC_FASTA="${CONCAT_QC_DEREP_FASTA/.fasta/_cc.fasta}"

"${vsearch}" \
--uchime_denovo "${CONCAT_QC_DEREP_FASTA}" \
--nonchimeras "${CONCAT_QC_DEREP_CC_FASTA}" \
--fasta_width 0 \
--abskew 2 \
--threads "${NSLOTS}"

if [[ $? != 0 ]]; then
  echo "chimera check ${vsearch} failed"
  exit 1
fi

###############################################################################
## 7 - count sequence number and length
###############################################################################

FILES_FASTQ="${R1},${R2},${MERGED},${UNMERGED_FORWARD},${UNMERGED_REVERSE},${MERGED_QC},\
${UNMERGED_FORWARD_QC},${UNMERGED_REVERSE_QC},${CONCAT_QC}"

COUNTS="${OUTDIR}"/seq_counts.tbl
printf "%s\t%s\t%s\n" "Sample" "n_seq" "aver_seq" > "${COUNTS}"

IFS=","

for F in $( echo "${FILES_FASTQ}" ); do

  NAME=$(basename "${F}")
  N=$( count_fastq "${F}" )

  L=$( infoseq "${F}" | awk ' NR > 1 { sum = sum + $6 } END { print sum }' )
  A=$( echo  "${L}" / "${N}" | bc -l )

  printf "%s\t%.3f\t%.3f\n" "${NAME}" "${N}" "${A}" >> "${COUNTS}"

done

if [[ $? != 0 ]]; then
  echo "seq counts fastq failed"
  exit 1
fi

FILES_FASTA="${CONCAT_QC_FASTA},${CONCAT_QC_DEREP_FASTA},${CONCAT_QC_DEREP_CC_FASTA}"

for F in $( echo "${FILES_FASTA}" ); do

  NAME=$(basename "${F}")
  N=$( count_fasta "${F}" )

  L=$( infoseq "${F}" | awk ' NR > 1 { sum = sum + $6 } END { print sum }' )
  A=$( echo "${L}" / "${N}" | bc -l )
  printf "%s\t%.3f\t%.3f\n" "${NAME}" "${N}" "${A}" >> "${COUNTS}"

done

if [[ $? != 0 ]]; then
  echo "seq counts fasta failed"
  exit 1
fi


###############################################################################
## 8 - get stats
###############################################################################

"${vsearch}" \
--fastq_stats "${R1}" \
--threads "${NSLOTS}" \
--log "${OUTDIR}"/sample_${SAMPLE_NAME}_R1_stats.tbl

"${vsearch}" \
--fastq_stats "${R2}" \
--threads "${NSLOTS}" \
--log "${OUTDIR}"/sample_${SAMPLE_NAME}_R2_stats.tbl

###############################################################################
## 9 - clean
###############################################################################

rm "${DISCARDED}" \
   "${MERGED}" \
   "${MERGED_QC}" \
   "${UNMERGED_FORWARD}" \
   "${UNMERGED_REVERSE}" \
   "${UNMERGED_FORWARD_QC}" \
   "${UNMERGED_REVERSE_QC}" \
   "${CONCAT_QC_FASTA}"
