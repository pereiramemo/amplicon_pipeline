###############################################################################
## 0 - define variables
###############################################################################

RUN_DIR="$(dirname "$(readlink -f "$0")")"
# RUN_DIR="/bioinf/home/epereira/workspace/16S_analyses/lagunas_16S_analysis2/scripts"
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
## 1 - Trimming by quality
###############################################################################

"${solexa}" dynamictrim "${R1}" "${R2}" \
  --phredcutoff 20 \
  --directory "${OUTDIR}"

if [[ $? != 0 ]]; then
  echo "trim with ${solexa} failed"
  exit 1
fi

###############################################################################
## 2 - Rename solexa quality check output 
###############################################################################

R1_QC_TMP="${OUTDIR}"/$(basename "${R1}" )\.trimmed
R2_QC_TMP="${OUTDIR}"/$(basename "${R2}" )\.trimmed
R1_QC="${OUTDIR}"/$(basename "${R1/.fastq/}" )\_qc.fastq
R2_QC="${OUTDIR}"/$(basename "${R2/.fastq/}" )\_qc.fastq

R1_QC_SEG_TMP="${OUTDIR}"/$(basename "${R1}" )\_trimmed.segments
R2_QC_SEG_TMP="${OUTDIR}"/$(basename "${R2}" )\_trimmed.segments
R1_QC_SEG="${OUTDIR}"/$(basename "${R1/.fastq/}" )\_qc.segments
R2_QC_SEG="${OUTDIR}"/$(basename "${R2/.fastq/}" )\_qc.segments

R1_QC_HIST_TMP="${OUTDIR}"/$(basename "${R1}" )\_trimmed.segments_hist.pdf
R2_QC_HIST_TMP="${OUTDIR}"/$(basename "${R2}" )\_trimmed.segments_hist.pdf
R1_QC_HIST="${OUTDIR}"/$(basename "${R1/.fastq/}" )\_qc.segments_hist.pdf
R2_QC_HIST="${OUTDIR}"/$(basename "${R2/.fastq/}" )\_qc.segments_hist.pdf

mv "${R1_QC_TMP}" "${R1_QC}"
mv "${R2_QC_TMP}" "${R2_QC}"
mv "${R1_QC_SEG_TMP}" "${R1_QC_SEG}"
mv "${R2_QC_SEG_TMP}" "${R2_QC_SEG}"
mv "${R1_QC_HIST_TMP}" "${R1_QC_HIST}"
mv "${R2_QC_HIST_TMP}" "${R2_QC_HIST}"

###############################################################################
## 3 - Length filtering
###############################################################################

"${solexa}" lengthsort "${R1_QC}" "${R2_QC}" \
  --length 50 \
  --directory "${OUTDIR}"

if [[ $? != 0 ]]; then
  echo "length filter with ${bbduk} failed"
  exit 1
fi

###############################################################################
## 4 - Rename solexa length filtering
###############################################################################

R1_QC_LF_TMP="${R1_QC}".paired
R2_QC_LF_TMP="${R2_QC}".paired
SE_QC_LF_TMP="${R1_QC}".single
DISC_QC_LF_TMP="${R1_QC}".discard
SUMMARY_QC_LF_TMP="${R1_QC}".summary.txt
SUMMARY_PDF_QC_LF_TMP="${R1_QC}".summary.txt.pdf

R1_QC_LF="${R1_QC/.fastq/}"_lf_paired.fastq
R2_QC_LF="${R2_QC/.fastq/}"_lf_paired.fastq
SE_QC_LF="${R1_QC/.fastq/}"_lf_single.fastq
DISC_QC_LF="${R1_QC/.fastq/}"_lf_discarded.fastq
SUMMARY_QC_LF="${R1_QC/.fastq/}"_lf_summary.txt
SUMMARY_PDF_QC_LF="${R1_QC/.fastq/}"_lf_summary.txt.pdf

mv "${R1_QC_LF_TMP}" "${R1_QC_LF}"
mv "${R2_QC_LF_TMP}" "${R2_QC_LF}"
mv "${SE_QC_LF_TMP}" "${SE_QC_LF}"
mv "${DISC_QC_LF_TMP}" "${DISC_QC_LF}"
mv "${SUMMARY_QC_LF_TMP}" "${SUMMARY_QC_LF}"
mv "${SUMMARY_PDF_QC_LF_TMP}" "${SUMMARY_PDF_QC_LF}"

###############################################################################
## 5 - trim adapters: PE
###############################################################################

R1_QC_LF_CLIPPED_PAIRED="${OUTDIR}"/$(basename "${R1_QC_LF/.fastq/}" )\_clipped_paired.fastq
R2_QC_LF_CLIPPED_PAIRED="${OUTDIR}"/$(basename "${R2_QC_LF/.fastq/}" )\_clipped_paired.fastq
R1_QC_LF_CLIPPED_UNPAIRED="${OUTDIR}"/$(basename "${R1_QC_LF/.fastq/}" )\_clipped_unpaired.fastq
R2_QC_LF_CLIPPED_UNPAIRED="${OUTDIR}"/$(basename "${R2_QC_LF/.fastq/}" )\_clipped_unpaired.fastq

java -jar "${trimmomatic}" PE \
  -threads "${NSLOTS}" \
  -trimlog "${OUTDIR}"/pe_trimmomatic.log \
  "${R1_QC_LF}" \
  "${R2_QC_LF}" \
  "${R1_QC_LF_CLIPPED_PAIRED}" \
  "${R1_QC_LF_CLIPPED_UNPAIRED}" \
  "${R2_QC_LF_CLIPPED_PAIRED}" \
  "${R2_QC_LF_CLIPPED_UNPAIRED}" \
  ILLUMINACLIP:"${ADAPTERS}"/TruSeq3-PE.fa:2:30:10:3:true

if [[ $? != 0 ]]; then
  echo "pe trimmomatic failed"
  exit 1
fi

###############################################################################
## 6 - trim adapters: SE
###############################################################################

SE_QC_LF_CLIPPED="${OUTDIR}"/$(basename "${SE_QC_LF/.fastq/}" )\_clipped.fastq

if [[ -s  "${SE_QC_LF}" ]]; then

  java -jar "${trimmomatic}" SE \
    -threads "${NSLOTS}" \
    -trimlog "${OUTDIR}"/se_trimmomatic.log \
    "${SE_QC_LF}" \
    "${SE_QC_LF_CLIPPED}" \
    ILLUMINACLIP:"${ADAPTERS}"/TruSeq3-SE.fa:2:30:10

  if [[ $? != 0 ]]; then
    echo "se trimmomatic failed"
    exit 1
  fi
fi

###############################################################################
## 7 - Merge with Flash
###############################################################################

MERGED=$(basename "${R1/R1/merged}" .fastq)

cd "${OUTDIR}"

"${flash}" \
  --threads "${NSLOTS}" \
  --output-prefix "${MERGED}" \
  --max-overlap 200 \
  "${R1_QC_LF_CLIPPED_PAIRED}" \
  "${R2_QC_LF_CLIPPED_PAIRED}"

if [[ $? != 0 ]]; then
  echo "merge with ${flash} failed"
  exit 1
fi

cd "${RUN_DIR}"

###############################################################################
## 8 - count sequence number and length
###############################################################################

FILES_FASTQ="${R1},${R2},\
${R1_QC},${R2_QC},\
${R1_QC_LF},${R2_QC_LF},${SE_QC_LF},${DISC_QC_LF},\
${R1_QC_LF_CLIPPED_PAIRED},${R1_QC_LF_CLIPPED_UNPAIRED},\
${R2_QC_LF_CLIPPED_PAIRED},${R2_QC_LF_CLIPPED_UNPAIRED},\
${SE_QC_LF_CLIPPED},\
${OUTDIR}/${MERGED}.notCombined_1.fastq,\
${OUTDIR}/${MERGED}.notCombined_2.fastq,\
${OUTDIR}/${MERGED}.extendedFrags.fastq"

COUNTS="${OUTDIR}"/seq_counts.tbl
printf "%s\t%s\t%s\n" "Sample" "n_seq" "aver_seq" > "${COUNTS}"

IFS=","

for F in $( echo "${FILES_FASTQ}" ); do

    NAME=$(basename "${F}")

    if [[ -s "${F}" ]]; then
    N=$( count_fastq "${F}" )
    L=$( infoseq "${F}" | awk ' NR > 1 { sum = sum + $6 } END { print sum }' )
    A=$( echo  "${L}" / "${N}" | bc -l )

    printf "%s\t%.3f\t%.3f\n" "${NAME}" "${N}" "${A}" >> "${COUNTS}"
  else
    printf "%s\t%.3f\t%.3f\n" "${NAME}" "0" "0" >> "${COUNTS}"
  fi

done

if [[ $? != 0 ]]; then
  echo "seq counts fastq failed"
  exit 1
fi

###############################################################################
## 9 - Concatenate all reads
###############################################################################

cat "${SE_QC_LF_CLIPPED}"  \
     "${R1_QC_LF_CLIPPED_UNPAIRED}" \
     "${R2_QC_LF_CLIPPED_UNPAIRED}" \
     "${OUTDIR}/${MERGED}".notCombined_1.fastq \
     "${OUTDIR}/${MERGED}".notCombined_2.fastq >> \
     "${OUTDIR}/${MERGED}".extendedFrags.fastq

if [[ $? != 0 ]]; then
  echo "concatenate files failed"
  exit 1
fi

###############################################################################
## 10 - Convert to fasaa
###############################################################################

"${fq2fa}" \
"${OUTDIR}/${MERGED}".extendedFrags.fastq > \
"${OUTDIR}/${MERGED}".extendedFrags.fasta

if [[ $? != 0 ]]; then
  echo " convert to fasta with ${fq2fa} faied"
  exit 1
fi

###############################################################################
## 11 - Rename sequences: add sample name
###############################################################################

sed -i "s/^>/>${SAMPLE_NAME}_/" "${OUTDIR}/${MERGED}".extendedFrags.fasta

if [[ $? != 0 ]]; then
  echo "rename headers failed"
  exit 1
fi

###############################################################################
## 12 - Clean
###############################################################################

rm "${R1_QC}" \
   "${R2_QC}" \
   "${R1_QC_LF}" \
   "${R2_QC_LF}" \
   "${SE_QC_LF}"  \
   "${DISC_QC_LF}" \
   "${SE_QC_LF_CLIPPED}" \
   "${R1_QC_LF_CLIPPED_PAIRED}" \
   "${R2_QC_LF_CLIPPED_PAIRED}" \
   "${R1_QC_LF_CLIPPED_UNPAIRED}" \
   "${R2_QC_LF_CLIPPED_UNPAIRED}" \
   "${OUTDIR}/${MERGED}".notCombined_1.fastq \
   "${OUTDIR}/${MERGED}".notCombined_2.fastq \
   "${OUTDIR}/${MERGED}".extendedFrags.fastq


if [[ $? != 0 ]]; then
  echo "clean failed"
  exit 1
fi