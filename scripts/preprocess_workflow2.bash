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
## 1 - Trimming by quality and length
###############################################################################

R1_QC="${OUTDIR}"/$(basename "${R1}" .fastq)\_qc.fastq
R2_QC="${OUTDIR}"/$(basename "${R2}" .fastq )\_qc.fastq
SE_QC="${OUTDIR}"/$(basename "${R1/R1/SE}" .fastq )\_qc.fastq

"${bbduk}" -Xmx1g \
  in="${R1}" \
  in2="${R2}" \
  out="${R1_QC}" \
  out2="${R2_QC}" \
  outs="${SE_QC}" \
  qtrim=rl \
  minlength=50 \
  overwrite=true \
  trimq=20 \
  threads="${NSLOTS}"

if [[ $? != 0 ]]; then
  echo "quality check with ${bbduk} failed"
  exit 1
fi

###############################################################################
## 2 - trim adapters: PE
###############################################################################

R1_QC_CLIPPED_PAIRED="${OUTDIR}"/$(basename "${R1_QC/.fastq/}" )\_clipped_paired.fastq
R2_QC_CLIPPED_PAIRED="${OUTDIR}"/$(basename "${R2_QC/.fastq/}" )\_clipped_paired.fastq
R1_QC_CLIPPED_UNPAIRED="${OUTDIR}"/$(basename "${R1_QC/.fastq/}" )\_clipped_unpaired.fastq
R2_QC_CLIPPED_UNPAIRED="${OUTDIR}"/$(basename "${R2_QC/.fastq/}" )\_clipped_unpaired.fastq

java -jar "${trimmomatic}" PE \
  -threads "${NSLOTS}" \
  -trimlog "${OUTDIR}"/pe_trimmomatic.log \
  "${R1_QC}" \
  "${R2_QC}" \
  "${R1_QC_CLIPPED_PAIRED}" \
  "${R1_QC_CLIPPED_UNPAIRED}" \
  "${R2_QC_CLIPPED_PAIRED}" \
  "${R2_QC_CLIPPED_UNPAIRED}" \
  ILLUMINACLIP:"${ADAPTERS}"/TruSeq3-PE.fa:2:30:10:3:true

if [[ $? != 0 ]]; then
  echo "pe trimmomatic failed"
  exit 1
fi

###############################################################################
## 3 - trim adapters: SE
###############################################################################

SE_QC_CLIPPED="${OUTDIR}"/$(basename "${SE_QC/.fastq/}" )\_clipped.fastq

if [[ -s  "${SE_QC}" ]]; then

  java -jar "${trimmomatic}" SE \
    -threads "${NSLOTS}" \
    -trimlog "${OUTDIR}"/se_trimmomatic.log \
    "${SE_QC}" \
    "${SE_QC_CLIPPED}" \
    ILLUMINACLIP:"${ADAPTERS}"/TruSeq3-SE.fa:2:30:10

  if [[ $? != 0 ]]; then
    echo "se trimmomatic failed"
    exit 1
  fi
fi

###############################################################################
## 4 - Merge with Flash
###############################################################################

MERGED=$(basename "${R1/R1/merged}" .fastq)

cd "${OUTDIR}"

"${flash}" \
  --threads "${NSLOTS}" \
  --output-prefix "${MERGED}" \
  --max-overlap 200 \
  "${R1_QC_CLIPPED_PAIRED}" \
  "${R2_QC_CLIPPED_PAIRED}"

if [[ $? != 0 ]]; then
  echo "merge with ${flash} failed"
  exit 1
fi

cd "${RUN_DIR}"

###############################################################################
## 5 - Concatenate all preprocessed reads
###############################################################################

cat "${R1_QC_CLIPPED_UNPAIRED}" \
    "${R2_QC_CLIPPED_UNPAIRED}" \
    "${OUTDIR}/${MERGED}".notCombined_1.fastq \
    "${OUTDIR}/${MERGED}".notCombined_2.fastq \
    "${OUTDIR}/${MERGED}".extendedFrags.fastq > \
    "${OUTDIR}/${MERGED}".merged_n_se.fastq

if [[ $? != 0 ]]; then
  echo "concatenate 1 files failed"
  exit 1
fi

if [[ -s  "${SE_QC_CLIPPED}" ]]; then
  cat "${SE_QC_CLIPPED}" >> "${OUTDIR}/${MERGED}".merged_n_se.fastq

  if [[ $? != 0 ]]; then
    echo "concatenate 1 files failed"
    exit 1
  fi

fi

###############################################################################
## 6 - Convert to fasaa
###############################################################################

"${fq2fa}" \
  "${OUTDIR}/${MERGED}".merged_n_se.fastq > \
  "${OUTDIR}/${MERGED}".merged_n_se.fasta

if [[ $? != 0 ]]; then
  echo " convert to fasta with ${fq2fa} faied"
  exit 1
fi

###############################################################################
## 7 - chimera check
###############################################################################

"${vsearch}" \
--uchime_denovo  "${OUTDIR}/${MERGED}".merged_n_se.fasta \
--nonchimeras   "${OUTDIR}/${MERGED}".merged_n_se_cc.fasta \
--fasta_width 0 \
--abskew 2 \
--threads "${NSLOTS}"

###############################################################################
## 8 - count sequence number and length
###############################################################################

FILES="${R1},${R2},\
${R1_QC},${R2_QC},${SE_QC},\
${R1_QC_CLIPPED_PAIRED},${R1_QC_CLIPPED_UNPAIRED},\
${R2_QC_CLIPPED_PAIRED},${R2_QC_CLIPPED_UNPAIRED},\
${SE_QC_CLIPPED},\
${OUTDIR}/${MERGED}.notCombined_1.fastq,\
${OUTDIR}/${MERGED}.notCombined_2.fastq,\
${OUTDIR}/${MERGED}.extendedFrags.fastq,\
"${OUTDIR}/${MERGED}".merged_n_se.fasta,\
"${OUTDIR}/${MERGED}".merged_n_se_cc.fasta"

COUNTS="${OUTDIR}"/seq_counts.tbl
printf "%s\t%s\t%s\n" "Sample" "n_seq" "aver_seq" > "${COUNTS}"

IFS=","

for F in $( echo "${FILES}" ); do

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
## 9 - Rename sequences: add sample name
###############################################################################

sed -i "s/^>/>${SAMPLE_NAME}_/" "${OUTDIR}/${MERGED}".merged_n_se_cc.fasta

if [[ $? != 0 ]]; then
  echo "rename headers failed"
  exit 1
fi

###############################################################################
## 10 - Clean
###############################################################################

rm "${R1_QC}" \
   "${R2_QC}" \
   "${SE_QC}"  \
   "${R1_QC_CLIPPED_PAIRED}" \
   "${R2_QC_CLIPPED_PAIRED}" \
   "${R1_QC_CLIPPED_UNPAIRED}" \
   "${R2_QC_CLIPPED_UNPAIRED}" \
   "${OUTDIR}/${MERGED}".notCombined_1.fastq \
   "${OUTDIR}/${MERGED}".notCombined_2.fastq \
   "${OUTDIR}/${MERGED}".extendedFrags.fastq \
   "${OUTDIR}/${MERGED}".merged_n_se.fastq \
   "${OUTDIR}/${MERGED}".merged_n_se.fasta


if [[ $? != 0 ]]; then
  echo "clean 1 failed"
  exit 1
fi

if [[ -a  "${SE_QC_CLIPPED}"  ]]; then
  rm "${SE_QC_CLIPPED}"

  if [[ $? != 0 ]]; then
    echo "clean 2 failed"
    exit 1
  fi
fi