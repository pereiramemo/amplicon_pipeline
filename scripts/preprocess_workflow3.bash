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
# 
# R1="${DATA}"/osd/OSD1_R1_16S_raw.fastq 
# R2="${DATA}"/osd/OSD1_R2_16S_raw.fastq
# OUTDIR="${RESULTS}"/test

mkdir "${OUTDIR}"

if [[ $? != 0 ]]; then
  echo "mkdir ${OUTDIR} failed"
  exit 1
fi

SAMPLE_NAME=$( basename "${OUTDIR}")
R1_NAME=$(basename "${R1/.fastq/}")
R2_NAME=$(basename "${R2/.fastq/}")

###############################################################################
## 1 - Clip adapters
###############################################################################

R1_CLIPPED_PAIRED="${OUTDIR}/${R1_NAME}"_clipped_paired.fastq
R2_CLIPPED_PAIRED="${OUTDIR}/${R2_NAME}"_clipped_paired.fastq
R1_CLIPPED_UNPAIRED="${OUTDIR}/${R1_NAME}"_clipped_unpaired.fastq
R2_CLIPPED_UNPAIRED="${OUTDIR}/${R2_NAME}"_clipped_unpaired.fastq

java -jar "${trimmomatic}" PE \
  -threads "${NSLOTS}" \
  -trimlog "${OUTDIR}"/pe_trimmomatic.log \
  "${R1}" \
  "${R2}" \
  "${R1_CLIPPED_PAIRED}" \
  "${R1_CLIPPED_UNPAIRED}" \
  "${R2_CLIPPED_PAIRED}" \
  "${R2_CLIPPED_UNPAIRED}" \
  ILLUMINACLIP:"${ADAPTERS}"/TruSeq3-PE.fa:2:30:10:8:true \
  MINLEN:50

if [[ $? != 0 ]]; then
  echo "adapter clipping with ${trimmomatic} failed"
  exit 1
fi

###############################################################################
## 2 - merge with pear
###############################################################################

OUTPREFIX="${OUTDIR}/${R1_NAME/R1_/}"_adapter_clipped

"${pear}" \
-f "${R1_CLIPPED_PAIRED}" \
-r "${R2_CLIPPED_PAIRED}" \
-o "${OUTPREFIX}" \
-j "${NSLOTS}"

if [[ $? != 0 ]]; then
  echo "merge with ${pear} failed"
  exit 1
fi

MERGED="${OUTPREFIX}".assembled.fastq
UNMERGED_FORWARD="${OUTPREFIX}".unassembled.forward.fastq
UNMERGED_REVERSE="${OUTPREFIX}".unassembled.reverse.fastq
DISCARDED="${OUTPREFIX}".discarded.fastq

###############################################################################
## 3 - quality trimming mergeed
###############################################################################

MERGED_QC="${MERGED/.fastq/_qc.fastq}"

"${bbduk}" -Xmx1g \
in="${MERGED}" \
out="${MERGED_QC}" \
qtrim=rl \
minlength=50 \
overwrite=true \
trimq=20 \
threads="${NSLOTS}"

if [[ $? != 0 ]]; then
  echo "quality check 1 with ${bbduk} failed"
  exit 1
fi

###############################################################################
## 4 - convert to fasta
###############################################################################

MERGED_QC_FASTA="${MERGED_QC/.fastq/.fasta}"

"${fq2fa}" "${MERGED_QC}" > "${MERGED_QC_FASTA}"

if [[ $? != 0 ]]; then
  echo "fasta conversion failed"
  exit 1
fi

###############################################################################
## 5 - chimera check
###############################################################################

MERGED_QC_CC_FASTA="${MERGED_QC_FASTA/.fasta/_cc.fasta}"

"${vsearch}" \
--uchime_denovo "${MERGED_QC_FASTA}" \
--nonchimeras "${MERGED_QC_CC_FASTA}" \
--fasta_width 0 \
--abskew 2 \
--threads "${NSLOTS}"

if [[ $? != 0 ]]; then
  echo "chimera check ${vsearch} failed"
  exit 1
fi

###############################################################################
## 6 - Rename sequences: add sample name
###############################################################################

awk -v s="${SAMPLE_NAME}" '{
  if ( $0 ~ ">" ) {
    i++
    sub(">",">Sample_"s"_id_"i"_",$0)
    print $0;
  } else {
    print $0;
  }
}' "${MERGED_QC_CC_FASTA}" > "${MERGED_QC_CC_FASTA/.fasta/_tmp.fasta}"

if [[ $? != 0 ]]; then
  echo "rename headers failed"
  exit 1
fi

mv "${MERGED_QC_CC_FASTA/.fasta/_tmp.fasta}" "${MERGED_QC_CC_FASTA}"

###############################################################################
## 7 - count sequence number and length
###############################################################################

FILES="${R1},${R2},${R1_CLIPPED_PAIRED},${R2_CLIPPED_PAIRED},\
${R1_CLIPPED_UNPAIRED},${R2_CLIPPED_UNPAIRED},${MERGED},${UNMERGED_FORWARD},\
${UNMERGED_REVERSE},${DISCARDED},${MERGED_QC},${MERGED_QC_CC_FASTA}"

COUNTS="${OUTDIR}"/seq_counts.tbl
printf "%s\t%s\t%s\n" "Sample" "n_seq" "aver_seq" > "${COUNTS}"

IFS=","
for F in $( echo "${FILES}" ); do

    NAME=$(basename "${F}")

    if [[ -s "${F}" ]]; then

      if [[ "${F}" =~ ".fastq" ]]; then
        N=$( count_fastq "${F}" )
      fi

      if [[ "${F}" =~ ".fasta" ]]; then
        N=$( count_fasta "${F}" )
      fi

    L=$( "${infoseq}" "${F}" | awk ' NR > 1 { sum = sum + $6 } END { print sum }' )
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
## 8 - clean
###############################################################################

# rm "${R1_CLIPPED_PAIRED}" \
#    "${R2_CLIPPED_PAIRED}" \
#    "${R1_CLIPPED_UNPAIRED}" \
#    "${R2_CLIPPED_UNPAIRED}" \
#    "${MERGED}" \
#    "${UNMERGED_FORWARD}" \
#    "${UNMERGED_REVERSE}" \
#    "${DISCARDED}" \
#    "${MERGED_QC}" \
#    "${MERGED_QC_FASTA}" \
#    "${MERGED_QC_CC_FASTA}"


