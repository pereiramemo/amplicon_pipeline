INPUT_DIR="/bioinf/home/epereira/workspace/pipeline_comparison/\
preprocessed_data_texas/"

INPUT_RAW="/bioinf/home/epereira/workspace/pipeline_comparison/Cecilia/3363Raw/"

OUTPUT_DIR="/bioinf/home/epereira/workspace/pipeline_comparison/\
preprocessed_data_texas/"

SAMPLES=$( echo "$(seq -s " "  1 6) $(seq -s " " 8 32)" )

for i in $( echo ${SAMPLES} ); do
  NSEQ=$( egrep -c ">" \
 <( zcat "${INPUT_DIR}/Bertoglio_texas_sample_${i}.fna.gz" ) )

  MEAN_LENGTH=$( infoseq \
  <(zcat "${INPUT_DIR}/Bertoglio_texas_sample_${i}.fna.gz") | \
  awk 'NR > 1 { sum = sum + $6 } END { print sum / (NR -1 ) }' )

  echo -e "total_seqs\tsample_${i}\t${NSEQ}\t${MEAN_LENGTH}"

  NSEQ_MERGED=$( infoseq \
  <(zcat "${INPUT_DIR}/Bertoglio_texas_sample_${i}.fna.gz") | \
   awk 'BEGIN { OFS="\t" }
    NR > 1 && $6 > 350 { n++; sum = $6 + sum; } END { print n,sum/n }' )

  echo -e "merged_seqs\tsample_${i}\t${NSEQ_MERGED}"

done >  "${OUTPUT_DIR}"/seq_counts.tbl



if [[ -d "${OUTPUT_DIR}/alignments_dir" ]]; then
  mkdir "${OUTPUT_DIR}/alignments_dir"
fi

for i in 1 10 30; do

  egrep ">" <( zcat "${INPUT_DIR}/Bertoglio_texas_sample_${i}.fna.gz" ) | \
  head -10 | sed "s/>//" | \
  while read LINE; do

    FILE=$( echo "${LINE}" | sed "s/\:/_/g" );
    echo "${FILE}"

    $filterbyname \
    in="${INPUT_DIR}/Bertoglio_texas_sample_${i}.fna.gz" \
    out="${OUTPUT_DIR}/alignments_dir/${FILE}.fasta" \
    names=<( echo "${LINE}" ) \
    include=t \
    overwrite=t

    R1="${INPUT_RAW}/3363-${i}-MS28F_R1.fastq"
    R2="${INPUT_RAW}/3363-${i}-MS28F_R2.fastq"

    $filterbyname \
    in="${R1}" \
    in2="${R2}" \
    out="${OUTPUT_DIR}/alignments_dir/${FILE}_r1.fasta" \
    out2="${OUTPUT_DIR}/alignments_dir/${FILE}_r2.fasta" \
    names=<( echo "${LINE/${i}-MS28F::/}" ) \
    include=t \
    overwrite=t


   revseq \
   --sequence "${OUTPUT_DIR}/alignments_dir/${FILE}_r2.fasta" \
   --osformat FASTA \
   --outseq "${OUTPUT_DIR}/alignments_dir/${FILE}_r2_rev.fasta"

   unset MAFFT_BINARIES  
  "${mafft}" <( cat "${OUTPUT_DIR}/alignments_dir/${FILE}.fasta" \
                    "${OUTPUT_DIR}/alignments_dir/${FILE}_r1.fasta" \
                    "${OUTPUT_DIR}/alignments_dir/${FILE}_r2_rev.fasta" ) > \
  "${OUTPUT_DIR}/alignments_dir/${FILE}.aligned"

  done

done