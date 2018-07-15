################################################################################
## 0 - define variables
###############################################################################

RUN_DIR="$(dirname "$(readlink -f "$0")")"
source "${RUN_DIR}"/config

INPUT="${1}"
TMP_FASTA1="$(mktemp)".fasta
FINAL_FASTA="${INPUT/.fasta/_final.fasta}"

OUTPUT=$( basename "${INPUT}" .fasta)
OUTDIR=$( dirname "${INPUT}")

##############################################
## 1 - derep and format
##############################################

### removes seqs. with Ns
awk '{ if ( $0 ~ /^>/ ) {
         header=$0;
       }
      if ( $0 !~ /^>/ && $0 !~ /[N,n]/ ) {
        printf header"\n"$0"\n";
      }
}' "${INPUT}" > "${TMP_FASTA1}"


### dereplicate
"${vsearch}" \
  --derep_fulllength "${TMP_FASTA1}" \
  --sizeout \
  --fasta_width 0 \
  --output  "${FINAL_FASTA}"

sed -i 's/;size=\(.*\)\;/_\1/' "${FINAL_FASTA}" 

rm "${TMP_FASTA1}" 

##############################################
## 2 - run swarm
##############################################

swarm \
  --differences 1 \
  --fastidious \
  --threads ${NSLOTS} \
  --uclust-file  ${FINAL_FASTA/.fasta/.uclust} \
  --internal-structure ${FINAL_FASTA/.fasta/_1f.struct} \
  --statistics-file ${FINAL_FASTA/.fasta/_1f.stats} \
  --seeds ${FINAL_FASTA/.fasta/_representatives.fasta} \
  --output-file  ${FINAL_FASTA/.fasta/_1f.swarms} < ${FINAL_FASTA}


##############################################
## 3 - make abundance table
##############################################

uclust \
  --uc2clstr ${FINAL_FASTA/.fasta/.uclust} \
  --output ${FINAL_FASTA/.fasta/.clstr}


"${MODULES}"/cd-hit_clstr_parser.awk \
${FINAL_FASTA/.fasta/.clstr} > \
${FINAL_FASTA/.fasta/_abund.tbl}


