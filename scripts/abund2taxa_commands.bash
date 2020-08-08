###############################################################################
## 0 - define variables
###############################################################################

RUN_DIR="/bioinf/home/epereira/workspace/16S_analyses/lagunas_16S_analysis2/scripts"
source "${RUN_DIR}"/config 

PROJECT_NAME="test284"
SUBSET="subset_61_84"

SWARM_ABUND="${LAGUNAS}/swarm_clustering/${SUBSET}/all_samples_final_abund.tbl"
TAXA_ANNOT="${LAGUNAS}/taxa_annot/${SUBSET}/${PROJECT_NAME}---ssu---otus.csv"
CLSTR="${LAGUNAS}/taxa_annot/${SUBSET}/\
${PROJECT_NAME}---ssu---sequence_cluster_map---all_samples_final_representatives.clstr"

###############################################################################
## 1 - reduce TAXA_ANNOT, modify headers and taxonomic annot
###############################################################################

# The first two rows are filtered out, and keep the sequence header and taxa annot
TAXA_TMP=$(mktemp)
HEADER_TMP=$(mktemp)
TAXA_ANNOT_REDU="${TAXA_ANNOT/.csv/_redu.csv}"

# remove tabs from taxonomic annot
# note: there were several spaces at the end of line???
tail -n +3 "${TAXA_ANNOT}" | \
cut -f8- | \
tr "\t" " " | \
sed 's/[[:space:]]*$//' > "${TAXA_TMP}"

# simplify header
tail -n+3 "${TAXA_ANNOT}" | \
cut -f3 | \
awk 'BEGIN {FS="\t"; OFS="\t"} {
    $1=gensub(/(.*id_[0-9]+).*/,"\\1","g");
    print $0}' > "${HEADER_TMP}"

# cbind tables
  paste "${HEADER_TMP}" "${TAXA_TMP}" > "${TAXA_ANNOT_REDU}"

# clean
rm  "${HEADER_TMP}" "${TAXA_TMP}" 

###############################################################################
## 2 - prepare CLSTR file for taxa annot propagation
###############################################################################

SEQS2REP="${LAGUNAS}/taxa_annot/${SUBSET}/${PROJECT_NAME}_seqs2rep.tsv"

awk 'BEGIN {OFS ="\t" } {
  if ($1 ~ ">") {
    rep_seq=gensub(/>(.*id_[0-9]+).*/,"\\1","g",$1);
  } else {
    clust_seq = gensub(/>(.*id_[0-9]+).*/,"\\1","g",$3);
    printf "%s\t%s\n", clust_seq,rep_seq;
  }
  
}' "${CLSTR}" > "${SEQS2REP}"

###############################################################################
## 3 - propagate taxa annot to all sequences
###############################################################################

SEQS2TAXA="${LAGUNAS}/taxa_annot/${SUBSET}/${PROJECT_NAME}_seqs2taxa.tsv"

"${MODULES}"/left_joiner2.perl "${SEQS2REP}" "${TAXA_ANNOT_REDU}" 2 1 | \
egrep -v "NA$" | awk 'BEGIN {OFS="\t"; FS="\t"}; {print $1,$4}' > "${SEQS2TAXA}"

###############################################################################
## 4 - modify headers in SWARM_ABUND
###############################################################################

# note: there was a space at the end of line???
SWARM_ABUND_REDU="${LAGUNAS}/taxa_annot/${SUBSET}/\
$(basename ${SWARM_ABUND/.tbl/_header_redu.tsv} )"

cat "${SWARM_ABUND}" |
awk 'BEGIN {FS="\t"; OFS="\t"} {
    $4=gensub(/(.*id_[0-9]+).*/,"\\1","g",$4);
    print $0}' > "${SWARM_ABUND_REDU}"

###############################################################################
## 5 - check number of shared sequences
###############################################################################

head -3 "${SEQS2TAXA}"
head -3 "${SWARM_ABUND_REDU}"

NABUND=$(cut -f4 "${SWARM_ABUND_REDU}" | sort | uniq | wc -l)
NTAXA=$(cut -f1 "${SEQS2TAXA}" | sort | uniq | wc -l)

echo "${NABUND}" "${NTAXA}"
echo "diff" $(echo "${NABUND}" - "${NTAXA}" | bc)

cat <(cut -f4 "${SWARM_ABUND_REDU}" | sort | uniq) \
    <(cut -f1 "${SEQS2TAXA}" | sort | uniq) | \
    sort | uniq -d | wc -l

# note: all the taxa annotated reads are in the abundance table (not so obvious)

###############################################################################
## 6 - cross tables
###############################################################################

CLUST2ABUND2TAXA="${LAGUNAS}/taxa_annot/${SUBSET}/\
${PROJECT_NAME}_${SUBSET}_clust2abund2taxa.tsv"

# Define header
echo -e "cluster\tsample\tabund\tseq_rep\ttaxonomy" > "${CLUST2ABUND2TAXA}"

# Cross tables and format output
"${MODULES}"/left_joiner2.perl "${SWARM_ABUND_REDU}" "${SEQS2TAXA}" 4 1 | \
awk 'BEGIN {OFS="\t"; FS="\t"}; { print $1,$2,$3,$4,$6}' | \
sed "s/\ /_/" >> "${CLUST2ABUND2TAXA}"

###############################################################################
## 7 - clean
###############################################################################

rm "${SEQS2TAXA}" \
   "${SWARM_ABUND_REDU}" \
   "${SEQS2REP}"





