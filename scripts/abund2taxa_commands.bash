###############################################################################
## 0 - define variables
###############################################################################

RUN_DIR="/bioinf/home/epereira/workspace/16S_analyses/lagunas_16S_analysis2/scripts"
source "${RUN_DIR}"/config 

PROJECT_NAME="test411"
SWARM_ABUND="${LAGUNAS}/taxa_annot/amp2compare/tables_swarm.tsv"
TAXA_ANNOT="${LAGUNAS}/taxa_annot/amp2compare/tables_silva.tsv"
CLSTR="${LAGUNAS}/taxa_annot/amp2compare/results/ssu/stats/sequence_cluster_map/\
data/${PROJECT_NAME}---ssu---sequence_cluster_map---all_samples_final_representatives.clstr"

###############################################################################
## 1 - reduce TAXA_ANNOT, modify headers and taxonomic annot
###############################################################################

# The first two rows are filtered out, and keep the sequence header and taxa annot
TAXA_TMP=$(mktemp)
HEADER_TMP=$(mktemp)
TAXA_ANNOT_REDU="${TAXA_ANNOT/.tsv/_redu.tsv}"

# remove tabs from taxonomic annot
# note: there were several spaces at the end of line???
tail -n+3 "${TAXA_ANNOT}" | \
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
## 2 - prepare .clster file for taxa annot propagation
###############################################################################

SEQS2REP="${LAGUNAS}/taxa_annot/amp2compare/seqs2rep.tsv"

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

SEQS2TAXA="${LAGUNAS}/taxa_annot/amp2compare/seqs2taxa.tsv"

"${MODULES}"/left_joiner2.perl "${SEQS2REP}" "${TAXA_ANNOT_REDU}" 2 1 | \
egrep -v "NA$" | awk 'BEGIN {OFS="\t"; FS="\t"}; {print $1,$4}' > "${SEQS2TAXA}"

###############################################################################
## 4 - modify headers in SWARM_ABUND
###############################################################################

# note: there was a space at the end of line???
SWARM_ABUND_REDU="${SWARM_ABUND/.tsv/_header_redu.tsv}"
cat "${SWARM_ABUND}" |
awk 'BEGIN {FS="\t"; OFS="\t"} {
    $5=gensub(/(.*id_[0-9]+).*/,"\\1","g",$5);
    print $0}' > "${SWARM_ABUND_REDU}"

###############################################################################
## 5 - check number of shared sequences
###############################################################################

head -3 "${SEQS2TAXA}"
head -3 "${SWARM_ABUND_REDU}"

NABUND=$(cut -f5 "${SWARM_ABUND_REDU}" | sort | uniq | wc -l)
NTAXA=$(cut -f1 "${SEQS2TAXA}" | sort | uniq | wc -l)

echo "${NABUND}" "${NTAXA}"
echo "diff" $(echo "${NABUND}" - "${NTAXA}" | bc)

cat <(cut -f5 "${SWARM_ABUND_REDU}" | sort | uniq) \
    <(cut -f1 "${SEQS2TAXA}" | sort | uniq) | \
    sort | uniq -d | wc 

# note: all the taxa annotated reads are in the abundance table (not so obvious)

###############################################################################
## 6 - cross tables
###############################################################################

CLUST2ABUND2TAXA="${LAGUNAS}/taxa_annot/amp2compare/clust2abund2taxa.tsv"

"${MODULES}"/left_joiner2.perl "${SWARM_ABUND_REDU}" "${SEQS2TAXA}"  5 1 | \
awk 'BEGIN {OFS="\t"; FS="\t"}; { print $1,$2,$3,$4,$5,$7}' > "${CLUST2ABUND2TAXA}"


