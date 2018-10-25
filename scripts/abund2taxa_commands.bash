###############################################################################
## 0 - define variables
###############################################################################

RUN_DIR="/bioinf/home/epereira/workspace/16S_analyses/lagunas_16S_analysis2/scripts"
source "${RUN_DIR}"/config 

SWARM_ABUND="${LAGUNAS}/taxa_annot/amp2compare/tables_swarm.tsv"
TAXA_ANNOT="${LAGUNAS}/taxa_annot/amp2compare/tables_silva.tsv"

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
## 2 - modify headers in SWARM_ABUND
###############################################################################

# note: there was a space at the end of line???
SWARM_ABUND_REDU="${SWARM_ABUND/.tsv/_header_redu.tsv}"
cat "${SWARM_ABUND}" |
awk 'BEGIN {FS="\t"; OFS="\t"} {
    $5=gensub(/(.*id_[0-9]+).*/,"\\1","g",$5);
    print $0}' > "${SWARM_ABUND_REDU}"

###############################################################################
## 3 - check number of shared sequences
###############################################################################

head -3 "${TAXA_ANNOT_REDU}"
head -3 "${SWARM_ABUND_REDU}"

NABUND=$(cut -f5 "${SWARM_ABUND_REDU}" | sort | uniq | wc -l)
NTAXA=$(cut -f1 "${TAXA_ANNOT_REDU}" | sort | uniq | wc -l)

echo "${NABUND}" "${NTAXA}"
echo "diff" $(echo "${NABUND}" - "${NTAXA}" | bc)

cat <(cut -f5 "${SWARM_ABUND_REDU}" | sort | uniq) \
    <(cut -f1 "${TAXA_ANNOT_REDU}" | sort | uniq) | \
    sort | uniq -d | wc 

# note: all the taxa annotated reads are in the abundance table (not so obvious)

###############################################################################
## 4 - cross tables
###############################################################################

"${MODULES}"/left_joiner2.perl "${SWARM_ABUND_REDU}" "${TAXA_ANNOT_REDU}" 5 1 | \
egrep -v "NA$" > "${LAGUNAS}/taxa_annot/amp2compare/abund2taxa.tsv"


