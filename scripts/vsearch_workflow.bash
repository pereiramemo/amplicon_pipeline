###############################################################################
## define variables
###############################################################################

RUN_DIR="$(dirname "$(readlink -f "$0")")"
source "${RUN_DIR}"/config

show_usage(){
  cat <<EOF
  Usage: ${0##*/}

-h, --help  print this help
-i, --input input fasta
-id, --identity clustering identity
-o, --outdir output file
-t, --nslots number of threads
EOF
}

while :; do
  case "${1}" in

    -h|-\?|--help) # Call a "show_help" function to display a synopsis, then
                   # exit.
    show_usage
    exit 1;
    ;;
#############
  -i|--input)
  if [[ -n "${2}" ]]; then
   INPUT="${2}"
   shift
  fi
  ;;
  --input=?*)
  INPUT="${1#*=}" # Delete everything up to "=" and assign the remainder.
  ;;
  --input=) # Handle the empty case
  printf "ERROR: --input requires a non-empty option argument.\n"  >&2
  exit 1
  ;;
#############
  -id|--identity)
  if [[ -n "${2}" ]]; then
   ID="${2}"
   shift
  fi
  ;;
  --identity=?*)
  ID="${1#*=}" # Delete everything up to "=" and assign the remainder.
  ;;
  --identity=) # Handle the empty case
  printf "ERROR: --identity requires a non-empty option argument.\n"  >&2
  exit 1
  ;;
#############
  -o|--outdir)
   if [[ -n "${2}" ]]; then
     OUTDIR="${2}"
     shift
   fi
  ;;
  --outdir=?*)
  OUTDIR="${1#*=}" # Delete everything up to "=" and assign the
# remainder.
  ;;
  --outdir=) # Handle the empty case
  printf 'Using default environment.\n' >&2
  ;;
#############
  -t|--nslots)
   if [[ -n "${2}" ]]; then
     NSLOTS="${2}"
     shift
   fi
  ;;
  --nslots=?*)
  NSLOTS="${1#*=}" # Delete everything up to "=" and assign the remainder.
  ;;
  --nslots=) # Handle the empty case
  printf 'Using default environment.\n' >&2
  ;;
#############
-s|--sizein)
   SIZEIN="1"
  ;;
############
    --)              # End of all options.
    shift
    break
    ;;
    -?*)
    printf 'WARN: Unknown option (ignored): %s\n' "$1" >&2
    ;;
    *) # Default case: If no more options then break out of the loop.
    break
    esac
    shift
done


##############################################################################
## define outputs
##############################################################################


DEREP_READS="${OUTDIR}/derep${ID/0./}_vsearch_out.fasta"
DEREP_CLUST="${OUTDIR}/derep${ID/0./}_vsearch_out.clust"

##############################################################################
## OUTs clustering 
##############################################################################

if [[ -z "${SIZEIN}" ]]; then

  "${vsearch}" --cluster_fast \
  "${INPUT}" \
  --id "${ID}" \
  --centroids "${DEREP_READS}" \
  --uc "${DEREP_CLUST}" \
  --otutabout "${DEREP_CLUST}.tbl" \
  --threads "${NSLOTS}"

  if [[ $? -ne "0" ]]; then
    echo  "vsearch 1 failed"
    exit
  fi
fi

if [[ -n "${SIZEIN}" ]]; then

  "${vsearch}" --cluster_fast \
  "${INPUT}" \
  --id "${ID}" \
  --centroids "${DEREP_READS}" \
  --uc "${DEREP_CLUST}" \
  --otutabout "${DEREP_CLUST}.tbl" \
  --threads "${NSLOTS}" \
  --sizein

  if [[ $? -ne "0" ]]; then
    echo  "vsearch 2 failed"
    exit
  fi

fi
