## OSD1
R1=$( readlink -m ../tests/OSD1_R1_16S_raw.fastq )
R2=$( readlink -m ../tests/OSD1_R2_16S_raw.fastq )
OUTDIR=$(readlink -m preprocess_data/OSD1_test )

rm -r "${OUTDIR}"

bash ./scripts/preprocess_workflow.bash $R1 $R2 $OUTDIR

## OSD4
R1=$( readlink -m ../tests/OSD4_R1_16S_raw.fastq )
R2=$( readlink -m ../tests/OSD4_R2_16S_raw.fastq )
OUTDIR=$(readlink -m preprocess_data/OSD4_test )

rm -r "${OUTDIR}"

bash ./scripts/preprocess_workflow.bash $R1 $R2 $OUTDIR

## OSD15
R1=$( readlink -m ../tests/OSD15-50m-depth_R1_16S_raw.fastq )
R2=$( readlink -m ../tests/OSD15-50m-depth_R2_16S_raw.fastq )
OUTDIR=$(readlink -m preprocess_data/OSD15_test )

rm -r "${OUTDIR}"

bash ./scripts/preprocess_workflow.bash $R1 $R2 $OUTDIR




