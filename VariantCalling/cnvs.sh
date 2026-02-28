#!/bin/bash
set -euo pipefail

# Load required modules
#module load jdk/8u172 gcc/12.2.0  R/4.2.0

# tumor_bam=$(find "$tumorbamdir" -name "*marked_duplicates.bam")

outdir="$1/cnvs"
sample=$2
tumor_bam=$3
refgenome=$4

mkdir -p "$outdir"

scriptdir=/omics/groups/OE0132/tandem/mguenther/Scripts
subscriptdir="${scriptdir}/WGS/VariantCalling/subcripts_cnvs"
readcount_file="${outdir}/tumor_readcounts.wig"

logfile="${sample} cnvs.log"
# shellcheck disable=SC1091
. "${scriptdir}/Logging/logging.sh" "$outdir" "$logfile"
init_logging CNV

###################################################################################################
# READCOUNTS
###################################################################################################

task="readcounts"
if already_done "${task}" ; then
    skip "${task}"
else
    write_progress "Starting ${task}"
    bsub \
        -L /bin/bash \
        -o "${outdir}/hmmcopy_readcounter.log" \
        -n 1 \
        -q verylong \
        -W 10:00 \
        -M 1G \
        -R "rusage[mem=1G]" \
        -J "${sample}_CNVs" \
        "$subscriptdir/readcounts.sh" "${tumor_bam}" "${readcount_file}"
    bwait -w "ended(${sample}_CNVs)" -t 35000
    eval_outcome $? "${task}"
fi

###################################################################################################
# CNV SEGMENTATION
###################################################################################################

task="segmentation"
if already_done "${task}" ; then
    skip "${task}"
else
    write_progress "Starting ${task}"
    Rscript "$subscriptdir/HMMcopy.R" "$outdir" "$refgenome"
    eval_outcome $? "${task}"
    gzip -v "$outdir"/TumorCopy.csv
fi

###################################################################################################

finalize_logging CNV
