#! /bin/bash

set -euo pipefail

print_usage() {
    cat << 'EOF'

Usage:

    variant_calling.sh [-ahnt] outdir sampleID normalID samplebam normalbam refgenome

Description:

    Perform variant calling for matched-pair `sampleID` and `normalID` using
    `refgenome`. The paths to the BAM files are `samplebam` and `normalbam`,
    and results will be stored in subdirectories of `outdir`.

    The following options are avilable:
      -c    enable calling of CNVs (default)
      -h    display help and exit
      -l    enable calling of LOH (default)
      -n    negates all subsequent options
      -s    enable calling of SNVs (default)
      -v    enable calling of SVs (default)

Examples:

    * Perform both trimming and alignment

        mapping.sh sample fqdir outdir refgenome

    * Omit trimming and only perform alignment

        mapping.sh -nt sample fqdir outdir refgenome

EOF
}

# PARSING FLAGS ###################################################################################

run_cnv=true
run_snv=true
run_sv=true
run_loh=true
negate=false
while getopts 'nclsvh' flag; do
    case "${flag}" in
    n)
    negate=true ;;
    c)
    if [ "$negate" = true ]; then
        run_cnv=false
    fi ;;
    l)
    if [ "$negate" = true ]; then
        run_loh=false
    fi ;;
    s)
    if [ "$negate" = true ]; then
        run_snv=false
    fi ;;
    v)
    if [ "$negate" = true ]; then
        run_sv=false
    fi ;;
    h)  print_usage
        exit 0 ;;
    *)  print_usage
        exit 1 ;;
    esac
done

shift "$((OPTIND - 1))"

# SETUP VARIABLES #################################################################################

outdir=$1
sampleID=$2
normalID=$3
samplebam=$4
normalbam=$5
refgenome=$6
scriptdir=/omics/groups/OE0132/tandem/mguenther/Scripts/WGS/VariantCalling
#scriptdir=/omics/groups/OE0132/internal/SequenceAlignment/ReferenceGenome/hg19/grhc37 ...fa


mkdir -p "$outdir"

# SUBMIT JOBS #####################################################################################

if [ "$run_cnv" = true ]; then
    bsub "${scriptdir}/cnvs.sh" \
        "$outdir" "$sampleID" "$samplebam" "$refgenome"
fi

if [ "$run_snv" = true ]; then
    bsub "${scriptdir}/snvs.sh" \
        "$outdir" "$sampleID" "$normalID" "$samplebam" "$normalbam" "$refgenome"
fi

if [ "$run_loh" = true ]; then
    bsub "${scriptdir}/loh.sh" \
        "$outdir" "$sampleID" "$normalID" "$samplebam" "$normalbam" "$refgenome"
fi

if [ "$run_sv" = true ]; then
    bsub "${scriptdir}/svs.sh" \
        "$outdir" "$sampleID" "$normalID" "$samplebam" "$normalbam" "$refgenome"
fi
