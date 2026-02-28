#!/bin/bash

set -euo pipefail

module load \
    GATK/4.6.1.0-Temurin-21.0.3_9 \
    BEDTools/2.31.1-GCC-14.1.0 \
    BCFtools/1.20-GCC-14.1.0 \
    SAMtools/1.20-GCC-14.1.0
# shellcheck disable=SC1091
. "$HOME/.bashrc"

outdir=$1/snvs
export SAMPLE_ID=$2
export CONTROL_ID=$3
export SAMPLEBAM=$4
export CONTROLBAM=$5
refgenome=$6

###################################################################################################

export MUTECT_DIR="${outdir}/Mutect"
export STRELKA_DIR="${outdir}/Strelka"

mkdir -p "${outdir}"
mkdir -p "$MUTECT_DIR"
mkdir -p "$STRELKA_DIR"

###################################################################################################
#/omics/odcf/project/ODCF/reference_genomes/bwa06_1KGRef_PhiX
#project/ODCF/reference_genomes/bwa06_1KGRef_PhiX/hs37d5_PhiX.fa
export REFGENOME_DIR="/omics/odcf/project/ODCF/reference_genomes/bwa06_1KGRef_PhiX"

REFFASTA=$(find "$REFGENOME_DIR" \
    -maxdepth 1 \
    -type f \
    -name "hs37d5_PhiX.fa")

export REFFASTA

# export REFGENOME_DIR="/omics/groups/OE0132/internal/SequenceAlignment/Reference_Genome/${refgenome}"
# REFFASTA=$(find "$REFGENOME_DIR" -name "*.primary_assembly.fa") #TODO changed regex
# export REFFASTA
export PROGRAM_DIR="/omics/groups/OE0132/internal/SequenceAlignment"

if [[ "$refgenome" == hg* ]]; then
    CHR_MAX=22
    species=human
elif [[ "$refgenome" == mm* ]]; then
    CHR_MAX=19
    species=mouse
fi
export CHR_MAX

###################################################################################################

scriptdir=/omics/groups/OE0132/tandem/kdammer/Scripts
subscriptdir="${scriptdir}/VariantCalling/subscripts_snvs"

unfiltered_vcf="${MUTECT_DIR}/unfiltered.vcf.gz"
filtered_vcf="${MUTECT_DIR}/filtered.vcf"
repeats_removed_vcf="${MUTECT_DIR}/annovar.remove.repeats.vcf"
repeats_removed_simple_vcf="${MUTECT_DIR}/annovar.remove.repeats.simple.vcf"
mutect_final_vcf="${outdir}/${SAMPLE_ID}_${refgenome}_mutect_final.vcf"

strelka_snv_vcf="${STRELKA_DIR}/results/variants/somatic.snvs.vcf.gz"
strelka_snv_PASS_vcf="${STRELKA_DIR}/results/variants/somatic.snvs.PASS.vcf"
strelka_indel_vcf="${STRELKA_DIR}/results/variants/somatic.indels.vcf.gz"
strelka_indel_PASS_vcf="${STRELKA_DIR}/results/variants/somatic.indels.PASS.vcf"
strelka_final_vcf="${outdir}/${SAMPLE_ID}_${refgenome}_strelka_final.vcf.gz"

###################################################################################################

logfile="${SAMPLE_ID} SNVs.log"
# shellcheck disable=SC1091
. "${scriptdir}/Logging/logging.sh" "$outdir" "$logfile"
init_logging "SNV ANALYSIS"

###################################################################################################
# 1. Mutect2
###################################################################################################

write_header "Mutect2"

# -------------------------------------------------------------------------------------------------
# 1.1. Mutect2 calls per chromosome
# -------------------------------------------------------------------------------------------------

task="mutect2-varcall"
if already_done "${task}" ; then
    skip "${task}"
else
    write_progress "Starting ${task}"
    sh "${subscriptdir}/mutect2.sh"
    eval_outcome $? "${task}"
fi

# -------------------------------------------------------------------------------------------------
# 1.2. Merge unfiltered VCFs
# -------------------------------------------------------------------------------------------------

task="merge-unfiltered-vcfs"
if already_done "${task}" ; then
    skip "${task}"
else
    write_progress "Starting ${task}"
    bcftools concat -Oz -o "$unfiltered_vcf" "${MUTECT_DIR}/"*_unfiltered.vcf.gz
    tabix -p vcf "$unfiltered_vcf"
    eval_outcome $? "${task}"
fi

# -------------------------------------------------------------------------------------------------
# 1.3. Learn Read Orientation Model
# -------------------------------------------------------------------------------------------------

task="read-orientation"
read_orientation_model="${MUTECT_DIR}/read-orientation-model.tar.gz"
if already_done "${task}" ; then
    skip "${task}"
else
    write_progress "Starting ${task}"
    all_f1r2=()
    for chr in $(seq 1 "$CHR_MAX") X Y; do
        file="${MUTECT_DIR}/${chr}_f1r2.tar.gz"
        all_f1r2+=( -I "$file" )
    done
    gatk LearnReadOrientationModel "${all_f1r2[@]}" -O "$read_orientation_model"
    eval_outcome $? "${task}"
fi

# -------------------------------------------------------------------------------------------------
# 1.4. Filter Mutect2 calls per chromosome
# -------------------------------------------------------------------------------------------------

task="mutect2-filtering"
if already_done "${task}" ; then
    skip "${task}"
else
    write_progress "Starting ${task}"
    sh "${subscriptdir}/mutect2_filter.sh" "$read_orientation_model"
    wait
    eval_outcome $? "${task}"
fi

# -------------------------------------------------------------------------------------------------
# 1.5. Merge filtered VCFs
# -------------------------------------------------------------------------------------------------

task="merge-chr-vcfs"
if already_done "${task}" ; then
    skip "${task}"
else
    write_progress "Starting ${task}"
    all_filtered_vcfs=()
    for chr in $(seq 1 "$CHR_MAX") X Y; do
        file="${MUTECT_DIR}/${chr}_filtered.vcf"
        all_filtered_vcfs+=( -I "$file" )
    done
    gatk MergeVcfs "${all_filtered_vcfs[@]}" -O "$filtered_vcf"
    eval_outcome $? "${task}"
fi

###################################################################################################
# 2. Annotate with Annovar
###################################################################################################

write_header "Annovar"

annovar_out_prefix="${MUTECT_DIR}/filtered"
annovar_vcf="${annovar_out_prefix}.${refgenome}_multianno.vcf"

task="annovar"
if already_done "${task}" ; then
    skip "${task}"
else
    write_progress "Starting ${task}"
    "${PROGRAM_DIR}/annovar/table_annovar.pl" \
        "$filtered_vcf" \
        "${PROGRAM_DIR}/annovar/${species}db/" \
        -remove \
        -buildver "$refgenome" \
        -outfile "$annovar_out_prefix" \
        -protocol refGene \
        -operation g \
        -nastring . \
        -vcfinput
    eval_outcome $? "${task}"
fi

###################################################################################################
# 3. VCF filtering
###################################################################################################

write_header "VCF filtering"

# -------------------------------------------------------------------------------------------------
# 3.1. Remove repeats
# -------------------------------------------------------------------------------------------------

task="remove-repeats"
if already_done "${task}" ; then
    skip "${task}"
else
    write_progress "Starting ${task}"
    bash "${subscriptdir}/add_chr_prefix.sh" "$annovar_vcf" "${MUTECT_DIR}/tmp.vcf"
    bedtools intersect \
        -a "${MUTECT_DIR}/tmp.vcf" \
        -b "${REFGENOME_DIR}/Repeats/Repeats_${refgenome}.dms" \
        -v -header \
        > "$repeats_removed_vcf"
    eval_outcome $? "${task}"
fi

# -------------------------------------------------------------------------------------------------
# 3.2. Remove simple repeats
# -------------------------------------------------------------------------------------------------

task="remove-repeats_simple"
if already_done "${task}" ; then
    skip "${task}"
else
    write_progress "Starting ${task}"
    bedtools intersect \
        -a "$repeats_removed_vcf" \
        -b "${REFGENOME_DIR}/Repeats/Simple_Repeats_${refgenome}.dms" \
        -v -header \
        > "$repeats_removed_simple_vcf"
    eval_outcome $? "${task}"
fi

# -------------------------------------------------------------------------------------------------
# 3.3. Finalize Mutect
# -------------------------------------------------------------------------------------------------

task="select-PASS"
if already_done "${task}" ; then
    skip "${task}"
else
    write_progress "Starting ${task}"
    bcftools view -f PASS "$repeats_removed_simple_vcf" > "$mutect_final_vcf"
    eval_outcome $? "${task}"
fi

# -------------------------------------------------------------------------------------------------
# 3.4. Compress VCF files
# -------------------------------------------------------------------------------------------------

task="compress-VCFs"
if already_done "${task}" ; then
    skip "${task}"
else
    write_progress "Starting ${task}"
    bgzip "$mutect_final_vcf"
    eval_outcome $? "${task}"
fi

# -------------------------------------------------------------------------------------------------
# 3.5. Clean-up
# -------------------------------------------------------------------------------------------------

rm -f "${MUTECT_DIR}/temporaryfile"*
rm -f "${MUTECT_DIR}/"*"_unfiltered.vcf.gz"*
rm -f "${MUTECT_DIR}/"*"_filtered.vcf"*
rm -f "${MUTECT_DIR}/"*"_f1r2.tar.gz"
rm -f "${MUTECT_DIR}/tmp.vcf"

###################################################################################################
# 4. SNV calling with strelka
###################################################################################################

write_header "Strelka"

# -------------------------------------------------------------------------------------------------
# 4.1. Define call regions based on Mutect calls
# -------------------------------------------------------------------------------------------------

task="define-callregions"
if already_done "${task}" ; then
    skip "${task}"
else
    write_progress "Starting ${task}"
    sh "${subscriptdir}/strelka_callregions.sh"
    eval_outcome $? "${task}"
fi

# -------------------------------------------------------------------------------------------------
# 4.2. Strelka varcall
# -------------------------------------------------------------------------------------------------

task="strelka-varcall"
if already_done "${task}" ; then
    skip "${task}"
else
    write_progress "Starting ${task}"
    bash "${subscriptdir}/strelka.sh"
    eval_outcome $? "${task}"
fi

# -------------------------------------------------------------------------------------------------
# 4.3. Finalize Strelka
# -------------------------------------------------------------------------------------------------

task="finalize-strelka"
if already_done "${task}" ; then
    skip "${task}"
else
    write_progress "Starting ${task}"
    bash "${subscriptdir}/add_chr_prefix.sh" "$strelka_snv_vcf" "${STRELKA_DIR}/tmp_snv.vcf"
    bash "${subscriptdir}/add_chr_prefix.sh" "$strelka_indel_vcf" "${STRELKA_DIR}/tmp_indel.vcf"
    bcftools view -f PASS "${STRELKA_DIR}/tmp_snv.vcf" > "$strelka_snv_PASS_vcf"
    bcftools view -f PASS "${STRELKA_DIR}/tmp_indel.vcf" > "$strelka_indel_PASS_vcf"
    rm -f "${STRELKA_DIR}/tmp"*
    bgzip "$strelka_snv_PASS_vcf"
    bgzip "$strelka_indel_PASS_vcf"
    bcftools index "${strelka_snv_PASS_vcf}.gz"
    bcftools index "${strelka_indel_PASS_vcf}.gz"
    # merge SNVs & indels
    bcftools concat \
        -a \
        "${strelka_snv_PASS_vcf}.gz" \
        "${strelka_indel_PASS_vcf}.gz" \
        -Oz \
        -o "$strelka_final_vcf"
    #bcftools index "$strelka_final_vcf"
    #bgzip "$strelka_final_vcf"

    eval_outcome $? "${task}"
fi

###################################################################################################
# 5. Mutect-Strelka intersect
###################################################################################################

task="mutect-strelka-intersect"
if already_done "${task}" ; then
    skip "${task}"
else
    write_progress "Starting ${task}"
    bcftools index "${mutect_final_vcf}.gz"
    bcftools index "$strelka_final_vcf"
    bcftools isec -p "$outdir" -Oz "${mutect_final_vcf}.gz" "$strelka_final_vcf"
    mv "${outdir}/0002.vcf.gz" "${outdir}/${SAMPLE_ID}_${refgenome}_intersect.vcf.gz"
    mv "${outdir}/0002.vcf.gz.tbi" "${outdir}/${SAMPLE_ID}_${refgenome}_intersect.vcf.gz.tbi"
    rm -f "${outdir}/000"*".vcf.gz"*
    rm -f "${outdir}/"*.txt
    eval_outcome $? "${task}"
fi

###################################################################################################

finalize_logging "SNV ANALYSIS"
