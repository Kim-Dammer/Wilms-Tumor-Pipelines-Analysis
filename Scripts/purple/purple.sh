#!/bin/bash
set -e

cleanup() {
    echo "Job finished or interrupted. Checking for error logs in /tmp..."

    find /tmp -maxdepth 5 -user "$USER" -mmin -60 \( -name "*err*" -o -name "*.R" \) \
         -exec cp {} /./jobout/ \; 2>/dev/null || true

    echo "Cleanup complete."
}

trap cleanup EXIT

module load GCC/14.1.0
module load R/4.4.2-GCCcore-14.1.0
export R_LIBS_USER="4.4"
export R_LIBS_SITE=""

patientid=2LB-053
sampleid=$1
refid='0'

ID=${patientid}_${sampleid}
normalID=${patientid}_${refid}

basedir=/./Wilms_Tumor


tumor_bam="xxx_merged.mdup.bam"
normal_bam="/controlxxx_merged.mdup.bam"

reference_fasta=/omics/odcf/project/ODCF/reference_genomes/bwa06_1KGRef_PhiX/hs37d5_PhiX.fa


purpledir=$basedir/${patientid}/$sampleid/CNVs/$refid/Purple_re/
cobaltdir=$basedir/${patientid}/$sampleid/CNVs/$refid/Cobalt/
amberdir=$basedir/${patientid}/$sampleid/CNVs/$refid/Amber/

ensembl_data_cache=/./ensembl_data
loci_file=/./AmberGermlineSites.37.vcf.gz
gc_profile=/./GC_profile.1000bp.37.cnp

if [ ! -d $basedir/${patientid}/$sampleid ]; then
    mkdir -p $basedir/${patientid}/$sampleid
    echo "Created sample directory: $basedir/${patientid}/$sampleid"
    echo "For this patient the following samples currently exist: "
    ls $basedir/${patientid}/
fi

if [ ! -d $purpledir ]; then

if [ ! -d $basedir/${patientid}/$sampleid/CNVs ]; then

mkdir $basedir/${patientid}/$sampleid/CNVs

fi

if [ ! -d $basedir/${patientid}/$sampleid/CNVs/$refid ]; then

mkdir $basedir/${patientid}/$sampleid/CNVs/$refid

fi

mkdir $purpledir

fi

if [ ! -d $cobaltdir ]; then
mkdir $cobaltdir
fi

if [ ! -d $amberdir ]; then
mkdir $amberdir
fi

echo $tumor_bam
echo $normal_bam

## I. generate tumor BAF file with Amber

if [ ! -f ${amberdir}/${patientid}_${sampleid}.amber.baf.pcf ]; then

# -min_depth_percent 0.1 \
# -max_depth_percent 10 \
# -min_het_af_percent 0.2 \
/omics/./jdk-22.0.2/bin/java -Xmx32G \
-jar /omics/./amber-3.9.jar \
-ref_genome_version 37 \
-reference ${normalID} \
-reference_bam $normal_bam \
-tumor ${ID} \
-tumor_bam ${tumor_bam} \
-output_dir $amberdir \
-threads 8 \
-loci $loci_file

fi

## II. estimate coverage ratio with Cobalt
#/omics/groups/OE0132/internal/SequenceAlignment/jdk-22.0.2/bin/java -jar -Xmx8G /omics/groups/OE0132/internal/SequenceAlignment/cobalt-1.14.1.jar -reference ${patientid}_${normalid} -reference_bam $normal_bam -tumor ${patientid}_${sampleid} -tumor_bam $tumor_bam -output_dir $cobaltdir -threads 8 -gc_profile /omics/groups/OE0132/internal/SequenceAlignment/Reference_Genome/hg38/hmf_dna_pipeline_resources.38_v5.34/ref/38/copy_number/GC_profile.1000bp.38.cnp

if [ ! -f ${cobaltdir}/${patientid}_${sampleid}.cobalt.ratio.pcf ]; then

/omics/./jdk-22.0.2/bin/java -Xmx8G \
-jar /omics/./cobalt-1.9.jar \
-reference ${normalID} \
-reference_bam $normal_bam \
-tumor ${ID} \
-tumor_bam $tumor_bam \
-output_dir $cobaltdir \
-threads 8 \
-gc_profile $gc_profile

fi

## III. estimate purity and ploidy with Purple
#/omics/groups/OE0132/internal/SequenceAlignment/jdk-22.0.2/bin/java -jar /omics/groups/OE0132/internal/SequenceAlignment/purple_v3.8.4.jar -reference ${patientid}_${normalid} -tumor ${patientid}_${sampleid} -amber $amberdir -cobalt $cobaltdir -gc_profile /omics/groups/OE0132/internal/SequenceAlignment/Reference_Genome/hg38/hmf_dna_pipeline_resources.38_v5.34/ref/38/copy_number/GC_profile.1000bp.38.cnp -ref_genome_version 38 -ref_genome /omics/groups/OE0132/internal/SequenceAlignment/Reference_Genome/hg38/Homo_sapiens.GRCh38.dna_rm.primary_assembly.fa -ensembl_data_dir $ensembl_data_cache -output_dir $purpledir
/omics/./jdk-22.0.2/bin/java \
-jar /omics/./purple_v3.8.4.jar \
-reference ${normalID} \
-tumor ${ID} \
-amber $amberdir \
-cobalt $cobaltdir \
-gc_profile $gc_profile \
-ref_genome_version 37 \
-ref_genome $reference_fasta \
-ensembl_data_dir $ensembl_data_cache \
-output_dir $purpledir
#-min_ploidy 3 \
#-min_purity 0.4 \
#-max_purity 0.8

echo "Purple job completed for sample $sampleid"
