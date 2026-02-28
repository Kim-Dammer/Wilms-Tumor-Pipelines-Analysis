#!/bin/bash

read_orientation_model=$1

dir="${REFGENOME_DIR}/Single_chromosomes"

for chromosome in $(seq 1 "$CHR_MAX") X Y; do
(
    refData=
    for f in "$dir"/*."$chromosome".fa "$dir"/"$chromosome".fa; do
        if [ -e "$f" ]; then
            refData=$f
        fi
    done

    file="${MUTECT_DIR}/temporaryfile_${chromosome}_filter.sh"
    {
        printf "module load GATK/4.6.1.0-Temurin-21.0.3_9\n"
        printf "gatk FilterMutectCalls"
        printf " -V %s" "${MUTECT_DIR}/${chromosome}_unfiltered.vcf.gz"
        printf " --ob-priors %s" "${read_orientation_model}"
        printf " --min-median-mapping-quality 30"
        printf " -R %s" "$refData"
        printf " -O %s " "${MUTECT_DIR}/${chromosome}_filtered.vcf"
    } > "$file"
    chmod u+x "$file"

    jobname="${SAMPLE_ID}_filtermutect_chr${chromosome}"
    bsub \
        -L /bin/bash \
        -o /home/k632r/jobout \
        -n 2 \
        -q short \
        -W 0:10 \
        -M 5G \
        -R "rusage[mem=5G]" \
        -J "${jobname}" \
        "$file"
    bwait -w "ended(${jobname})"
)&
done
wait
