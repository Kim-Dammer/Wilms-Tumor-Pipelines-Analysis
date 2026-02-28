#!/bin/bash

dir="${REFGENOME_DIR}/Single_chromosomes"

for chromosome in $(seq 1 "$CHR_MAX") X Y; do
(
    refData=
    for f in "$dir"/*."$chromosome".fa "$dir"/"$chromosome".fa; do
        if [ -e "$f" ]; then
            refData=$f
        fi
    done

    file="${MUTECT_DIR}/temporaryfile_${chromosome}.sh"
    {
        printf "module load GATK/4.6.1.0-Temurin-21.0.3_9\n"
        printf "gatk Mutect2"
        printf " -R %s" "$refData"
        printf " -I %s" "$SAMPLEBAM"
        printf " -I %s" "$CONTROLBAM"
        printf " -L %s" "${chromosome}"
        printf " -normal %s" "$CONTROL_ID"
        printf " --f1r2-tar-gz %s" "${MUTECT_DIR}/${chromosome}_f1r2.tar.gz"
        printf " -O %s" "${MUTECT_DIR}/${chromosome}_unfiltered.vcf.gz"
    } > "$file"
    chmod u+x "$file"

    jobname="${SAMPLE_ID}_mutect_chr${chromosome}"
    bsub \
        -L /bin/bash \
        -o /home/k632r/jobout \
        -n 8 \
        -q verylong \
        -W 240:00 \
        -M 15G \
        -R "rusage[mem=15G]" \
        -J "${jobname}" \
        "$file"
    bwait -w "ended(${jobname})" -t 20000
)&
done
wait
