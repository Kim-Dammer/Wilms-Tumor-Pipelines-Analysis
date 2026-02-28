#BSUB -L /bin/bash
#BSUB -n 1 -R "span[hosts=1]"
#BSUB -J Run_purple
#BSUB -q short
#BSUB -W 00:10
#BSUB -M 1G -R "rusage[mem=1G]"


export R_LIBS_USER="/./4.4"

export PATH="/software/R/4.4.2-GCCcore-14.1.0/bin:$PATH"

patientid="2LB-053"
samples=("tumor-01-01" "tumor-01-02" "tumor-01-03" "tumor-01-04" "plasma-01-01" "metastasis-01-01")
scriptdir=/omics/./Scripts/purple

logdir=/home/./jobout

for sample in "${samples[@]}"; do
    echo "Submitting Purple pipeline for $sample..."

    bsub \
        -L /bin/bash \
        -n 8 \
        -R "rusage[mem=100G]" \
        -q verylong \
        -W 240:00 \
        -J "Purple_${sample}" \
        -o ${logdir}/purple_${patientid}_${sample}_%J.out \
        -e ${logdir}/purple_${patientid}_${sample}_%J.err \
        "${scriptdir}/purple.sh" "$sample"
done
