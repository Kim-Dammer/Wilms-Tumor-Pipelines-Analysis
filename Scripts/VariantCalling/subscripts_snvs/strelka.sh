#! /bin/bash

### run strelka

"/software/strelka/2.9.10/bin/configureStrelkaSomaticWorkflow.py" \
    --normalBam "$CONTROLBAM" \
    --tumorBam "$SAMPLEBAM" \
    --referenceFasta "$REFFASTA" \
    --callRegions "${STRELKA_DIR}/mutect2_callRegions.bed.gz" \
    --runDir "$STRELKA_DIR" \
    &> "${STRELKA_DIR}/STRELKA.log"
"${STRELKA_DIR}/runWorkflow.py" -m local -j 8 &>> "${STRELKA_DIR}/STRELKA.log"
#gunzip "${STRELKA_DIR}/results/variants/somatic.snvs.vcf.gz"
