#! /bin/bash

# module load \
#     BEDTools/2.31.1-GCC-14.1.0 \
#     BCFtools/1.20-GCC-14.1.0

### define call regions for strelka based on (unfiltered) mutect output
pad_snv=100
pad_indel=200

bcftools query \
    -f '%CHROM\t%POS0\t%END\t%TYPE\n' \
    "${MUTECT_DIR}/unfiltered.vcf.gz" | \
awk -v ps=$pad_snv -v pi=$pad_indel 'BEGIN{OFS="\t"}{
    pad = ($4=="SNP" ? ps : pi);
    start = $2 - pad; if (start < 0) start = 0;
    end   = $3 + pad;
    print $1, start, end
}' | \
sort -k1,1 -k2,2n | \
bedtools merge > "${STRELKA_DIR}/mutect2_callRegions.bed"

bgzip "${STRELKA_DIR}/mutect2_callRegions.bed"
tabix -p bed "${STRELKA_DIR}/mutect2_callRegions.bed.gz"
