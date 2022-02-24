#!/bin/bash

echo "entering haplotypecaller.sh script"

SAMPLE_NAME=$1
BWA_GREF=$2
TWIST_SNPS=$3

if [ ! -f "${SAMPLE_NAME}_haplotypecaller.gvcf" ]; then
    echo "Calling germline variants with HaplotypeCaller..."
    module load gatk4
    gatk HaplotypeCaller \
    --input "${SAMPLE_NAME}.bam" \
    --output "${SAMPLE_NAME}_haplotypecaller.gvcf" \
    --reference "${BWA_GREF}" \
    --intervals "${TWIST_SNPS}" \
    --bamout "${SAMPLE_NAME}_haplotypecaller.bam" \
    --emit-ref-confidence GVCF

    module load samtools
    samtools index "${SAMPLE_NAME}_haplotypecaller.bam" "${SAMPLE_NAME}_haplotypecaller.bam.bai"

    echo "...germline variants called."
else
    echo "HaplotypeCaller gVCF already exists"
fi

if [ ! -f "${SAMPLE_NAME}_haplotypecaller_genotypes.vcf" ]; then
    echo "Genotyping germline variants in gVCF with GenotypeGVCF..."
    module load gatk4
    gatk GenotypeGVCFs \
    --variant "${SAMPLE_NAME}_haplotypecaller.gvcf" \
    echo "Germline variants already genotyped with GenotypeGVCFs"
fi

echo "haploypecaller.sh complete"
