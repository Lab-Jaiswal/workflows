#!/bin/bash

echo "entering haplotypecaller.sh script"

SAMPLE_NAME=$1
BWA_GREF=$2
TWIST_SNPS=$3
PARAMETER_FILE=$4
MODE=$5

echo "haplotypecaller command used the following parameters:
$0 $1 $2 $3 $4"

if [ $SLURM_ARRAY_TASK_ID -eq 1 ]; then
        echo "arguments used for the haplotypecaller.sh script:
            SAMPLE_NAME=$1
            BWA_GREF=$2
            TWIST_SNPS=$3
            PARAMETER_FILE=$4
            MODE=$5
             " >> $PARAMETER_FILE
fi

if [ ! -f "${SAMPLE_NAME}_haplotypecaller.gvcf" ]; then
    echo "Calling germline variants with HaplotypeCaller..."
    if [[ $MODE = "slurm" ]]; then
        module load gatk4
    fi
    gatk HaplotypeCaller \
    --input "${SAMPLE_NAME}.bam" \
    --output "${SAMPLE_NAME}_haplotypecaller.gvcf" \
    --reference "${BWA_GREF}" \
    --intervals "${TWIST_SNPS}" \
    --bamout "${SAMPLE_NAME}_haplotypecaller.bam" \
    --emit-ref-confidence GVCF

    if [[ $MODE = "slurm" ]]; then
        module load samtools
    fi
    samtools index "${SAMPLE_NAME}_haplotypecaller.bam" "${SAMPLE_NAME}_haplotypecaller.bam.bai"

    echo "...germline variants called."
else
    echo "HaplotypeCaller gVCF already exists"
fi

if [ ! -f "${SAMPLE_NAME}_haplotypecaller_genotypes.vcf" ]; then
    echo "Genotyping germline variants in gVCF with GenotypeGVCF..."
    if [[ $MODE = "slurm" ]]; then
        module load gatk4
    fi
    gatk GenotypeGVCFs \
    --variant "${SAMPLE_NAME}_haplotypecaller.gvcf" \
    echo "Germline variants already genotyped with GenotypeGVCFs"
fi

echo "haploypecaller.sh complete"
