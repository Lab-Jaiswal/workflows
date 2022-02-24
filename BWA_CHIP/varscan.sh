#!/bin/bash

echo "entering varscan.sh script"

SAMPLE_NAME=$1
BWA_GREF=$2
min_coverage=$3
min_var_freq=$4
p_value=$5

if [ ! -f "${SAMPLE_NAMEY}.pileup" ]; then
    module load samtools
    echo "Generating pileup from BAM..."
    samtools mpileup \
    -A \
    --max-depth 0 \
    -C50 \
    -f "${BWA_GREF}" \
    "${SAMPLE_NAME}.bam" > "${SAMPLE_NAME}.pileup"
    echo "...pileup generated"
else
    echo "Pileup already generated"
fi
if [ ! -f "${SAMPLE_NAME}_varscan2.vcf" ]; then
    echo "Calling variants from pileup..."
    module load varscan
    varscan mpileup2cns \
    "${SAMPLE_NAME}.pileup" \
    --min-coverage "${min_coverage}" \
    --min-var-freq "${min_var_freq}" \
    --p-value "${p_value}" \
    --output-vcf 1 > "${SAMPLE_NAME}_varscan2.vcf"
    echo "...variants called"
else
    echo "Variants already called"
fi

if [ ! -f "${SAMPLE_NAME}_varscan2_filter.vcf" ]; then
    echo "Filtering variants in VCF..."
    varscan filter \
    "${SAMPLE_NAME}_varscan2.vcf" \
    --output-file "${SAMPLE_NAME}_varscan2_filter.vcf" \
    --min-coverage "${min_coverage}" \
    --min-var-freq "${min_var_freq}" \
    --p-value "${p_value}"
    echo "...variants filtered"
else
    echo "Variants already filtered"
fi

if [ ! -f "${SAMPLE_NAME}_varscan2_filter_annovar.hg38_multianno.vcf" ]; then
    echo "Annotating VCF with Annovar..."
    ASSEMBLY_REFGENE="hg38"
    ANNOVARROOT="/labs/sjaiswal/tools/annovar"
    "${ANNOVARROOT}/table_annovar.pl" \
    "${SAMPLE_NAME}_varscan2_filter.vcf" \
    "${ANNOVARROOT}/humandb" \
    --buildver "${ASSEMBLY_REFGENE}" \
    --remove \
    --outfile "${SAMPLE_NAME}_varscan2_filter_annovar" \
    --protocol ensGene \
    --operation g \
    --nastring '.' \
    --vcfinput \
    --thread 1
    echo "...VCF annotated"
else
    echo "VCF already annotated"
fi

echo "varscan analysis complete"
