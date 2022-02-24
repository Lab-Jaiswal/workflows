#!/bin/bash

echo "running mutect unpaired and with intervals"

PAIRED=$1
NORMAL_SAMPLE=$2
NORMAL_NAME=$3
INTERVALS_FILE=$4
MUTECT_INPUT=$5
SAMPLE_NAME=$6
BWA_GREF=$7
FUNCOTATOR_SOURCES=$8
TRANSCRIPT_LIST=$9

 if [ $PAIRED == false ]; then
        INPUTS="--input "$MUTECT_INPUT.bam""
else 
        INPUTS="--input "$MUTECT_INPUT.bam" --input $NORMAL_SAMPLE.bam --normal $NORMAL_NAME"
fi

if [ $INTERVALS_FILE != false ]; then
        OPTIONAL_ARGS="--intervals $INTERVALS_FILE --dont-use-soft-clipped-bases"
    else
        OPTIONAL_ARGS="--dont-use-soft-clipped-bases"      
fi
    

if [ ! -f "${SAMPLE_NAME}_mutect2.vcf" ]; then
    echo "Calling somatic variants with Mutect2..."
    module load gatk4
    gatk Mutect2 \
    $INPUTS \
    --output "${SAMPLE_NAME}_mutect2.vcf" \
    --reference "${BWA_GREF}" \
    $OPTIONAL_ARGS \
    --bamout "${SAMPLE_NAME}_mutect2.bam"

    #--input "${PREFIX}_${ASSEMBLY}.bam" \
    #-tumor "${PREFIX}" \
    #when have panel of normals, we need to use an argument fo tumor normal mode
    #have --input with normal bam and another with the other bams?
    #one extra args --panel-of-normals (a vcf file we get once made panel of normals)

    module load samtools
    samtools index "${SAMPLE_NAME}_mutect2.bam" "${SAMPLE_NAME}_mutect2.bam.bai"

    echo "...somatic variants called."
else
    echo "Mutect2 somatic variants already called"
fi

if [ ! -f "${SAMPLE_NAME}_mutect2_filter_header.vcf" ]; then
    echo "Filtering somatic variants with FilterMutectCalls..."
    module load gatk4
    gatk FilterMutectCalls \
    --variant "${SAMPLE_NAME}_mutect2.vcf" \
    --output "${SAMPLE_NAME}_mutect2_filter.vcf" \
    --reference "${BWA_GREF}"
    echo "...somatic variants filtered."

    bcftools view -h "${SAMPLE_NAME}_mutect2_filter.vcf" | \
        grep -vP "ID=chr.*_|HLA|EBV" > \
        "${SAMPLE_NAME}_mutect2_filter_header"
    bcftools reheader "${SAMPLE_NAME}_mutect2_filter.vcf" \
        -h "${SAMPLE_NAME}_mutect2_filter_header" \
        -o "${SAMPLE_NAME}_mutect2_filter_reheader.vcf"
else
    echo "Mutect2 somatic variants already filtered"
fi

if [ ! -f "${SAMPLE_NAME}_mutect2_filtr_funcotator.vcf" ]; then
    #FUNCOTATOR_SOURCES="/labs/sjaiswal/tools/funcotator/funcotator_dataSources.v1.7.20200521s"

    echo "Annotating Mutect2 VCF with Funcotator..."
    module load gatk4
    gatk Funcotator \
    --variant "${SAMPLE_NAME}_mutect2_filter_reheader.vcf" \
    --reference "${BWA_GREF}" \
    --ref-version hg38 \
    --data-sources-path "${FUNCOTATOR_SOURCES}" \
    --transcript-list "${TRANSCRIPT_LIST}" \
    --output "${SAMPLE_NAME}_mutect2_filter_funcotator.vcf" \
    --output-file-format VCF 
    echo "...VCF annotated."
else
    echo "Mutect2 VCF already annotated"
fi

if  [ ! -f "${SAMPLE_NAME}_mutect2_filter_funcotator_coding.vcf" ]; then
    grep -E "^#|FRAME_SHIFT_DEL|FRAME_SHIFT_INS|MISSENSE|NONSENSE|SPLICE_SITE" <  "${SAMPLE_NAME}_mutect2_filter_funcotator.vcf" > "${SAMPLE_NAME}_mutect2_filter_funcotator_coding.vcf"
    number_nonsyn_vcf=$(grep -E "FRAME_SHIFT_DEL|FRAME_SHIFT_INS|MISSENSE|NONSENSE|SPLICE_SITE" < "${SAMPLE_NAME}_mutect2_filter_funcotator_coding.vcf" | wc -l)        

    if [ $number_nonsyn_vcf -lt 1 ]; then
        mv "${SAMPLE_NAME}_mutect2_filter_funcotator_coding.vcf" "${SAMPLE_NAME}_mutect2_filter_funcotator_coding_null.vcf"
    fi

fi
    
if  [ ! -f "${SAMPLE_NAME}_mutect2_filter_funcotator.maf" ]; then
    echo "Annotating VCF with Funcotator (MAF output)..."
    module load gatk4
    gatk Funcotator \
    --variant "${SAMPLE_NAME}_mutect2_filter_reheader.vcf" \
    --reference "${BWA_GREF}" \
    --ref-version hg38 \
    --data-sources-path "${FUNCOTATOR_SOURCES}" \
    --transcript-list "${TRANSCRIPT_LIST}" \
    --output "${SAMPLE_NAME}_mutect2_filter_funcotator.maf" \
    --output-file-format MAF 
    echo "...VCF annotated (MAF ouput)."
else
    echo "Mutect2 VCF already annotated (MAF output)"
fi

if  [ ! -f "${SAMPLE_NAME}_mutect2_filter_funcotator_coding.maf" ]; then
    grep -E "^#|^Hugo_Symbol|Frame_Shift_Del|Frame_Shift_Ins|Missense_Mutation|Nonsense_Mutation|Splice_Site" <  "${SAMPLE_NAME}_mutect2_filter_funcotator.maf" > "${SAMPLE_NAME}_mutect2_filter_funcotator_coding.maf" 
    number_nonsyn_maf=$(grep -E "Frame_Shift_Del|Frame_Shift_Ins|Missense_Mutation|Nonsense_Mutation|Splice_Site" < "${SAMPLE_NAME}_mutect2_filter_funcotator_coding.maf" | wc -l)        
    if [ $number_nonsyn_maf -lt 1 ]; then
        mv "${SAMPLE_NAME}_mutect2_filter_funcotator_coding.maf" "${SAMPLE_NAME}_mutect2_filter_funcotator_coding_null.maf"
    fi
fi

echo "mutect analysis complete"
