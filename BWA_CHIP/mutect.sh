#!/bin/bash

echo "entering mutect script"

TUMOR_NORMAL=$1
NORMAL_SAMPLE=$2
NORMAL_NAME=$3
INTERVALS_FILE=$4
MUTECT_INPUT=$5
SAMPLE_NAME=$6
BWA_GREF=$7
FUNCOTATOR_SOURCES=$8
TRANSCRIPT_LIST=$9
PARAMETER_FILE="${10}"
FILTERED=${11}


if [ $SLURM_ARRAY_TASK_ID -eq 1 ]; then
         echo "arguments used for the mutect.sh script:
               TUMOR_NORMAL=$1
               NORMAL_SAMPLE=$2
               NORMAL_NAME=$3
               INTERVALS_FILE=$4
               MUTECT_INPUT=$5
               SAMPLE_NAME=$6
               BWA_GREF=$7
               FUNCOTATOR_SOURCES=$8
               TRANSCRIPT_LIST=$9
               PARAMETER_FILE=${10}
                " >> $PARAMETER_FILE
fi

INPUTS="--input "$MUTECT_INPUT.bam""
if [ $TUMOR_NORMAL == true ]; then
        INPUTS="${INPUTS} --input ${NORMAL_SAMPLE} --normal $NORMAL_NAME"
        echo "$NORMAL_SAMPLE: $NORMAL_SAMPLE"
fi

OPTIONAL_ARGS="--dont-use-soft-clipped-bases"
if [ $INTERVALS_FILE != false ]; then
        OPTIONAL_ARGS="--intervals $INTERVALS_FILE ${OPTIONAL_ARGS}"
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



if [[ $FILTERED -eq 1 ]] && [[ ! -f "${SAMPLE_NAME}_mutect2_filter_funcotator_coding.vcf" ]] ; then
    grep -E "^#|FRAME_SHIFT_DEL|FRAME_SHIFT_INS|MISSENSE|NONSENSE|SPLICE_SITE" <  "${SAMPLE_NAME}_mutect2_filter_funcotator.vcf" > "${SAMPLE_NAME}_mutect2_filter_funcotator_coding.vcf"

    number_nonsyn_vcf=$(grep -E "FRAME_SHIFT_DEL|FRAME_SHIFT_INS|MISSENSE|NONSENSE|SPLICE_SITE" < "${SAMPLE_NAME}_mutect2_filter_funcotator_coding.vcf" | wc -l)        
    if [ $number_nonsyn_vcf -lt 1 ]; then
        mv "${SAMPLE_NAME}_mutect2_filter_funcotator_coding.vcf" "${SAMPLE_NAME}_mutect2_filter_funcotator_coding_null.vcf"
    fi
elif [ ! -f "${SAMPLE_NAME}_mutect2_filter_funcotator.vcf" ]; then
    
    number_nonsyn_vcf=$(grep -E "FRAME_SHIFT_DEL|FRAME_SHIFT_INS|MISSENSE|NONSENSE|SPLICE_SITE" < "${SAMPLE_NAME}_mutect2_filter_funcotator.vcf" | wc -l)        
    if [ $number_nonsyn_vcf -lt 1 ]; then
        mv "${SAMPLE_NAME}_mutect2_filter_funcotator.vcf" "${SAMPLE_NAME}_mutect2_filter_funcotator_null.vcf"
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

  

if  [[ $FILTERED -eq 1 ]] && [[ ! -f "${SAMPLE_NAME}_mutect2_filter_funcotator_coding.maf" ]]; then
    grep -E "^#|^Hugo_Symbol|Frame_Shift_Del|Frame_Shift_Ins|Missense_Mutation|Nonsense_Mutation|Splice_Site" <  "${SAMPLE_NAME}_mutect2_filter_funcotator.maf" > "${SAMPLE_NAME}_mutect2_filter_funcotator_coding.maf" 
    number_nonsyn_maf=$(grep -E "Frame_Shift_Del|Frame_Shift_Ins|Missense_Mutation|Nonsense_Mutation|Splice_Site" < "${SAMPLE_NAME}_mutect2_filter_funcotator_coding.maf" | wc -l)        
    number_nonsyn_maf=$(grep -E "Frame_Shift_Del|Frame_Shift_Ins|Missense_Mutation|Nonsense_Mutation|Splice_Site" < "${SAMPLE_NAME}_mutect2_filter_funcotator_coding.maf" | wc -l) 
    if [ $number_nonsyn_maf -lt 1 ]; then
        mv "${SAMPLE_NAME}_mutect2_filter_funcotator_coding.maf" "${SAMPLE_NAME}_mutect2_filter_funcotator_coding_null.maf"
    fi

elif [ ! -f "${SAMPLE_NAME}_mutect2_filter_funcotator.maf" ]; then
    number_nonsyn_maf=$(grep -E "Frame_Shift_Del|Frame_Shift_Ins|Missense_Mutation|Nonsense_Mutation|Splice_Site" < "${SAMPLE_NAME}_mutect2_filter_funcotator.maf" | wc -l) 
    if [ $number_nonsyn_maf -lt 1 ]; then
        mv "${SAMPLE_NAME}_mutect2_filter_funcotator.maf" "${SAMPLE_NAME}_mutect2_filter_funcotator_null.maf"
    fi

fi

echo "mutect analysis complete"
