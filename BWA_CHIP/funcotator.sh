#!/bin/bash

echo "entering mutect script"

SAMPLE_NAME=$1
BWA_GREF=$2
FUNCOTATOR_SOURCES=$3
TRANSCRIPT_LIST=$4
PARAMETER_FILE="${5}"
FILTERED=${6}
OUTPUTS=${7}
RUN_FUNCOTATOR=${8}
OUTPUT_DIRECTORY=${9}

if [ $SLURM_ARRAY_TASK_ID -eq 1 ]; then
         echo "arguments used for the mutect.sh script:
               SAMPLE_NAME=$1
               BWA_GREF=$2
               FUNCOTATOR_SOURCES=$3
               TRANSCRIPT_LIST=$4
               PARAMETER_FILE=${5}
               FILTERED=${6}
               OUTPUTS=${7}
               RUN_FUNCOTATOR=${8}
               OUTPUT_DIRECTORY=${9}
                " >> $PARAMETER_FILE
fi

cd $OUTPUTS

OUTPUT_NAME="${OUTPUT_DIRECTORY}/${SAMPLE_NAME}"

# Move to funcotator.sh
if [ $RUN_FUNCOTATOR = true ]; then
    if [ ! -f "${SAMPLE_NAME}_mutect2_filter_funcotator.vcf" ]; then
        echo "Annotating Mutect2 VCF with Funcotator..."
        module load gatk4
        gatk Funcotator \
        --variant "${OUTPUT_NAME}_mutect2_filter.vcf" \
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
    else 
        number_nonsyn_vcf=$(grep -E "FRAME_SHIFT_DEL|FRAME_SHIFT_INS|MISSENSE|NONSENSE|SPLICE_SITE" < "${SAMPLE_NAME}_mutect2_filter_funcotator.vcf" | wc -l)        
        if [ $number_nonsyn_vcf -lt 1 ]; then
            mv "${SAMPLE_NAME}_mutect2_filter_funcotator.vcf" "${SAMPLE_NAME}_mutect2_filter_funcotator_null.vcf"
        fi

    fi
    
    if  [ ! -f "${SAMPLE_NAME}_mutect2_filter_funcotator.maf" ]; then
        echo "Annotating VCF with Funcotator (MAF output)..."
        module load gatk4
        gatk Funcotator \
        --variant "${OUTPUT_NAME}_mutect2_filter.vcf" \
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
        if [ $number_nonsyn_maf -lt 1 ]; then
            mv "${SAMPLE_NAME}_mutect2_filter_funcotator_coding.maf" "${SAMPLE_NAME}_mutect2_filter_funcotator_coding_null.maf"
        fi

    else
        number_nonsyn_maf=$(grep -E "Frame_Shift_Del|Frame_Shift_Ins|Missense_Mutation|Nonsense_Mutation|Splice_Site" < "${SAMPLE_NAME}_mutect2_filter_funcotator.maf" | wc -l) 
        if [ $number_nonsyn_maf -lt 1 ]; then
            mv "${SAMPLE_NAME}_mutect2_filter_funcotator.maf" "${SAMPLE_NAME}_mutect2_filter_funcotator_null.maf"
        fi

    fi
else 
    echo "Funcotator analysis not requested"
fi

echo "funcotator.sh complete"
