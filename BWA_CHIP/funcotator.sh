#!/bin/bash

echo "entering funcotator script"

SAMPLE_NAME=$1
BWA_GREF=$2
FUNCOTATOR_SOURCES=$3
TRANSCRIPT_LIST=$4
PARAMETER_FILE="${5}"
FILTERED=${6}
OUTPUTS=${7}
RUN_FUNCOTATOR=${8}
OUTPUT_DIRECTORY=${9}
MODE=${10}
LINE_NUMBER=${11}
GATK_COMMAND=${12}

if [ $LINE_NUMBER = "1" ]; then
         echo "arguments used for the funcotator.sh script:
               SAMPLE_NAME=$1
               BWA_GREF=$2
               FUNCOTATOR_SOURCES=$3
               TRANSCRIPT_LIST=$4
               PARAMETER_FILE=${5}
               FILTERED=${6}
               OUTPUTS=${7}
               RUN_FUNCOTATOR=${8}
               OUTPUT_DIRECTORY=${9}
               MODE=${10}
               LINE_NUMBER=${11}
               GATK_COMMAND=${12}
                " >> $PARAMETER_FILE
fi

echo "!!!!!!!RUN_FUNCOTATOR: $RUN_FUNCOTATOR !!!!!!!!!!!!!!!^^^^^^^^^"

cd $OUTPUTS

#change outputs to output temp
OUTPUT_NAME="${OUTPUTS}/${SAMPLE_NAME}"
echo " OUTPUT DIRECTORY"

# Move to funcotator.sh
if [ $RUN_FUNCOTATOR = true ]; then
    if [ ! -f "${OUTPUT_NAME}_mutect2_filter_funcotator.vcf" ]; then
        MUTECT_FILTER="${OUTPUT_NAME}_mutect2_filter.vcf"
        NEW_HEADER="${OUTPUT_NAME}_new_header"
        MUTECT_REHEADERED="${OUTPUT_NAME}_mutect2_filter_reheadered.vcf"
                  
        grep "^#" < ${MUTECT_FILTER} | grep -v -E "chrM|chr.*_|HLA|chrEBV" > ${NEW_HEADER} 
        
        echo "incorrect header filtered"

        ${GATK_COMMAND} FixVcfHeader --INPUT ${MUTECT_FILTER} --OUTPUT ${MUTECT_REHEADERED} --HEADER ${NEW_HEADER}

echo "incorrect header replaced"
        echo "Annotating Mutect2 VCF with Funcotator..."
        ${GATK_COMMAND} Funcotator \
        --variant "${MUTECT_REHEADERED}" \
        --reference "${BWA_GREF}" \
        --ref-version hg38 \
        --data-sources-path "${FUNCOTATOR_SOURCES}" \
        --transcript-list "${TRANSCRIPT_LIST}" \
        --output "${OUTPUT_NAME}_mutect2_filter_funcotator.vcf" \
        --output-file-format VCF 
        
        echo "...VCF annotated."
    else
        echo "Mutect2 VCF already annotated"
    fi



    if [[ $FILTERED -eq 1 ]] && [[ ! -f "${OUTPUT_NAME}_mutect2_filter_funcotator_coding.vcf" ]] ; then
        grep -E "^#|FRAME_SHIFT_DEL|FRAME_SHIFT_INS|MISSENSE|NONSENSE|SPLICE_SITE" <  "${OUTPUT_NAME}_mutect2_filter_funcotator.vcf" > "${OUTPUT_NAME}_mutect2_filter_funcotator_coding.vcf"

        number_nonsyn_vcf=$(grep -E "FRAME_SHIFT_DEL|FRAME_SHIFT_INS|MISSENSE|NONSENSE|SPLICE_SITE" < "${OUTPUT_NAME}_mutect2_filter_funcotator_coding.vcf" | wc -l)        
        if [ $number_nonsyn_vcf -lt 1 ]; then
            mv "${OUTPUT_NAME}_mutect2_filter_funcotator_coding.vcf" "${OUTPUT_NAME}_mutect2_filter_funcotator_coding_null.vcf"
        fi
    else 
        number_nonsyn_vcf=$(grep -E "FRAME_SHIFT_DEL|FRAME_SHIFT_INS|MISSENSE|NONSENSE|SPLICE_SITE" < "${OUTPUT_NAME}_mutect2_filter_funcotator.vcf" | wc -l)        
        if [ $number_nonsyn_vcf -lt 1 ]; then
            mv "${OUTPUT_NAME}_mutect2_filter_funcotator.vcf" "${OUTPUT_NAME}_mutect2_filter_funcotator_null.vcf"
        fi

    fi
    
    if  [ ! -f "${OUTPUT_NAME}_mutect2_filter_funcotator.maf" ]; then
        echo "Annotating VCF with Funcotator (MAF output)..."
        ${GATK_COMMAND} Funcotator \
        --variant "${MUTECT_REHEADERED}" \
        --reference "${BWA_GREF}" \
        --ref-version hg38 \
        --data-sources-path "${FUNCOTATOR_SOURCES}" \
        --transcript-list "${TRANSCRIPT_LIST}" \
        --output "${OUTPUT_NAME}_mutect2_filter_funcotator.maf" \
        --output-file-format MAF
            echo "...VCF annotated (MAF ouput)."
    else
        echo "Mutect2 VCF already annotated (MAF output)"
    fi

      

    if  [[ $FILTERED -eq 1 ]] && [[ ! -f "${OUTPUT_NAME}_mutect2_filter_funcotator_coding.maf" ]]; then
        grep -E "^#|^Hugo_Symbol|Frame_Shift_Del|Frame_Shift_Ins|Missense_Mutation|Nonsense_Mutation|Splice_Site" <  "${OUTPUT_NAME}_mutect2_filter_funcotator.maf" > "${OUTPUT_NAME}_mutect2_filter_funcotator_coding.maf" 
        number_nonsyn_maf=$(grep -E "Frame_Shift_Del|Frame_Shift_Ins|Missense_Mutation|Nonsense_Mutation|Splice_Site" < "${OUTPUT_NAME}_mutect2_filter_funcotator_coding.maf" | wc -l)        
        if [ $number_nonsyn_maf -lt 1 ]; then
            mv "${OUTPUT_NAME}_mutect2_filter_funcotator_coding.maf" "${OUTPUT_NAME}_mutect2_filter_funcotator_coding_null.maf"
        fi

    else
        number_nonsyn_maf=$(grep -E "Frame_Shift_Del|Frame_Shift_Ins|Missense_Mutation|Nonsense_Mutation|Splice_Site" < "${OUTPUT_NAME}_mutect2_filter_funcotator.maf" | wc -l) 
        if [ $number_nonsyn_maf -lt 1 ]; then
            mv "${OUTPUT_NAME}_mutect2_filter_funcotator.maf" "${OUTPUT_NAME}_mutect2_filter_funcotator_null.maf"
        fi

    fi
else 
    echo "Funcotator analysis not requested"
fi

echo "funcotator.sh complete"
