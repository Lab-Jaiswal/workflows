#!/bin/bash

echo "entering mutect script"

NORMAL_SAMPLE=$1
INTERVALS_FILE=$2
MUTECT_INPUT=$3
SAMPLE_NAME=$4
BWA_GREF=$5
PARAMETER_FILE="${6}"
OUTPUTS=${7}
BAM_OUT=${8}
OUTPUT_DIRECTORY=${9}
MODE=${10}
SPLIT_BY_CHR=${11}
LINE_NUMBER=${12}
CHR_INTERVALS=${13}
GNOMAD_GENOMES=${14}
RUN_MUTECT=${15}
GATK_COMMAND="${16}"
INTERVAL_NUMBER=${17}

if [ $LINE_NUMBER = "1" ]; then
         echo "###########################################################
               arguments used for the mutect_and_pileups.sh script: 
                    NORMAL_SAMPLE=$1
                    INTERVALS_FILE=$2
                    MUTECT_INPUT=$3
                    SAMPLE_NAME=$4
                    BWA_GREF=$5
                    PARAMETER_FILE=${6}
                    OUTPUTS=${7}
                    BAM_OUT=${8}
                    OUTPUT_DIRECTORY=${9}
                    MODE=${10}
                    SPLIT_BY_CHR=${11}
                    LINE_NUMBER=${12}
                    CHR_INTERVALS=${13}
                    GNOMAD_GENOMES=${14}
                    RUN_MUTECT=${15}
                    GATK_COMMAND="${16}"
                    INTERVAL_NUMBER=${17}
                ########################################################
                " >> $PARAMETER_FILE
fi

cd $OUTPUTS

INPUTS="--input "$MUTECT_INPUT.bam""

if [ $NORMAL_SAMPLE != false ]; then
    NORMAL_NAME=$(samtools samples -h $NORMAL_SAMPLE | tail -n 1 | cut -f 1)
    echo "normal name $NORMAL_NAME"

    INPUTS_MUTECT="${INPUTS} --input ${NORMAL_SAMPLE} --normal $NORMAL_NAME"
    echo "$NORMAL_SAMPLE: $NORMAL_SAMPLE"
else 
    INPUTS_MUTECT="$INPUTS"
fi

if [ $SPLIT_BY_CHR = true ]; then
    mkdir -p Intervals
    interval_line=$(grep -v "@" < $CHR_INTERVALS | sed "${INTERVAL_NUMBER}q; d" ) # Remove header from interval list and choose appropriate line
    interval_name=$(echo "${interval_line}" | cut -f 1) # Extract chromosome name from interval line
    new_interval_file=Intervals/${interval_name}.interval_list
    grep "@" < $CHR_INTERVALS > $new_interval_file # Copy header from original interval list into new interval list
    echo "${interval_line}" >> $new_interval_file # Append line after header in new interval list
    SAMPLE_NAME="${interval_name}_${SAMPLE_NAME}"
    OPTIONAL_ARGS="--intervals $new_interval_file --dont-use-soft-clipped-bases" # Tell Mutect2 to use new interval list
else
    OPTIONAL_ARGS="--intervals $INTERVALS_FILE --dont-use-soft-clipped-bases"
fi

if [ $BAM_OUT = true ]; then
    OPTIONAL_ARGS="$OPTIONAL_ARGS --bamout ${SAMPLE_NAME}_mutect2.bam"
else
    echo "BAM output not requested"
fi

if [ $INTERVALS_FILE = false ]; then
    if [ ! -d ${OUTPUTS}/vcfs ]; then
        mkdir -p vcfs
    fi
    OUTPUT_TEMP_NAME="$OUTPUTS/vcfs/${SAMPLE_NAME}"
    OUTPUT_NAME="$OUTPUT_DIRECTORY/vcfs/${SAMPLE_NAME}"
    
        if [ ! -d ${OUTPUTS}/pileups ]; then
                mkdir -p pileups
        fi
        PILEUP_TEMP_NAME="$OUTPUTS/pileups/${SAMPLE_NAME}"
        PILEUP_NAME="${OUTPUT_DIRECTORY}/pileups/${SAMPLE_NAME}"

else
    OUTPUT_TEMP_NAME=${SAMPLE_NAME}
    OUTPUT_NAME="${OUTPUT_DIRECTORY}/${SAMPLE_NAME}"
    
        PILEUP_TEMP_NAME="${SAMPLE_NAME}"
        PILEUP_NAME="${OUTPUT_DIRECTORY}/${SAMPLE_NAME}"

fi

if [[ ! -f "${OUTPUT_NAME}_mutect2.vcf" ]] && [[ $RUN_MUTECT == true ]]; then
    echo "output name: ${OUTPUT_NAME}_mutect2.vcf"
    echo "mutect analysis requested"
    
    F1R2="${OUTPUTS}/f1r2"
    mkdir -p $F1R2
    F1R2_TEMP_NAME="${OUTPUTS}/f1r2/${SAMPLE_NAME}"
    echo "F1R2 TEMP NAME: $F1R2_TEMP_NAME"

    echo "Calling somatic variants with Mutect2 with the following command:
            gatk Mutect2 \
            $INPUTS_MUTECT \
            --output "${OUTPUT_TEMP_NAME}_mutect2.vcf" \
            --reference "${BWA_GREF}" \
            "${OPTIONAL_ARGS}" \
            --f1r2-tar-gz "${F1R2_TEMP_NAME}_f1r2.tar.gz" \
            --annotation OrientationBiasReadCounts"

    # Next time sleep inf, log in to node.
    
    #gatk Mutect2 \
    ${GATK_COMMAND} Mutect2 \
        $INPUTS_MUTECT \
        --output "${OUTPUT_TEMP_NAME}_mutect2.vcf" \
        --reference "${BWA_GREF}" \
        $OPTIONAL_ARGS \
        --f1r2-tar-gz "${F1R2_TEMP_NAME}_f1r2.tar.gz" \
        --annotation OrientationBiasReadCounts

    if [ $BAM_OUT = true ]; then
        ${GATK_COMMAND} BuildBamIndex \
        --INPUT "${SAMPLE_NAME}_mutect2.bam"
    fi
        echo "...somatic variants called."
else
    echo "output name: ${OUTPUT_TEMP_NAME}_mutect2.vcf"
    echo "Mutect2 somatic variants already called"
fi


if  [[ ! -f "${PILEUP_NAME}_pileups.table" ]]; then
    echo "$PILEUP_NAME: $PILEUP_NAME"
    if [ $INTERVALS_FILE != "false" ]; then
        pileup_intervals="${INTERVALS_FILE}"
    else
        pileup_intervals="${new_interval_file}"
    fi
    echo "Getting pileup summaries with GetPileupSummaries..."
    # Need to make gnome_common_variants.vcf.bgz a variable
    echo "command: ${GATK_COMMAND} GetPileupSummaries \
        $INPUTS \
        -V ${GNOMAD_GENOMES} \
        -L "${pileup_intervals}" \
        -O "${PILEUP_TEMP_NAME}_pileups.table"
"
    #gatk GetPileupSummaries \
   ${GATK_COMMAND} GetPileupSummaries \
        $INPUTS \
        -V ${GNOMAD_GENOMES} \
        -L "${pileup_intervals}" \
        -O "${PILEUP_TEMP_NAME}_pileups.table"
elif [ -f "${PILEUP_NAME}_pileups.table" ]; then
    echo "$PILEUP_NAME: $PILEUP_NAME"
    echo "Pileup summaries already calculated."
else
    echo "Pileup summaries not requested"
fi

echo "mutect_and_pileups.sh analysis complete"
