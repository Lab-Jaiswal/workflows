#!/bin/bash

echo "entering mutect script"

SAMPLE_NAME=$1
BWA_GREF=$2
PARAMETER_FILE="${3}"
OUTPUTS=${4}
OUTPUT_DIRECTORY=${5}
NORMAL_SAMPLE=${6}
NORMAL_PILEUPS=${7}
MODE=${8}
LINE_NUMBER=${9}
CONTAINER_ENGINE=${10}
GATK_COMMAND=${11}

if [ $LINE_NUMBER -eq 1 ]; then
         echo "arguments used for the mutect.sh script:
               SAMPLE_NAME=$1
               BWA_GREF=$2
               PARAMETER_FILE=${3}
               OUTPUTS=${4}
               OUTPUT_DIRECTORY=${5}
               NORMAL_SAMPLE=${6}
               NORMAL_PILEUPS=${7}
               MODE=${8}
               LINE_NUMBER=${9}
               CONTAINER_ENGINE=${10}
               GATK_COMMAND=${11}
                " >> $PARAMETER_FILE
fi

cd $OUTPUTS

if [[ $MODE = "slurm" ]]; then
    OUTPUT_NAME="${OUTPUTS}/${SAMPLE_NAME}"
    echo "OUTPUT_NAME" $OUTPUT_NAME
else
    OUTPUT_NAME="${OUTPUT_DIRECTORY}/${SAMPLE_NAME}"
    echo "OUTPUT_NAME: $OUTPUT_NAME"
fi

if [[ $CONTAINER_ENGINE = "singularity" ]]; then
         #INPUT_NAME="${SAMPLE_NAME}"
         INPUT_NAME="${OUTPUT_NAME}"
else
         INPUT_NAME="${OUTPUT_DIRECTORY}/${SAMPLE_NAME}/${SAMPLE_NAME}"
fi

# Move these two sections to mutect_filter.sh
if [[ ! -f "${OUTPUT_NAME}_contamination.table" ]] ; then
    INPUTS="${INPUT_NAME}_pileups.table"

    # Uncomment this when testing tumor-normal
    #if [ "${NORMAL_SAMPLE}" != "false" ]; then
        #INPUTS="${INPUTS} -matched ${NORMAL_PILEUPS}"
    #fi

    echo "Getting contamination rate with CalculateContamination"
    ${GATK_COMMAND} CalculateContamination \
        --input "${INPUTS}" \
        --output "${INPUT_NAME}_contamination.table"
else
    echo "Contamination rate already calculated."
    echo "Contamination table: ${INPUT_NAME}_contamination.table "
fi

if [ ! -f "${OUTPUT_NAME}_mutect2_filter.vcf" ]; then
    echo "Filtering somatic variants with FilterMutectCalls..."
    ${GATK_COMMAND} FilterMutectCalls \
        --variant "${INPUT_NAME}_mutect2.vcf" \
        --output "${INPUT_NAME}_mutect2_filter.vcf" \
        --contamination-table ${INPUT_NAME}_contamination.table \
        --ob-priors  ${INPUT_NAME}_mutect2_artifact_prior.tar.gz \
        --reference "${BWA_GREF}"
    #need to change so is dependent on if need contamination, CANNOT USE QOUTE TRICK

    echo "...somatic variants filtered."
else
    echo "Mutect2 somatic variants already filtered"
fi

echo "mutect_filter.sh complete"
