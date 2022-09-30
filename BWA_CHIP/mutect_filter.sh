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
GATK_COMMAND=${10}

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
               GATK_COMMAND=${10}
                " >> $PARAMETER_FILE
fi

cd $OUTPUTS

OUTPUT_NAME="${OUTPUT_DIRECTORY}/${SAMPLE_NAME}"

# Move these two sections to mutect_filter.sh
if [[ ! -f "${OUTPUT_NAME}_contamination.table" ]] ; then
    INPUTS="${SAMPLE_NAME}_pileups.table"

    # Uncomment this when testing tumor-normal
    #if [ "${NORMAL_SAMPLE}" != "false" ]; then
        #INPUTS="${INPUTS} -matched ${NORMAL_PILEUPS}"
    #fi

    echo "Getting contamination rate with CalculateContamination"
    ${GATK_COMMAND} CalculateContamination \
        --input "${INPUTS}" \
        --output "${SAMPLE_NAME}_contamination.table"
else
    echo "Contamination rate already calculated."
    echo "Contamination table: ${SAMPLE_NAME}_contamination.table "
fi

if [ ! -f "${OUTPUT_NAME}_mutect2_filter.vcf" ]; then
    echo "Filtering somatic variants with FilterMutectCalls..."
    ${GATK_COMMAND} FilterMutectCalls \
        --variant "${SAMPLE_NAME}_mutect2.vcf" \
        --output "${SAMPLE_NAME}_mutect2_filter.vcf" \
        --contamination-table ${SAMPLE_NAME}_contamination.table \
        --ob-priors  ${SAMPLE_NAME}_mutect2_artifact_prior.tar.gz \
        --reference "${BWA_GREF}"
    #need to change so is dependent on if need contamination, CANNOT USE QOUTE TRICK

    echo "...somatic variants filtered."
else
    echo "Mutect2 somatic variants already filtered"
fi

echo "mutect_filter.sh complete"
