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
INTERVAL_NUMBER=${14}


#CHR_INTERVALS="/oak/stanford/groups/smontgom/maurertm/ADRC/workflows/Interval_filtering/whole_genome_intervals.interval_list"
#CHR_INTERVALS="/oak/stanford/groups/smontgom/maurertm/ADRC/workflows/Interval_filtering/chr21_chr22.interval_list"

if [ $LINE_NUMBER -eq 1 ]; then
         echo "arguments used for the mutect.sh script:
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
               INTERVAL_NUMBER=${14}
                " >> $PARAMETER_FILE
fi

cd $OUTPUTS

INPUTS="--input "$MUTECT_INPUT.bam""

if [ $NORMAL_SAMPLE != false ]; then
    if [[ $MODE = "slurm" ]]; then
        module load samtools
    fi
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

if [ ! -f "${OUTPUT_NAME}_mutect2.vcf" ]; then
    echo "output name: ${OUTPUT_NAME}_mutect2.vcf"
    echo "mutect analysis requested"
    
    if [[ $MODE = "slurm" ]]; then
        module load gatk4
    fi

    echo "Calling somatic variants with Mutect2 with the following command:
            gatk Mutect2 \
            $INPUTS_MUTECT \
            --output "${OUTPUT_TEMP_NAME}_mutect2.vcf" \
            --reference "${BWA_GREF}" \
            "${OPTIONAL_ARGS}" \
            --annotation OrientationBiasReadCounts"
    
    gatk Mutect2 \
    $INPUTS_MUTECT \
    --output "${OUTPUT_TEMP_NAME}_mutect2.vcf" \
    --reference "${BWA_GREF}" \
    $OPTIONAL_ARGS \
    --annotation OrientationBiasReadCounts

    if [ $BAM_OUT = true ]; then
        if [[ $MODE = "slurm" ]]; then
            module load samtools
        fi
        samtools index "${SAMPLE_NAME}_mutect2.bam" "${SAMPLE_NAME}_mutect2.bam.bai"
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
    if [[ $MODE = "slurm" ]]; then
        module load gatk4
    fi
    # Need to make gnome_common_variants.vcf.bgz a variable
    echo "command: gatk GetPileupSummaries \
        $INPUTS \
        -V "/oak/stanford/groups/smontgom/dnachun/workflows/gnomad.genomes.v3.1.2.sites.maf05.vcf.bgz" \
        -L "${pileup_intervals}" \
        -O "${PILEUP_TEMP_NAME}_pileups.table"
"
    gatk GetPileupSummaries \
        $INPUTS \
        -V "/oak/stanford/groups/smontgom/dnachun/workflows/gnomad.genomes.v3.1.2.sites.maf05.vcf.bgz" \
        -L "${pileup_intervals}" \
        -O "${PILEUP_TEMP_NAME}_pileups.table"
elif [ -f "${PILEUP_NAME}_pileups.table" ]; then
    echo "$PILEUP_NAME: $PILEUP_NAME"
    echo "Pileup summaries already calculated."
else
    echo "Pileup summaries not requested"
fi

echo "mutect_and_pileups.sh analysis complete"
