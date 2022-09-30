#!/bin/bash

echo "entering mutect script"

NORMAL_SAMPLE=$1
INTERVALS_FILE=$2
MUTECT_INPUT=$3
SAMPLE_NAME=$4
BWA_GREF=$5
FUNCOTATOR_SOURCES=$6
TRANSCRIPT_LIST=$7
PARAMETER_FILE="${8}"
FILTERED=${9}
OUTPUTS=${10}
BAM_OUT=${11}
RUN_FUNCOTATOR=${12}
INTERVAL_NUMBER=${13}
#CHR_INTERVALS="/oak/stanford/groups/smontgom/maurertm/workflows/Interval_filtering/whole_genome_intervals.interval_list"
CHR_INTERVALS="/oak/stanford/groups/smontgom/maurertm/workflows/Interval_filtering/chr21_chr22.interval_list"

if [ $SLURM_ARRAY_TASK_ID -eq 1 ]; then
         echo "arguments used for the mutect.sh script:
               NORMAL_SAMPLE=$1
               INTERVALS_FILE=$2
               MUTECT_INPUT=$3
               SAMPLE_NAME=$4
               BWA_GREF=$5
               FUNCOTATOR_SOURCES=$6
               TRANSCRIPT_LIST=$7
               PARAMETER_FILE=${8}
               FILTERED=${9}
               OUTPUTS=${10}
               BAM_OUT=${11}
               RUN_FUNCOTATOR=${12}
               INTERVAL_NUMBER=${13}
                " >> $PARAMETER_FILE
fi

cd $OUTPUTS

INPUTS="--input "$MUTECT_INPUT.bam""
if [ $NORMAL_SAMPLE != false ]; then
    module load samtools
    NORMAL_NAME=$(samtools samples -h $NORMAL_SAMPLE | tail -n 1 | cut -f 1)
    echo "normal name $NORMAL_NAME"

    INPUTS="${INPUTS} --input ${NORMAL_SAMPLE} --normal $NORMAL_NAME"
    echo "$NORMAL_SAMPLE: $NORMAL_SAMPLE"
fi

if [ $INTERVALS_FILE != false ]; then
    OPTIONAL_ARGS="--intervals $INTERVALS_FILE --dont-use-soft-clipped-bases"
else
    mkdir -p Intervals
    interval_line=$(grep -v "@" < $CHR_INTERVALS | sed "${INTERVAL_NUMBER}q; d" ) # Remove header from interval list and choose appropriate line
    interval_name=$(echo "${interval_line}" | cut -f 1) # Extract chromosome name from interval line
    new_interval_file=Intervals/${interval_name}.interval_list
    grep "@" < $CHR_INTERVALS > $new_interval_file # Copy header from original interval list into new interval list
    echo "${interval_line}" >> $new_interval_file # Append line after header in new interval list
    SAMPLE_NAME="${interval_name}_${SAMPLE_NAME}"
    OPTIONAL_ARGS="--intervals $new_interval_file --dont-use-soft-clipped-bases" # Tell Mutect2 to use new interval list
     #OPTIONAL_ARGS="--intervals $new_interval_file --dont-use-soft-clipped-bases \
        #--use-jdk-deflater --use-jdk-inflater --pair-hmm-implementation LOGLESS_CACHING" 
   # OPTIONAL_ARGS="--intervals $new_interval_file --dont-use-soft-clipped-bases \
        #--use-jdk-deflater --use-jdk-inflater --pair-hmm-implementation LOGLESS_CACHING \
        #--java-options '-XX:ParallelGCThreads=1'" # Tell Mutect2 to use new interval list

fi

if [ $BAM_OUT = true ]; then
    OPTIONAL_ARGS="$OPTIONAL_ARGS --bamout ${SAMPLE_NAME}_mutect2.bam"
else
    echo "BAM output not requested"
fi

if [ ! -f "${SAMPLE_NAME}_mutect2.vcf" ]; then
    module load gatk4

    echo "Calling somatic variants with Mutect2 with the following command:
            gatk Mutect2 \
            $INPUTS \
            --output "${SAMPLE_NAME}_mutect2.vcf" \
            --reference "${BWA_GREF}" \
            "${OPTIONAL_ARGS}" \
            --annotation OrientationBiasReadCounts"
    
    gatk Mutect2 \
    $INPUTS \
    --output "${SAMPLE_NAME}_mutect2.vcf" \
    --reference "${BWA_GREF}" \
    $OPTIONAL_ARGS \
    --annotation OrientationBiasReadCounts

    if [ $BAM_OUT = true ]; then
        module load samtools
        samtools index "${SAMPLE_NAME}_mutect2.bam" "${SAMPLE_NAME}_mutect2.bam.bai"
    fi

    echo "...somatic variants called."
else
    echo "Mutect2 somatic variants already called"
fi

if [ $INTERVALS_FILE = false ]; then
    mkdir -p vcfs
    mkdir -p pileups
    OUTPUT_NAME="vcfs/${SAMPLE_NAME}"
    PILEUP_NAME="pileups/${SAMPLE_NAME}"
else
    OUTPUT_NAME=${SAMPLE_NAME}
    PILEUP_NAME="${SAMPLE_NAME}"
fi

if [ ! -f "${PILEUP_NAME}_pileups.table" ]; then
    if [ $INTERVALS_FILE != "false" ]; then
        pileup_intervals="${INTERVALS_FILE}"
    else
        pileup_intervals="${new_interval_file}"
    fi
    echo "Getting pileup summaries with GetPileupSummaries..."
    module load gatk4
    # Need to make gnome_common_variants.vcf.bgz a variable
    gatk GetPileupSummaries \
        $INPUTS \
        -V gnomad_common_variants.vcf.bgz \
        -I "${pileup_intervals}" \
        -O "${PILEUP_NAME}_pileups.table"
else
    echo "Pileup summaries already calculated."
fi

# Move these two sections to mutect_filter.sh
if [ ! -f "${SAMPLE_NAME}_contamination.table" ]; then
    echo "Getting contamination rate with CalculateContamination"
    module load gatk4
    gatk CalculateContamination \
        -I "${SAMPLE_NAME}_pileups.table" \
        -matched ${NORMAL_PILEUPS} \
        -O "${SAMPLE_NAME}_contamination.table"
else
    echo "Contamination rate already calculated."
fi

if [ ! -f "${SAMPLE_NAME}_mutect2_filter.vcf" ]; then
    echo "Filtering somatic variants with FilterMutectCalls..."
    module load gatk4
    gatk FilterMutectCalls \
    --variant "${SAMPLE_NAME}_mutect2.vcf" \
    --output "${SAMPLE_NAME}_mutect2_filter.vcf" \
    --contamination-table "${SAMPLE_NAME}_contamination.table" \
    --reference "${BWA_GREF}"
    echo "...somatic variants filtered."
else
    echo "Mutect2 somatic variants already filtered"
fi

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

#if [ $INTERVALS_FILE != false ]; then
    
#else



#optional contamination step (should it be optional)
#need to get the reference (from Gnomad - look at GATK example in Josh's email to find the right one- it is a vcf that it uses)
#need to get pileups from the bam files first (sum up reads at each position/ seq depth at the snp from the bam file)
#bam_file: ${SAMPLE_NAME}_mutect2.bam
echo "mutect analysis complete"
