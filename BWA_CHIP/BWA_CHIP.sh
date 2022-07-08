#!/bin/bash

#SBATCH --time=72:00:00
#SBATCH --account=sjaiswal
#SBATCH --cpus-per-task=8
#SBATCH --mem=124GB
#SBATCH --job-name=CHIP_variant_call

##################################################################################################################################
#############################################---STEP 1: SET UP PARAMETERS---###################################################### 
##################################################################################################################################
PARENT_DIRECTORY=$1
OUTPUT_DIRECTORY="$2"
MIN_COVERAGE="$3"
MIN_VAR_FREQ="$4"
P_VALUE="$5"
GET_MUTECT="$6"
GET_VARSCAN="$7"
GET_HAPLOTYPE="$8"
USE_BAM="$9"
INTERVALS_FILE="${10}"
NORMAL_SAMPLE_FILENAME="${11}"
CODE_DIRECTORY="${12}"
PARAMETER_FILE=${13}
BWA_GREF_FILENAME=${14}
TWIST_SNPS=${15}
ASSEMBLY=${16}
FUNCOTATOR_SOURCES_FILENAME=${17}    
TRANSCRIPT_LIST_FILENAME=${18}
FILTERED=${19}
RUN_FUNCOTATOR=${20}
BAM_OUT=${21}
NORMAL_PILEUPS_FILENAME="${22}"
CALCULATE_CONTAMINATION=${23}

echo "PARENT_DIRECTORY=$1 \
    OUTPUT_DIRECTORY=$2 \
    MIN_COVERAGE=$3 \
    MIN_VAR_FREQ=$4 \
    P_VALUE=$5 \
    GET_MUTECT=$6 \
    GET_VARSCAN=$7 \
    GET_HAPLOTYPE=$8 \
    USE_BAM=$9 \
    INTERVALS_FILE=${10} \
    NORMAL_SAMPLE_FILENAME=${11} \
    CODE_DIRECTORY=${12} \
    PARAMETER_FILE=${13} \
    BWA_GREF_FILENAME=${14} \
    TWIST_SNPS=${15} \
    ASSEMBLY=${16} \
    FUNCOTATOR_SOURCES_FILENAME=${17} \ 
    TRANSCRIPT_LIST_FILENAME=${18} \
    FILTERED=${19} \
    RUN_FUNCOTATOR=${20} \
    BAM_OUT=${21}
    NORMAL_PILEUPS_FILENAME=${22}
    CALCULATE_CONTAMINATION=${23}" >> $PARAMETER_FILE


#NORMAL_NAME=$(basename $NORMAL_SAMPLE | sed -e 's/.bam//')

LINE_NUMBER=$SLURM_ARRAY_TASK_ID #get index of which file to process from $SLURM_ARRAY_TASK_ID provided by SLURM

if [ $USE_BAM = false ]; then
    ARRAY_FILE="${PARENT_DIRECTORY}/fastq_files" #provide path to file containing list of fastq files
else 
    ARRAY_FILE="${PARENT_DIRECTORY}/bam_files" #provide path to file containing list of fastq files
fi

ARRAY_PREFIX="$(sed "${LINE_NUMBER}q; d" "${ARRAY_FILE}")" #extract only the line number corresponding to $SLURM_ARRAY_TASK_ID

CHR_INTERVALS="/oak/stanford/groups/smontgom/maurertm/workflows/Interval_filtering/whole_genome_intervals.interval_list"
#CHR_INTERVALS="/oak/stanford/groups/smontgom/maurertm/workflows/Interval_filtering/chr21_chr22.interval_list"

#Run script and save files in the same location of the fastq files
cd "${OUTPUT_DIRECTORY}" || exit

FILENAME=$(basename "${ARRAY_PREFIX}")
SAMPLE_NAME="${FILENAME}_${ASSEMBLY}"
if [ $USE_BAM = false ]; then
   R1="${ARRAY_PREFIX}_R1_001.fastq.gz"
   R2="${ARRAY_PREFIX}_R2_001.fastq.gz"
   READGROUP="@RG\tID:${FILENAME}\tLB:${FILENAME}\tPL:illumina\tSM:${FILENAME}"
   echo "$READGROUP"
   MUTECT_INPUT="${FILENAME}_${ASSEMBLY}"

   echo "FILENAME: $FILENAME"

   ${CODE_DIRECTORY}/fastq_to_bam.sh $SAMPLE_NAME $READGROUP $BWA_GREF_FILENAME $R1 $R2 $PARAMETER_FILE
else
    MUTECT_INPUT_FILENAME="${ARRAY_PREFIX}"
fi

temp_dir="/tmp/sjaiswal"
INPUT="${TMPDIR}/Inputs"
OUTPUTS="${TMPDIR}/Outputs"
PARAMS="${temp_dir}/Params"


mkdir -p $INPUT
mkdir -p $OUTPUTS
mkdir -p $PARAMS


echo "temp_path is: " $temp_dir
echo "copying FASTQs from the data path..."                                                     #copy data from data_path to the temp_path

MUTECT_INPUT_DIRNAME=$(dirname "${MUTECT_INPUT_FILENAME}")
MUTECT_INPUT_BASENAME=$(basename "${MUTECT_INPUT_FILENAME}")
MUTECT_INPUT="${INPUT}/${MUTECT_INPUT_BASENAME}"
echo "MUTECT_INPUT_BASENAME: $MUTECT_INPUT_BASENAME
      MUTECT_INPUT: $MUTECT_INPUT" >> $PARAMETER_FILE
rsync -vurPhlt "${MUTECT_INPUT_FILENAME}.bam" "${MUTECT_INPUT}.bam"
MUTECT_INPUT_PATTERN=${MUTECT_INPUT_BASENAME%.*}
rsync -vurPhlt $MUTECT_INPUT_DIRNAME/$MUTECT_INPUT_PATTERN* $INPUT/

BWA_GREF_DIRNAME=$(dirname "${BWA_GREF_FILENAME}")
BWA_GREF_BASENAME=$(basename "${BWA_GREF_FILENAME}")
BWA_GREF="${PARAMS}/${BWA_GREF_BASENAME}"
echo "BWA_GREF_BASENAME: $BWA_GREF_BASENAME
      BWA_GREF: $BWA_GREF" >> $PARAMETER_FILE
rsync -vurPhlt "$BWA_GREF_FILENAME" $BWA_GREF
BWA_GREF_PATTERN=${BWA_GREF_BASENAME%.*}
rsync -vurPhlt $BWA_GREF_DIRNAME/$BWA_GREF_PATTERN* $PARAMS/



##################################################################################################################################
#################################################---STEP 2: MUTECT.sh---########################################################## 
##################################################################################################################################       
if [ $GET_MUTECT = true ]; then
    echo "Mutect analysis requested"
    echo "$PARAMETER_FILE"

    NORMAL="${temp_dir}/Normal"
    
    mkdir -p $NORMAL
    
    if [ $RUN_FUNCOTATOR = true ]; then
        FUNCTATOR_SOURCES_BASENAME=$(basename $FUNCOTATOR_SOURCES_FILENAME)
        rsync -vurPhlt $FUNCOTATOR_SOURCES_FILENAME $PARAMS
        FUNCOTATOR_SOURCES=$PARAMS/$FUNCOTATOR_SOURCES_BASENAME
    else 
        FUNCOTATOR_SOURCES=false
    fi
    
    if [ $NORMAL_SAMPLE_FILENAME != false ]; then 
        NORMAL_SAMPLE_DIRNAME=$(dirname "${NORMAL_SAMPLE_FILENAME}")
        NORMAL_SAMPLE_BASENAME=$(basename "${NORMAL_SAMPLE_FILENAME}")
        NORMAL_SAMPLE="${NORMAL}/${NORMAL_SAMPLE_BASENAME}"
        echo "NORMAL_SAMPLE_BASENAME: $NORMAL_SAMPLE_BASENAME
              NORMAL SAMPLE: $NORMAL_SAMPLE" >> $PARAMETER_FILE
        rsync -vurPhlt "$NORMAL_SAMPLE_FILENAME" $NORMAL_SAMPLE
        NORMAL_SAMPLE_PATTERN=${NORMAL_SAMPLE_BASENAME%.*}
        rsync -vurPhlt $NORMAL_SAMPLE_DIRNAME/$NORMAL_SAMPLE_PATTERN* $NORMAL/
        
        if [ $CALCULATE_CONTAMINATION != false ]; then
            NORMAL_PILEUPS_DIRNAME=$(dirname "${NORMAL_PILEUPS_FILENAME}")
            NORMAL_PILEUPS_BASENAME=$(basename "${NORMAL_PILEUPS_FILENAME}")
            NORMAL_PILEUPS="${NORMAL}/${NORMAL_PILEUPS_BASENAME}"
            echo "NORMAL_PILEUPS_BASENAME: $NORMAL_PILEUPS_BASENAME
                  NORMAL PILEUPS: $NORMAL_PILEUPS" >> $PARAMETER_FILE
            rsync -vurPhlt "$NORMAL_PILEUPS_FILENAME" $NORMAL_PILEUPS
            NORMAL_PILEUPS_PATTERN=${NORMAL_PILEUPS_BASENAME%.*}
            rsync -vurPhlt $NORMAL_PILEUPS_DIRNAME/$NORMAL_PILEUPS_PATTERN* $NORMAL/
        else
            NORMAL_PILEUPS=false
        fi

    else
        NORMAL_SAMPLE=false
        NORMAL_PILEUPS=false
    fi
    
    TRANSCRIPT_LIST_BASENAME=$(basename "${TRANSCRIPT_LIST_FILENAME}")
    TRANSCRIPT_LIST="${PARAMS}/${TRANSCRIPT_LIST_BASENAME}"
    echo "TRANSCRIPT_LIST_BASENAME: $TRANSCRIPT_LIST_BASENAME
          TRANSCRIPT_LIST: $TRANSCRIPT_LIST" >> $PARAMETER_FILE
    rsync -vurPhlt "$TRANSCRIPT_LIST_FILENAME" $TRANSCRIPT_LIST

    ls $OUTPUTS

    if [ $INTERVALS_FILE = false ]; then
        num_intervals=$(grep -v "@" < $CHR_INTERVALS | wc -l)
        seq 1 ${num_intervals} | parallel -j8 --progress --ungroup "${CODE_DIRECTORY}/mutect_and_pileups.sh  $NORMAL_SAMPLE $INTERVALS_FILE $MUTECT_INPUT $SAMPLE_NAME $BWA_GREF $PARAMETER_FILE  $OUTPUTS $BAM_OUT $OUTPUT_DIRECTORY $CALCULATE_CONTAMINATION {}"

        echo "arguments passed to mutect_and_pileups.sh: seq 1 ${num_intervals} | parallel -j8 --group --progress ${CODE_DIRECTORY}/mutect_and_pileups.sh  $NORMAL_SAMPLE $INTERVALS_FILE $MUTECT_INPUT $SAMPLE_NAME $BWA_GREF $PARAMETER_FILE  $OUTPUTS $BAM_OUT $OUTPUT_DIRECTORY $CALCULATE_CONTAMINATION {}">> $PARAMETER_FILE
    
    else
        ${CODE_DIRECTORY}/mutect_and_pileups.sh  $NORMAL_SAMPLE $INTERVALS_FILE $MUTECT_INPUT $SAMPLE_NAME $BWA_GRE  $PARAMETER_FILE  $OUTPUTS $BAM_OUT $OUTPUT_DIRECTORY $CALCULATE_CONTAMINATION
        echo "arguments passed to mutect_and_pileups.sh: ${CODE_DIRECTORY}/mutect_and_pileups.sh $NORMAL_SAMPLE $INTERVALS_FILE $MUTECT_INPUT $SAMPLE_NAME $BWA_GRE  $PARAMETER_FILE  $OUTPUTS $BAM_OUT $OUTPUT_DIRECTORY $CALCULATE_CONTAMINATION" >> $PARAMETER_FILE
    
    fi
   
   echo "mutect_and_pileups.sh analysis complete"
  
   if [ $INTERVALS_FILE = false ]; then
       module load bcftools
       module load gatk4
       chr_vcfs=$(find $OUTPUT_DIRECTORY/vcfs -type f | sort -V)

       #if [[ $(echo "${chr_vcfs}" | wc -l) -lt ${num_intervals} ]] && [[ -d "${OUTPUT_DIRECTORY}/vcfs" ]]; then
       if [ -d "${OUTPUT_DIRECTORY}/vcfs" ]; then
          rsync -vurPhlt "${OUTPUT_DIRECTORY}/vcfs" "${OUTPUTS}"
       fi
       chr_pileups=$(find $OUTPUT_DIRECTORY/pileups -type f | sort -V)
       #if [[ $(echo "${chr_pileups}" | wc -l) -lt ${num_intervals} ]] && [[ -d "${OUTPUT_DIRECTORY}/pileups" ]]; then
       if [ -d "${OUTPUT_DIRECTORY}/pileups" ]; then
            rsync -vurPhlt "${OUTPUT_DIRECTORY}/pileups" "${OUTPUTS}"
       fi
       # Change to MergeVCFs
       echo "bcftools concat $(find $OUTPUTS/vcfs -type f | grep -E ".*.vcf$" | sort -V) > ${OUTPUTS}/${SAMPLE_NAME}_mutect2.vcf"
       bcftools concat $(find $OUTPUTS/vcfs -type f | grep -E ".*.vcf$" | sort -V) > ${OUTPUTS}/${SAMPLE_NAME}_mutect2.vcf

       vcf_stats=$(find $OUTPUTS/vcfs -type f | grep -E ".*.vcf.stats$" | sed -e 's/^/--stats /g' | tr '\n' ' ')
       gatk MergeMutectStats ${vcf_stats} -O "${OUTPUTS}/${SAMPLE_NAME}_mutect2.vcf.stats"

       echo "gatk MergeMutectStats ${vcf_stats} -O ${OUTPUTS}/${SAMPLE_NAME}_mutect2.vcf.stats"
       # Add GatherPileupSummaries for pileups
       pileup_tables=$(find $OUTPUTS/pileups -type f | sed -e 's/^/-I /g' | tr '\n' ' ')
       echo "gatk GatherPileupSummaries --sequence-dictionary "/oak/stanford/groups/sjaiswal/Herra/CHIP_Panel_AmpliSeq/GRCh38.p12.genome.u2af1l5_mask.fa.dict" ${pileup_tables} -O ${OUTPUTS}/${SAMPLE_NAME}_pileups.table"

       gatk GatherPileupSummaries --sequence-dictionary "/oak/stanford/groups/sjaiswal/Herra/CHIP_Panel_AmpliSeq/GRCh38.p12.genome.u2af1l5_mask.fa.dict" ${pileup_tables} -O ${OUTPUTS}/${SAMPLE_NAME}_pileups.table
   fi
    echo "$OUTPUTS folder contains:"
    tree -h $OUTPUTS

    ${CODE_DIRECTORY}/mutect_filter.sh $SAMPLE_NAME $BWA_GREF $PARAMETER_FILE $OUTPUTS $OUTPUT_DIRECTORY $NORMAL_SAMPLE $CALCULATE_CONTAMINATION $NORMAL_PILEUPS

    ${CODE_DIRECTORY}/funcotator.sh $SAMPLE_NAME $BWA_GREF $FUNCOTATOR_SOURCES $TRANSCRIPT_LIST $PARAMETER_FILE $FILTERED $OUTPUTS $RUN_FUNCOTATOR $OUTPUT_DIRECTORY

   rsync -vurPhlt $OUTPUTS/ $OUTPUT_DIRECTORY


   echo "Results copied over to specified output directory"
else
    echo "No mutect analysis requested"
fi
            
##################################################################################################################################
############################################---STEP 3: HAPLOTYPECALLER.sh---######################################################
##################################################################################################################################    
if [ $GET_HAPLOTYPE = true ]; then
    echo "Haplotypecaller analysis requested"
    ${CODE_DIRECTORY}/haplotypecaller.sh $SAMPLE_NAME $BWA_GREF $TWIST_SNPS $PARAMETER_FILE
    echo "Haplotypecaller analysis complete"
else
    echo "No HaplotypeCaller analysis requested"
fi

##################################################################################################################################
################################################---STEP 4: VARSCAN.sh---########################################################## 
##################################################################################################################################    
if [ $GET_VARSCAN = true ]; then
    echo "Varscan analysis requested"
   ${CODE_DIRECTORY}/varscan.sh $SAMPLE_NAME $BWA_GREF $MIN_COVERAGE $MIN_VAR_FREQ $P_VALUE $PARAMETER_FILE
   echo "Varscan analysis complete"
else 
    echo "No Varscan analysis requested"
fi
