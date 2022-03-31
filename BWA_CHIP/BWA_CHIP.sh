#!/bin/bash

#SBATCH --time=3:00:00
#SBATCH --account=sjaiswal
#SBATCH --cpus-per-task=8
#SBATCH --mem=32GB
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
TUMOR_NORMAL="${11}"
NORMAL_SAMPLE="${12}"
CODE_DIRECTORY="${13}"
PARAMETER_FILE=${14}
BWA_GREF=${15}
TWIST_SNPS=${16}
ASSEMBLY=${17}
FUNCOTATOR_SOURCES=${18}    
TRANSCRIPT_LIST=${19}
FILTERED=${20}
NORMAL_NAME=${21}
echo "normal name: $NORMAL_NAME"

#NORMAL_NAME=$(basename $NORMAL_SAMPLE | sed -e 's/.bam//')

LINE_NUMBER=$SLURM_ARRAY_TASK_ID #get index of which file to process from $SLURM_ARRAY_TASK_ID provided by SLURM

if [ $USE_BAM == false ]; then
    ARRAY_FILE="${PARENT_DIRECTORY}/fastq_files" #provide path to file containing list of fastq files
else 
    ARRAY_FILE="${PARENT_DIRECTORY}/bam_files" #provide path to file containing list of fastq files
fi

ARRAY_PREFIX="$(sed "${LINE_NUMBER}q; d" "${ARRAY_FILE}")" #extract only the line number corresponding to $SLURM_ARRAY_TASK_ID
FILENAME=$(basename "${ARRAY_PREFIX}")
PREFIX=$FILENAME

echo "FILENAME: $FILENAME"
echo "PREFIX: $PREFIX"

R1="${ARRAY_PREFIX}_R1_001.fastq.gz"
R2="${ARRAY_PREFIX}_R2_001.fastq.gz"
READGROUP="@RG\tID:${PREFIX}\tLB:${PREFIX}\tPL:illumina\tSM:${PREFIX}"

#Run script and save files in the same location of the fastq files
cd "${OUTPUT_DIRECTORY}" || exit

SAMPLE_NAME="${PREFIX}_${ASSEMBLY}"
if [ $USE_BAM == false ]; then
   MUTECT_INPUT="${PREFIX}_${ASSEMBLY}"
   ${CODE_DIRECTORY}/fastq_to_bam.sh $SAMPLE_NAME $READGROUP $BWA_GREF $R1 $R2 $PARAMETER_FILE
else
    MUTECT_INPUT="${ARRAY_PREFIX}"
fi

##################################################################################################################################
#################################################---STEP 2: MUTECT.sh---########################################################## 
##################################################################################################################################       
if [ $GET_MUTECT = true ]; then
    echo "Mutect analysis requested"
    echo "$PARAMETER_FILE"

   ${CODE_DIRECTORY}/mutect.sh $TUMOR_NORMAL $NORMAL_SAMPLE $NORMAL_NAME $INTERVALS_FILE $MUTECT_INPUT $SAMPLE_NAME $BWA_GREF $FUNCOTATOR_SOURCES $TRANSCRIPT_LIST $PARAMETER_FILE $FILTERED
   echo "Mutect analysis complete"
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
