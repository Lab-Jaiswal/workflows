#!/bin/bash

#SBATCH --time=3:00:00
#SBATCH --account=sjaiswal
#SBATCH --cpus-per-task=8
#SBATCH --mem=32GB
#SBATCH --job-name=CHIP_variant_call

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
PAIRED="${11}"
NORMAL_SAMPLE="${12}"
CODE_DIRECTORY="${13}"
echo "intervals file : $INTERVALS_FILE"
echo "parameters used:
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
PAIRED="${11}"
NORMAL_SAMPLE="${12}"
CODE_DIRECTORY="${13}"
"

NORMAL_NAME=$(basename $NORMAL_SAMPLE | sed -e 's/.bam//')

line_number=$SLURM_ARRAY_TASK_ID #get index of which file to process from $SLURM_ARRAY_TASK_ID provided by SLURM
echo "use_bam: $USE_BAM"
if [ $USE_BAM == false ]; then
    array_file="${PARENT_DIRECTORY}/fastq_files" #provide path to file containing list of fastq files
else 
    array_file="${PARENT_DIRECTORY}/bam_files" #provide path to file containing list of fastq files
fi
array_prefix="$(sed "${line_number}q; d" "${array_file}")" #extract only the line number corresponding to $SLURM_ARRAY_TASK_ID

echo "array prefix:i $array_prefix"
FILENAME=$(basename "${array_prefix}")
PREFIX=$FILENAME

echo "FILENAME: $FILENAME"
echo "PREFIX: $PREFIX"

R1="${array_prefix}_R1_001.fastq.gz"
R2="${array_prefix}_R2_001.fastq.gz"
READGROUP="@RG\tID:${PREFIX}\tLB:${PREFIX}\tPL:illumina\tSM:${PREFIX}"
#BWA_GREF="/labs/sjaiswal/genomes/GRCh38/GRCh38.p13.genome.fa"
BWA_GREF="/oak/stanford/groups/sjaiswal/Herra/CHIP_Panel_AmpliSeq/GRCh38.p12.genome.u2af1l5_mask.fa" #reference genome
TWIST_SNPS="/labs/sjaiswal/workflows/BWA_mutect_twist/twist_snps.bed" #SNPs for germline calling
ASSEMBLY="GRCh38" #Genome version
TRANSCRIPT_LIST="/oak/stanford/groups/sjaiswal/Herra/CHIP_TWIST-PANEL_ATHEROMA/chip_transcript_list.txt" #Transcript list for Mutect
FUNCOTATOR_SOURCES="/labs/sjaiswal/tools/funcotator/funcotator_dataSources.v1.6.20190124s" #Reference for Funcotator

#Run script and save files in the same location of the fastq files
cd "${OUTPUT_DIRECTORY}" || exit

SAMPLE_NAME="${PREFIX}_${ASSEMBLY}"
if [ $USE_BAM == false ]; then
   MUTECT_INPUT="${PREFIX}_${ASSEMBLY}"
   ${CODE_DIRECTORY}/fastq_to_bam.sh $SAMPLE_NAME $READGROUP $BWA_GREF $R1 $R2
else
    MUTECT_INPUT="${array_prefix}"
fi

        
if [ $GET_MUTECT = true ]; then
   ${CODE_DIRECTORY}/mutect.sh $PAIRED $NORMAL_SAMPLE $NORMAL_NAME $INTERVALS_FILE $MUTECT_INPUT $SAMPLE_NAME $BWA_GREF $FUNCOTATOR_SOURCES $TRANSCRIPT_LIST
else
    echo "mutect analysis not requested"
fi
            

if [ $GET_HAPLOTYPE = true ]; then
    ${CODE_DIRECTORY}/haplotypecaller.sh $SAMPLE_NAME $BWA_GREF $TWIST_SNPS
else
    echo "No HaplotypeCaller analysis requested"
fi

if [ $GET_VARSCAN = true ]; then
   ${CODE_DIRECTORY}/varscan.sh $SAMPLE_NAME $BWA_GREF $MIN_COVERAGE $MIN_VAR_FREQ $P_VALUE
else 
    echo "No Varscan analysis requested"
fi
