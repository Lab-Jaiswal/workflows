#!/bin/bash
#SBATCH --job-name=sort
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=batch
#SBATCH --mem=256G
#SBATCH --time=24:00:00
#SBATCH --account=sjaiswal

echo "entering fastq_to_bam.sh"

bam_path=$1
output_path=$2

line_number=$SLURM_ARRAY_TASK_ID #get index of which file to process from $SLURM_ARRAY_TASK_ID provided by SLURM
bam_file="${bam_path}/bam_files" #provide path to file containing list of fastq files
bam_prefix="$(sed "${line_number}q; d" "${bam_file}")" #extract only the line number corresponding to $SLURM_ARRAY_TASK_ID

#Collect the sample name from the BAMs file
#Ex: If $bam_prefix is "/atac_seq/data/ATAC_tet2_KO_LDL", then the PREFIX is "ATAC_tet2_KO_LDL"
SAMPLE_NAME=$(basename "${bam_prefix}")
echo "SAMPLE_NAME: $SAMPLE_NAME"

fastq_path_R1="$output_path/${SAMPLE_NAME}_R1.fq.gz"
fastq_path_R2="$output_path/${SAMPLE_NAME}_R2.fq.gz"

echo "bam_to_fastq command used the following parameters: 
$0 $1 $2"

echo "variables used:
SAMPLE_NAME: $SAMPLE_NAME
R1: $fastq_path_R1
R2: $fastq_path_R2"

if [ ! -f "$fastq_path_R1" ]; then
    module load samtools
    samtools collate -u -O "${SAMPLE_NAME}.bam" | \
    samtools fastq -1 $fastq_path_R1 -2 $fastq_path_R2 -0 /dev/null -s /dev/null -n
    echo "Conversion of bam to fastq file complete"
else
    echo "Conversion of bam to fastq file already completed"
fi

echo "bam_to_fastq.sh complete"
