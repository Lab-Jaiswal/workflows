#!/bin/bash

#SBATCH --time=1:00:00
#SBATCH --account=sjaiswal
#SBATCH --cpus-per-task=8
#SBATCH --mem=64GB
#SBATCH --job-name=bam_to_fastq
#SBATCH --array=37

bam_location=$1
code_directory=$2

cd $bam_location
bam_file="/labs/sjaiswal/maurertm/bam_files"

line_number=$SLURM_ARRAY_TASK_ID
bam_path="$(sed "${line_number}q; d" "${bam_file}")"

fastq_path_R1="${bam_path}_R1_001.fastq.gz"
fastq_path_R2="${bam_path}_R2_001.fastq.gz"

samtools collate -u -O "${bam_path}.bam" | \\
samtools fastq -1 $fastq_path_R1 -2 $fastq_path_R2 -0 /dev/null -s /dev/null -n

