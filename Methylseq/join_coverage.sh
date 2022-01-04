#!/bin/bash

#SBATCH --time=1:00:00
#SBATCH --account=sjaiswal
#SBATCH --cpus-per-task=8
#SBATCH --mem=32GB
#SBATCH --job-name=methylseq

data_path=$1
seq_path="$1/fastq"
unmethyl_control=$2
hydroxymethyl_control=$3
cores=$4

#temp_path=$(mktemp -d /tmp/tmp.XXXXXXXXXX)
#echo "temp_path is: " $temp_path
#mkdir $temp_path
#copy fastq files to temp_path
#rsync -vur "$data_path/fastq/" $temp_path

line_number=$SLURM_ARRAY_TASK_ID #get index of which file to process from $SLURM_ARRAY_TASK_ID provided by SLURM
fastq_file="${data_path}/fastq/FASTQs" #provide path to file containing list of fastq files
fastq_path="$(sed "${line_number}q; d" $fastq_file)" #extract only the line number corresponding to $SLURM_ARRAY_TASK_ID
sample_name="${fastq_path##*/}"
echo "sample name: $sample_name"
output_bed="${sample_name}_R1_001.trimmed.fastq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1_bismark_bt2_pe.bam.bedGraph"
output_cov="${sample_name}_R1_001.trimmed.fastq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1_bismark_bt2_pe.bam.bismark.cov"

cd "$data_path/fastq/$unmethyl_control/$hydroxymethyl_control/genome_alignment" 
find split_bams -name "*${sample_name}.bedGraph.gz" | xargs -I % bash -c 'zcat % >>'"$output_bed" 

find split_bams -name "*${sample_name}.bismark.cov.gz" | xargs -I % bash -c 'zcat % >>'"$output_cov"
#rsync -vur $temp_path/ "$data_path/fastq/" 

gzip ${sample_name}_R1_001.trimmed.fastq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1_bismark_bt2_pe.bam.bismark.cov.gz
gzip ${sample_name}_R1_001.trimmed.fastq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1_bismark_bt2_pe.bam.bedGraph.gz 
