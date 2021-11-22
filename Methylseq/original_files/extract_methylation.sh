#!/bin/bash

#SBATCH --job-name=extract_methyl
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=14
#SBATCH --partition=interactive
#SBATCH --mem=128G
#SBATCH --time=23:00:00

module load bismark/0.22.3
#module load java
#module load R
#module load picard/2.9.5

#module load bowtie2/2.4.4

seq_path=$1/fastq
data_path=$1
unmethyl_control=$2
hydroxymethyl_control=$3
genome_path=$4
temp_path=$5
cores=$6

temp_genomes_path=/tmp/genomes
rsync -vur $genome_path $temp_genomes_path
genome_name=$(basename $genome_path)


#cd $data_path/fastq/$unmethyl_control/$hydroxymethyl_control/genome_alignment/
cd $temp_path/$unmethyl_control/$hydroxymethyl_control/genome_alignment/ 
for sample in *deduplicated_sorted.bam
do
    bismark_methylation_extractor --buffer_size 20G --gzip --cytosine_report --bedGraph --genome_folder $temp_genomes_path/$genome_name $sample --multicore $cores
done

#copy files back to seq_path directory
rsync -vur $temp_path/ $seq_path



