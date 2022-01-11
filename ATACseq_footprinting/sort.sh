#!/bin/bash

#SBATCH --job-name=sort
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --partition=batch
#SBATCH --mem=500G
#SBATCH --time=3:00:00
#SBATCH --account=sjaiswal
#cd /oak/stanford/groups/sjaiswal/kameronr/JG97/methylseq/miseq/213566405/FASTQ/

seq_path=$1
temp_path=$2

cd $temp_path

module load samtools/1.9



