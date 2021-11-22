#!/bin/bash

#SBATCH --job-name=methylseq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --partition=batch
#SBATCH --mem=500G
#SBATCH --time=3:00:00
#SBATCH --account=sjaiswal
#cd /oak/stanford/groups/sjaiswal/kameronr/JG97/methylseq/miseq/213566405/FASTQ/

seq_path=$1/fastq
temp_path=$2

cd $temp_path

#module load fastq
#fastqc *.fq
module load cutadapt/3.4

#should I do any filtering of reads or trimming reads based on their quality? yes. I probably should but I don't know by how much or by what thresholds! For now I'll just do -q 20

#keeps reads at least 15 or more base pairs long

#copy files to temp directory
#temp_path=$(mktemp -d /tmp/tmp.XXXXXXXXXX)
#echo "temp_path is: " $temp_path
#mkdir $temp_path
#rsync -vur $seq_path/ $temp_path

for read1 in *R1_001.fastq.gz
#for read1 in *_1.fq.gz
do 
    read2=$(echo $read1| sed 's/R1_001.fastq.gz/R2_001.fastq.gz/')
    read1_output=$(echo $read1| sed 's/fastq.gz/trimmed.fastq.gz/')
    echo $read1_output
    
    read2_output=$(echo $read2| sed 's/fastq.gz/trimmed.fastq.gz/')
    echo $read2_output

    cutadapt -q 20 -m 15 -a AGATCGGAAGAGC -A AAATCAAAAAAAC -o $read1_output -p $read2_output $read1 $read2 --cores=${SLURM_CPUS_PER_TASK}
done

#copy files back to seq_path directory
rsync -vur $temp_path/ $seq_path

