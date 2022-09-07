#!/bin/bash

#SBATCH --job-name=control_map
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=14
#SBATCH --partition=interactive
#SBATCH --mem=128G
#SBATCH --time=24:00:00

module load bismark/0.22.3

#module load bowtie2/2.4.4

seq_path=$1/fastq
unmethyl_control=$2
hydroxymethyl_control=$3
genome_path=$4
phix_path=$5
temp_path=$6
cores=$7

temp_genomes_path=/tmp/genomes
genome_name=$(basename $genome_path)
#temp_path=$(mktemp -d /tmp/tmp.XXXXXXXXXX)
#echo "temp_path is: " $temp_path

#copy reference genome files to temp directory
rsync -vur $genome_path $temp_genomes_path

#cd /oak/stanford/groups/sjaiswal/kameronr/JG97/methylseq/miseq/213566405/FASTQ/
#cd /oak/stanford/groups/sjaiswal/kameronr/sjaiswal_old/bisulfiteseq/bisulfite_run180925/

#for read1 in *R1_001.trimmed.fastq.gz
#do 
#    read2=$(echo $read1| sed 's/R1/R2/')
   
#    bismark --bam --maxins 800 /oak/stanford/groups/sjaiswal/kameronr/JG97/methylseq/miseq/fasta_unmethylated -1 $read1 -2 $read2 --parallel 8 -o ./unmethyl_control --unmapped 

    #bismark --bam --maxins 800 /oak/stanford/groups/sjaiswal/kameronr/JG97/methylseq/miseq/fasta_oxBS -1 $read1 -2 $read2 --parallel 8 -o ./oxBS_control --unmapped

#done

#cd /oak/stanford/groups/sjaiswal/kameronr/JG97/methylseq/miseq/213566405/FASTQ/unmethyl_control/
#cd $seq_path/fastq/$unmethyl_control/$hydroxymethyl_control/

#copy reads to temp directory
#rsync -vu $seq_path/fastq/$unmethyl_control/$hydroxymethyl_control/ $temp_path

cd $temp_path/$unmethyl_control/$hydroxymethyl_control

#cd /oak/stanford/groups/sjaiswal/kameronr/sjaiswal_old/bisulfiteseq/bisulfite_run180925/

for read1 in *_R1_001.trimmed_unmapped2.fq.gz
#for read1 in *R1_001.fastq.gz
do
    read2=$(echo $read1| sed 's/R1/R2/')
#    read2=$(echo $read1| sed 's/R1/R2/') 
    bismark --bam --maxins 1000 $temp_genomes_path/$genome_name -1 $read1 -2 $read2 -o $temp_path/$unmethyl_control/$hydroxymethyl_control/genome_alignment --unmapped --nucleotide_coverage --multicore $cores
#--parallel ${SLURM_CPUS_PER_TASK}
done

#copy all results back to #seq_path and updates to the genome
rsync -vur $temp_path/$unmethyl_control/$hydroxymethyl_control/ $seq_path/$unmethyl_control/$hydroxymethyl_control/
rsync -vur $temp_genomes_path/$genome_name $genome_path/..

