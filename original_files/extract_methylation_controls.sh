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



#cd /oak/stanford/groups/sjaiswal/kameronr/JG97/methylseq/miseq/213566405/FASTQ/
#cd /oak/stanford/groups/sjaiswal/kameronr/sjaiswal_old/bisulfiteseq/bisulfite_run180925/

#for read1 in *R1_001.trimmed.fastq.gz
#do 
#    read2=$(echo $read1| sed 's/R1/R2/')
   
#    bismark --bam --maxins 800 /oak/stanford/groups/sjaiswal/kameronr/JG97/methylseq/miseq/fasta_unmethylated -1 $read1 -2 $read2 --parallel 8 -o ./unmethyl_control --unmapped 

    #bismark --bam --maxins 800 /oak/stanford/groups/sjaiswal/kameronr/JG97/methylseq/miseq/fasta_oxBS -1 $read1 -2 $read2 --parallel 8 -o ./oxBS_control --unmapped

#done
seq_path=$1/fastq
unmethyl_control=$2
unmethyl_control_fasta=$3
hydroxymethyl_control=$4
hydroxymethyl_control_fasta=$5
temp_path=$6
cores=$7

temp_genomes_path=/tmp/genomes
rsync -vur /oak/stanford/groups/smontgom/kameronr/methylseq/JG97/$unmethyl_control_fasta $temp_genomes_path
rsync -vur /oak/stanford/groups/smontgom/kameronr/methylseq/JG97/$hydroxymethyl_control_fasta $temp_genomes_path

#cd /oak/stanford/groups/sjaiswal/kameronr/JG97/methylseq/miseq/213566405/FASTQ/unmethyl_control/
#cd $seq_path/fastq/$unmethyl_control
cd $temp_path/$unmethyl_control

#cd /oak/stanford/groups/sjaiswal/kameronr/sjaiswal_old/bisulfiteseq/bisulfite_run180925/

for sample in *.trimmed_bismark_bt2_pe.bam
do
    bismark_methylation_extractor --gzip --cytosine_report --bedGraph --genome_folder $temp_genomes_path/$unmethyl_control_fasta $sample --multicore $cores
done


cd $temp_path/$unmethyl_control/$hydroxymethyl_control

for sample in *.trimmed.fastq.gz_unmapped_reads_1_bismark_bt2_pe.bam
do
    bismark_methylation_extractor --gzip --cytosine_report --bedGraph --genome_folder $temp_genomes_path/$hydroxymethyl_control_fasta $sample --multicore $cores
done

#sync temp directory with main directory
rsync -vur $temp_path/ $seq_path

