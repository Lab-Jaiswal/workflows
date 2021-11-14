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

#cd /oak/stanford/groups/sjaiswal/kameronr/JG97/methylseq/miseq/213566405/FASTQ/unmethyl_control/

seq_path=$1
unmethyl_control=$2
hydroxymethyl_control=$3

cd $seq_path/fastq/$unmethyl_control/

#https://rawgit.com/FelixKrueger/Bismark/master/Docs/Bismark_User_Guide.html
bismark2report

cd $seq_path/fastq/$unmethyl_control/$hydroxymethyl_control

bismark2report

#cd /oak/stanford/groups/sjaiswal/kameronr/sjaiswal_old/bisulfiteseq/bisulfite_run180925/

#for sample in *.trimmed_bismark_bt2_pe.bam
#do
#    bismark_methylation_extractor --gzip --cytosine_report --bedGraph --genome_folder /home/kameronr/sjaiswal/kameronr/JG97/methylseq/miseq/fasta_unmethylated $sample
#done



