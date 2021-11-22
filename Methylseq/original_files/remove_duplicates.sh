#!/bin/bash

#SBATCH --job-name=deduplicate
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
unmethyl_control=$2
hydroxymethyl_control=$3
temp_path=$4

#cd /oak/stanford/groups/sjaiswal/kameronr/JG97/methylseq/miseq/213566405/FASTQ/
#cd /oak/stanford/groups/sjaiswal/kameronr/sjaiswal_old/bisulfiteseq/bisulfite_run180925/

#for read1 in *R1_001.trimmed.fastq.gz
#do 
#    read2=$(echo $read1| sed 's/R1/R2/')
   
#    bismark --bam --maxins 800 /oak/stanford/groups/sjaiswal/kameronr/JG97/methylseq/miseq/fasta_unmethylated -1 $read1 -2 $read2 --parallel 8 -o ./unmethyl_control --unmapped 

    #bismark --bam --maxins 800 /oak/stanford/groups/sjaiswal/kameronr/JG97/methylseq/miseq/fasta_oxBS -1 $read1 -2 $read2 --parallel 8 -o ./oxBS_control --unmapped

#done

#cd /oak/stanford/groups/sjaiswal/kameronr/JG97/methylseq/miseq/213566405/FASTQ/unmethyl_control/alignment

#cd $seq_path/fastq/$unmethyl_control/$hydroxymethyl_control/genome_alignment
cd $temp_path/$unmethyl_control/$hydroxymethyl_control/genome_alignment

#cd /oak/stanford/groups/sjaiswal/kameronr/sjaiswal_old/bisulfiteseq/bisulfite_run180925/

for sample in *bismark_bt2_pe.bam
do
    #java -Xmx4g -XX:+UseG1GC -XX:ParallelGCThreads=8 -jar ~/picard.jar CollectInsertSizeMetrics INPUT=$sample OUTPUT={$sample}_picard_insert_size_metrics.txt HISTOGRAM_FILE={$sample}_picard_insert_size_plot.pdf METRIC_ACCUMULATION_LEVEL=ALL_READS
    #picard CollectInsertSizeMetrics INPUT=$sample OUTPUT={$sample}_picard_insert_size_metrics.txt HISTOGRAM_FILE={$sample}_picard_insert_size_plot.pdf METRIC_ACCUMULATION_LEVEL=ALL_READS
    deduplicate_bismark -p --bam $sample
done

#copy files back to seq_path directory
rsync -vur $temp_path/ $seq_path

