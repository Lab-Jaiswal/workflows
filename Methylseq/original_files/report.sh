#!/bin/bash

#SBATCH --job-name=extract_methyl
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=14
#SBATCH --partition=interactive
#SBATCH --mem=128G
#SBATCH --time=23:00:00

module load bismark/0.22.3

#cd /oak/stanford/groups/sjaiswal/kameronr/JG97/methylseq/miseq/213566405/FASTQ/
#cd /oak/stanford/groups/sjaiswal/kameronr/sjaiswal_old/bisulfiteseq/bisulfite_run180925/

seq_path=$1/fastq
unmethyl_control=$2
hydroxymethyl_control=$3
temp_path=$4

#http://felixkrueger.github.io/Bismark/Docs/
#bismark report options:
#--alignment_report FILE
#--dedup_report FILE
#--splitting_report FILE
#--mbias_report FILE
#--nucleotide_report FILE



#https://rawgit.com/FelixKrueger/Bismark/master/Docs/Bismark_User_Guide.html


#cd $seq_path/fastq/$unmethyl_control/$hydroxymethyl_control/genome_alignment/

cd $temp_path/$unmethyl_control/$hydroxymethyl_control/genome_alignment

bismark2report
bismark2summary
#need to specify a nucleotide coverage report file in the above command! ^

#cd /oak/stanford/groups/sjaiswal/kameronr/sjaiswal_old/bisulfiteseq/bisulfite_run180925/

#for sample in *.trimmed_bismark_bt2_pe.bam
#do
#    bismark_methylation_extractor --gzip --cytosine_report --bedGraph --genome_folder /home/kameronr/sjaiswal/kameronr/JG97/methylseq/miseq/fasta_unmethylated $sample
#done

#copy files back to seq_path directory
rsync -vur $temp_path/ $seq_path



