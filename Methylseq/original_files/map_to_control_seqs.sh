#!/bin/bash

#SBATCH --job-name=control_map
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=14
#SBATCH --partition=interactive
#SBATCH --mem=32G
#SBATCH --time=01:00:00

module load bismark/0.22.3
#module load bowtie2/2.4.4

#cd /oak/stanford/groups/sjaiswal/kameronr/JG97/methylseq/miseq/213566405/FASTQ/
#cd /oak/stanford/groups/sjaiswal/kameronr/sjaiswal_old/bisulfiteseq/bisulfite_run180925/

data_path=$1
seq_path=$1/fastq
unmethyl_control=$2
unmethyl_control_fasta=$3
hydroxymethyl_control=$4
hydroxymethyl_control_fasta=$5
temp_path=$6
cores=$7

temp_genomes_path=/tmp/genomes
#temp_path=$(mktemp -d /tmp/tmp.XXXXXXXXXX)
#echo "temp_path is: " $temp_path
#mkdir $temp_path

#copy reference genome files to temp directory
#mkdir $temp_genomes_path/$unmethyl_control_fasta
#mkdir $temp_genomes_path/$hydroxymethyl_control_fasta
#cp /oak/stanford/groups/sjaiswal/kameronr/JG97/methylseq/miseq/$unmethyl_control_fasta $temp_path -R
#cp /oak/stanford/groups/sjaiswal/kameronr/JG97/methylseq/miseq/$hydroxymethyl_control_fasta $temp_path -R
rsync -vur /oak/stanford/groups/smontgom/kameronr/methylseq/JG97/$unmethyl_control_fasta $temp_genomes_path
rsync -vur /oak/stanford/groups/smontgom/kameronr/methylseq/JG97/$hydroxymethyl_control_fasta $temp_genomes_path

#copy reads to temp
#rsync -vur $seq_path/ $temp_path

cd $temp_path

for read1 in *R1_001.trimmed.fastq.gz
do 
    read2=$(echo $read1| sed 's/R1/R2/')   
    file_name=$(echo $read1| sed 's/.fastq.gz//')
    file_name_check=$file_name\_bismark_bt2_PE_report.txt
    if [ ! -f $temp_path/$unmethyl_control/$file_name_check ]; then
	echo $file_name_check does not exist yet
	echo Run bismark on $read1 and $read2
	bismark --bam --maxins 800 $temp_genomes_path/$unmethyl_control_fasta -1 $read1 -2 $read2 -o $temp_path/$unmethyl_control --unmapped --nucleotide_coverage --multicore $cores
    else
	echo Skipping Bismark mapping for $read1
    fi
    
done
#--p ${SLURM_CPUS_PER_TASK}

#cd /oak/stanford/groups/sjaiswal/kameronr/JG97/methylseq/miseq/213566405/FASTQ/unmethyl_control/
#unmapped_path=$seq_path/$unmethyl_control/

unmapped_path=$temp_path/$unmethyl_control
cd $unmapped_path
#cd /oak/stanford/groups/sjaiswal/kameronr/sjaiswal_old/bisulfiteseq/bisulfite_run180925/

for read1 in *_unmapped_reads_1.fq.gz
do
    read2=$(echo $read1| sed 's/R1/R2/' | sed 's/_unmapped_reads_1/_unmapped_reads_2/')
    file_name=$(echo $read1| sed 's/.fastq.gz//')
    file_name_check=$file_name\_bismark_bt2_PE_report.txt
    if [ ! -f $temp_path/$unmapped_path/$hydroxymethyl_control/$file_name_check ]; then
	echo $file_name_check does not exist yet
	echo Run bismark on $read1 and $read2
	bismark --bam --maxins 800 $temp_genomes_path/$hydroxymethyl_control_fasta -1 $read1 -2 $read2 -o $unmapped_path/$hydroxymethyl_control --unmapped --nucleotide_coverage --multicore $cores
#--parallel ${SLURM_CPUS_PER_TASK}
    else
	echo Skipping Bismark mapping for $read1
    fi
done

#exit 1

#copy all results back to $seq_path and updates to the genome
#cp -R -u -p $temp_path $seq_path
rsync -vur $temp_path/ $seq_path
rsync -vur $temp_genomes_path/ /oak/stanford/groups/smontgom/kameronr/methylseq/JG97/genomes/ 


