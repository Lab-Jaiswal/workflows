#!/bin/bash

#SBATCH --time=12:00:00
#SBATCH --account=sjaiswal
#SBATCH --cpus-per-task=16
#SBATCH --mem=128GB
#SBATCH --job-name=methylseq
#####SBATCH --partition=nih_s10

data_path=$1
seq_path="$1/fastq"
unmethyl_control_fasta=$2
unmethyl_control=$3
hydroxymethyl_control_fasta=$4
hydroxymethyl_control=$5
methyl_control_fasta=$6
methyl_control=$7
genome_path=$8
phix_path=$9
cores=${10}

line_number=$SLURM_ARRAY_TASK_ID #get index of which file to process from $SLURM_ARRAY_TASK_ID provided by SLURM
bam_file="$data_path/fastq/$unmethyl_control/$hydroxymethyl_control/$methyl_control/genome_alignment/bam_files" #provide path to file containing list of fastq files
bam_path="$(sed "${line_number}q; d" $bam_file)" #extract only the line number corresponding to $SLURM_ARRAY_TASK_ID

#sample_prefix=$(ls $bam_path | grep bam | head -n 1 | sed -e 's/_S[0-9]*_L[0-9]*_[IR][0-9]_[0-9]*.*$//g')

#echo $sample_prefix

filename=$(basename "${bam_path}")
bam_temp="${filename}"
PREFIX=$filename

seq_path="$data_path/fastq"

echo "copying bams..."
temp_path=$(mktemp -d /tmp/tmp.XXXXXXXXXX)
echo "temp_path is: " $temp_path
mkdir $temp_path
#copy fastq files to temp_path
rsync -vur "$data_path/fastq/" $temp_path

echo "bams have been copied to the temporary file"

echo "copying reference transcriptome..."

mkdir /tmp/genomes
mkdir -p /tmp/genomes/main_genome
mkdir -p /tmp/genomes/unmethyl_genome
mkdir -p /tmp/genomes/hydroxymethyl_genome
mkdir -p /tmp/genomes/methyl_genomes

rsync -vur "$genome_path/" "$temp_path/main_genome"
#genomes_path=$(dirname $genome_path)
#genome_name=$(basename $genome_path)
#echo $genomes_path
rsync -vur "$unmethyl_control_fasta/" "$temp_path/unmethyl_genome"
rsync -vur "$hydroxymethyl_control_fasta/" "$temp_path/hydroxymethyl_genome" 
rsync -vur "$methyl_control_fasta/" "$temp_path/methyl_genome" 
echo "Genome have been copied to the temporary file directory"

module load bismark/0.22.3
module load samtools
"Bismark and samtools modules have been loaded"

bam_sample_path=$(echo $bam_path| sed 's/.bam$//')




#####################previous extract_methylation.sh############################
################################################################################
if [ ! -f "$temp_path/$unmethyl_control/$hydroxymethyl_control/$methyl_control/genome_alignment/split_bams/${bam_sample_path}.bismark.cov.gz" ]; then 
    cd "$temp_path/$unmethyl_control/$hydroxymethyl_control/$methyl_control/genome_alignment/split_bams"
        echo "$temp_path/$unmethyl_control/$hydroxymethyl_control/$methyl_control/genome_alignment/split_bams/${bam_sample_path}.bismark.cov.gz"
    #remove --cytosine_report
    bismark_methylation_extractor --buffer_size 20G --gzip --bedGraph --genome_folder "$temp_path/main_genome" $bam_path --multicore $cores
        echo "extract_methylation complete"
    rsync -vur --exclude "main_genome" --exclude "unmethyl_genome" --exclude "methyl_genome" --exclude "hydroxymethyl_genome" $temp_path/ $seq_path
else 
    echo "methlyation extraction already completed"
fi
