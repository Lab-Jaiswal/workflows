#!/bin/bash

#SBATCH --job-name=methylseq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --partition=batch
#SBATCH --mem=500G
#SBATCH --time=150:00:00
#SBATCH --account=sjaiswal


####SBATCH --partititon=interactive
####SBATCH --partition=nih_s10
####SBATCH --account=sjaiswal

#data_path=/oak/stanford/groups/sjaiswal/kameronr/JG97/methylseq/miseq/214460250
data_path=$1
unmethyl_control=$2
unmethyl_control_fasta=$3
hydroxymethyl_control=$4
hydroxymethyl_control_fasta=$5
genome_path=$6
phix_path=$7
temp_path=$8

if ! [ -d "$temp_path" ]; then
    #${1:-$temp_path}
    temp_path=$(mktemp -d /tmp/tmp.XXXXXXXXXX)
    echo "temp_path is: " $temp_path
    mkdir $temp_path
    #copy fastq files to temp_path
    rsync -vur $data_path/fastq/ $temp_path

    temp_genomes_path=/tmp/genomes
    rsync -vur $genome_path $temp_genomes_path
fi

echo this is temp_path = $temp_path

#cores=$(expr ${SLURM_CPUS_PER_TASK}/2)
cores=24
echo core = $cores

bash control_ref_genome_prep.sh

bash phix_ref_genome_prep.sh

bash process_miseq.sh $data_path
bash trim.sh $data_path $temp_path
bash map_to_control_seqs.sh $data_path $unmethyl_control $unmethyl_control_fasta $hydroxymethyl_control $hydroxymethyl_control_fasta $temp_path $cores
bash extract_methylation_controls.sh $data_path $unmethyl_control $unmethyl_control_fasta $hydroxymethyl_control $hydroxymethyl_control_fasta $temp_path $cores
bash report_controls.sh $data_path $unmethyl_control $hydroxymethyl_control $temp_path

bash map_to_genome_seqs.sh $data_path $unmethyl_control $hydroxymethyl_control $genome_path $phix_path $temp_path $cores
bash remove_duplicates.sh $data_path $unmethyl_control $hydroxymethyl_control $temp_path
bash insert_size_analysis.sh $data_path $unmethyl_control $hydroxymethyl_control $temp_path
bash index_bam.sh $data_path $unmethyl_control $hydroxymethyl_control $temp_path
bash extract_methylation.sh $data_path $unmethyl_control $hydroxymethyl_control $genome_path $temp_path $cores
bash report.sh $data_path $unmethyl_control $hydroxymethyl_control $temp_path

