#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH --account=smontgom
#SBATCH --cpus-per-task=8
#SBATCH --mem=32GB
#SBATCH --job-name=bcl2fastq

input_directory=$1
output_directory=$2
cores=8

#script_directory=$(dirname "$0")

bash "bcl_to_fastq.sh" -i "${input_directory}" -o "${output_directory}" -p $cores


