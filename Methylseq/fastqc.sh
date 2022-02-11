#!/bin/bash

#SBATCH --job-name=fastqc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=interactive
#SBATCH --mem=64G
#SBATCH --time=16:00:00

experiment_directory=$1

module load fastqc/0.11.9 

echo "Starting fastqc for $experiment_directory"

for file in *.fastq.gz ; do
	fastqc $file ;
	done

echo "Fastqc finished for $experiment_directory" 


