#!/bin/bash

#SBATCH --job-name=ref_genome_prep
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=interactive
#SBATCH --mem=64G
#SBATCH --time=16:00:00

module load bismark/0.22.3

#cd /oak/stanford/groups/sjaiswal/kameronr/JG97/methylseq/miseq/

cd /oak/stanford/groups/sjaiswal/kameronr/sjaiswal_old/genomes/GRCh38/

bismark_genome_preparation ./pUC19_fasta
