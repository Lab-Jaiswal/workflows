#!/bin/bash

#SBATCH --job-name=ref_genome_prep
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --partition=interactive
#SBATCH --mem=4G
#SBATCH --time=1:00:00

module load bismark/0.22.3

cd /oak/stanford/groups/sjaiswal/kameronr/JG97/methylseq/miseq/

bismark_genome_preparation ./fasta_active_motif_sequence

bismark_genome_preparation ./fasta_zymo_sequence







