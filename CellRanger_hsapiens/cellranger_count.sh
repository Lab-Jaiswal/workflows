#!/bin/bash

#SBATCH --job-name=Cell_Ranger_Count
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=batch
#SBATCH --account=sjaiswal
#SBATCH --time=16:00:00
#SBATCH --mem=128G
#SBATCH --mail-type=ALL

parent_directory=$1
output_directory=$2

line_number=$SLURM_ARRAY_TASK_ID
fastq_file="${parent_directory}/fastq_files" #provide path to file containing list of fastq files
fastq_path="$(sed "${line_number}q; d" $fastq_file)" #extract only the line number corresponding to $SLURM_ARRAY_TASK_ID

module load cellranger
transcriptome="/labs/sjaiswal/genomes/CellRanger_GRCh38/refdata-gex-GRCh38-2020-A"
transcriptome_temp="$(basename ${transcriptome})"

sample_prefix=$(ls $fastq_path | grep fastq | head -n 1 | sed -e 's/_S[0-9]*_L[0-9]*_[IR][0-9]_[0-9]*.fastq.gz//g')
filename=$(basename $fastq_path)
fastq_temp="${filename}_FASTQ"

cd /tmp
rm -rf ${fastq_temp} #delete any data from previous runs
rm -rf ${sample_prefix} #delete any data from previous runs

echo "copying FASTQs..."
rsync -PrhLtv $fastq_path ${fastq_temp} #copy data to /tmp
echo "...done"

echo "copying reference transcriptome..."
rsync -PrhLtv ${transcriptome} .
echo "...done"

cellranger count \
    --localmem=128 \
    --localcores=16 \
    --id="${filename}" \
    --transcriptome=${transcriptome_temp} \
    --fastqs=${fastq_temp} \
    --sample=${sample_prefix} 

rm -rf ${fastq_temp}
cd ${filename} 
tar --remove-files -zvcf SC_RNA_COUNTER_CS.tar.gz SC_RNA_COUNTER_CS
cd ..

rsync -PrhLtv ${filename} ${output_directory}
