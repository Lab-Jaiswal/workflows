#!/bin/bash

min_coverage="10"
min_var_freq="0.001"
p_value="0.1"

if [ -z "$1" ] || [ -z "$2" ]; then
    echo "run_BWA_mutect [fastq_directory] [output_directory]"
    echo "fastq_directory: path to raw .fastq or .fastq.gz files"
    echo "output_directory: path for BWA and mutect output"
else
    fastq_directory=$1 #get directory path from second argument (first argument $0 is the path of this script)
    output_directory=$2
    parent_directory=$(dirname $fastq_directory) #get parent directory of $fastq_directory
    code_directory="/labs/sjaiswal/workflows/BWA_CHIP" #specify location of star_align_and_qc.sh
    fastq_list="${parent_directory}/fastq_files" #give a path to a file to store the paths to the fastq files in $fastq_directory

    mkdir -p $output_directory
    mkdir -p "${output_directory}/Logs"

    find "${fastq_directory}/" -type f `#list all files in ${fastq_directory}` | \
        grep "fastq" `#only keep files with FASTQ in name (case insensitive)` | \
        grep -v "Undetermined" `#remove Undetermined FASTQs` | \
        sed -e 's/_R1.*$//g' | sed -e 's/_R2.*$//g' `#remove _R1/2_fastq.gz file extension`| \
        sort -u > ${fastq_list} #sort and remove duplicate names, to generate list of unique FASTQs
    array_length=$(wc -l < ${fastq_list}) #get the number of FASTQs 

    sbatch -o "${output_directory}/Logs/%A_%a.log" `#put into log` \
        -a "1-${array_length}" `#initiate job array equal to the number of fastq files` \
        "${code_directory}/BWA_CHIP.sh" \
            ${parent_directory} \
            ${output_directory} \
            ${min_coverage} \
            ${min_var_freq} \
            ${p_value}
fi
