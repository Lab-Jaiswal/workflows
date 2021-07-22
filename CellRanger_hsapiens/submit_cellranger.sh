#!/bin/bash

if [ -z $1 ] || [ -z $2 ]; then
    echo "submit_cellranger [fastq_directory] [max_jobs]"
    echo "fastq_directory: path to raw .fastq or .fastq.gz files"
    echo "output_directory: path for CellRanger output"
else
    fastq_directory=$1 #get directory path from second argument (first argument $0 is the path of this script)
    parent_directory=$(dirname ${fastq_directory}) #get parent directory of $fastq_directory
    fastq_list="${parent_directory}/fastq_files" #give a path to a file to store the paths to the fastq files in $fastq_directory

    output_directory=$2 
    mkdir -p ${output_directory} #Creates output directory if it does not exist
    mkdir -p ${output_directory}/Logs #Creates output directory if it does not exist

    code_directory="/labs/sjaiswal/workflows/CellRanger_hsapiens" #specify location of star_align_and_qc.sh

    dirname $(find -L ${fastq_directory} -type f | grep fastq.gz | grep -v Undetermined) | sort -u > ${fastq_list}
    array_length=$(wc -l < ${fastq_list}) #get the number of files 

    #Note - do NOT delete the the '\' or '`' characters - they allow the command to have multiple lines with comments!
    sbatch -o "${output_directory}/Logs/%A_%a.log" `#put output into logs` \
        -a "1-${array_length}" `#initiate job array equal to the number of fastq files` \
        "${code_directory}/cellranger_count.sh" `#specify STAR alignment script` \
        ${parent_directory} `#provide absolute path to fastq file`\
        ${output_directory} `#provide path for output`
fi
