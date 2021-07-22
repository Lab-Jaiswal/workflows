#!/bin/bash

min_coverage="10"
min_var_freq="0.001"
p_value="0.1"
amplicon_annot="/oak/stanford/groups/sjaiswal/Herra/CHIP_TWIST-PANEL_ATHEROMA/chip_submitted_targets_Twist.xls"

if [ -z $1 ] || [ -z $2 ]; then
    echo "run_BWA_varscan [fastq_directory] [output_directory]"
    echo "fastq_directory: path to raw .fastq or .fastq.gz files"
    echo "output_directory: path for BWA and varscan output"
else
    fastq_directory=$1 #get directory path from second argument (first argument $0 is the path of this script)
    output_directory=$2
    parent_directory=$(dirname $fastq_directory) #get parent directory of $fastq_directory
    code_directory="/labs/sjaiswal/workflows/BWA_varscan_hsapiens" #specify location of star_align_and_qc.sh
    fastq_list="${parent_directory}/fastq_files" #give a path to a file to store the paths to the fastq files in $fastq_directory

    mkdir -p $output_directory

    find "${fastq_directory}/" -type f | grep fastq | grep -v "Undetermined" | sed -e 's/_R1.*$//g' | sed -e 's/_R2.*$//g' | sort -u > ${fastq_list} #generate list of full paths to fastq files and save to the file in $fastq_list
    array_length=$(wc -l < ${fastq_list}) #get the number of files 

    sbatch -o "${output_directory}/variant_calling_log" `#put standard output into log` \
        -e "${output_directory}/variant_calling_log" `#put standard error into log` \
        -a "1-${array_length}" `#initiate job array equal to the number of fastq files` \
        "${code_directory}/BWA_varscan.sh" \
            ${parent_directory} \
            ${output_directory} \
            ${min_coverage} \
            ${min_var_freq} \
            ${p_value} \
            ${amplicon_annot}
fi
