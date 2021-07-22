#!/bin/bash

if [ -z $1 ]; then
    echo "submit_star_align_and_qc [fastq_directory] [max_jobs]"
    echo "fastq_directory: path to raw .fastq or .fastq.gz files"
else
    fastq_directory=$1 #get directory path from second argument (first argument $0 is the path of this script)
    output_directory=$(dirname $fastq_directory) #get parent directory of $fastq_directory
    code_directory="/labs/sjaiswal/workflows/RNASeq_mmusculus" #specify location of star_align_and_qc.sh
    fastq_list="${output_directory}/fastq_files" #give a path to a file to store the paths to the fastq files in $fastq_directory

    mkdir "${output_directory}/Logs"

    ### Do NOT make this value too large. Running too many simultaneous jobs WILL lock up the server.  
    #if [ -z $2 ]; then
        #max_jobs=3 #specify the number of intervals to split your jobs into.
    #else
        #max_jobs=$2    
    #fi

    find "${fastq_directory}/" -type f | grep fastq | grep -v "Undetermined" | sed -e 's/_R1.*$//g' | sed -e 's/_R2.*$//g' | sort -u > ${fastq_list} #generate list of full paths to fastq files and save to the file in $fastq_list
    array_length=$(wc -l < ${fastq_list}) #get the number of files 

    #Note - do NOT delete the the '\' or '`' characters - they allow the command to have multiple lines with comments!
    sbatch -o "${output_directory}/Logs/%A_%a.log" `#put into log` \
        -a "1-${array_length}" `#initiate job array equal to the number of fastq files` \
        "${code_directory}/star_align_and_qc.sh" `#specify STAR alignment script` \
        ${output_directory} `#provide absolute path to fastq file`
fi


