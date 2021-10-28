#!/bin/bash

if [ -z $1 ] || [ -z $2 ]; then
    echo "submit_cellranger [fastq_directory] [max_jobs]"
    echo "fastq_directory: path to raw .fastq or .fastq.gz files"
    echo "output_directory: path for CellRanger output"
else
    TEMP=`getopt -o vdm: --long human,hsapiens,hsapien,mouse,mmusculus,musculus,human_nuclei,hsapiens_nuclei,hsapien_nuclei,nuclei \
    -n 'submit_cellranger.sh' -- "$@"`

    if [ $? != 0 ] ; then echo "Unrecognized argument. Possible arguments: human, human_nuclei, human_mkfastq, mouse, and mouse_cre. The argument represents the origin of the reference genome you would like aligned to your FASTQ files." >&2 ; exit 1 ; fi
        eval set -- "$TEMP"
        
        get_human=false
        get_human_nuclei=false
        get_mouse=false
                        
    while true; do
        case "$1" in
            --human | --hsapiens | --hsapien ) get_human=true; shift ;;
            --human_nuclei | hsapiens_nuclei | --hsapien_nulcei | nuclei ) get_human_nuclei=true; shift ;;
            --mouse | --mmusculus | musculus ) get_mouse=true; shift ;;
            -- ) shift; break ;;
            * ) break ;;
        esac
    done

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
        ${parent_directory} `#provide absolute path to fastq file` \
        ${output_directory} `#provide path for output` \
        ${get_human} \
        ${get_human_nuclei} \
        ${get_mouse} \
        
fi
