#!/bin/bash

#SBATCH --time=3:00:00
#SBATCH --account=sjaiswal
#SBATCH --cpus-per-task=8
#SBATCH --mem=32GB
#SBATCH --job-name=amplican

if [ -z "$1" ] || [ -z "$2" ]; then
    echo "./submit_amplican.sh [fastq_directory] [output_directory]"
    echo "fastq_directory: path to raw .fastq or .fastq.gz files"
    echo "output_directory: desired output path"
else
     TEMP=`getopt -o vdm: --long config: \
    -n 'submit_BWA_CHIP.sh' -- "$@"`

    if [ $? != 0 ] ; then echo "Unrecognized argument. Possible arguments: --config. Example: type --config = FILE_NAME, where FILE_NAME is the name of your config csv." >&2 ; exit 1 ; fi
        eval set -- "$TEMP"

        specify_config="NA"
        
    while true; do
        case "$1" in
            --config ) specify_config="$2"; shift 2 ;;
            -- ) shift; break ;;
            * ) break ;;
        esac
    done

    FASTQ_directory=$1
    output_directory=$2
    rm -r $output_directory
    mkdir $output_directory
    mkdir -p "$output_directory/Logs"

    module load R

    Rscript amplican.R  $FASTQ_directory $output_directory $specify_config \
        "$output_directory" > "$output_directory/Logs/amplicanOutFile.Rout" 2>&1
fi 
