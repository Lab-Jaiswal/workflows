!/bin/bash

#script to convert bcl files to fastq file with error messaging also.

#TO DO: test all the paths in this script!

while getopts i:o:p: flag
do
    case "${flag}" in
        i) input_directory=${OPTARG};;
        o) output_directory=${OPTARG};;
	    p) processing_threads=${OPTARG};;
    esac
done

echo ""
echo "bcl_to_fastq.sh command used:"
echo "$0 -i $input_directory -o $output_directory -p $processing_threads"
echo ""

#TO DO: add code to check all inputs are present and if not to give error message with help information about all the script inputs

#TO DO: check for sample sheet and if no sample sheet then quit

module load bcl2fastq2/2.20.0.422

echo "Expected output directory is '$output_directory'"

if [ ! -d "${output_directory}" ]
    then
    echo "Starting to convert bcl files to fastq files for $input_directory"
    starting_directory=$(pwd)
    cd $input_directory
    bcl2fastq -o $output_directory -p $processing_threads
    cd $starting_directory
    echo "Finished converting bcl files to fastq files for $input_directory"
    else
    echo "Fastq directory (and likely files) already exists, so skipped converting bcl files to fastq files for $input_directory"
fi
