#!/bin/bash

#script to extract methylation with error checking code wrapped around it. 

#TO DO: test all the paths in this script!

bismark_input=$1
bismark_output=$2
output_directory=$3
genome_fasta_path=$4
cores=$5

echo "bismark_input=$bismark_input
bismark_output=$bismark_output
output_directory=$output_directory
genome_fasta_path=$genome_fasta_path
cores=$cores
"

module load bismark/0.22.3


if [ ! -f $bismark_output ]; then  
        echo "$bismark_output does not exist yet"
    bismark_methylation_extractor --gzip --cytosine_report --bedGraph --genome_folder "$genome_fasta_path" $bismark_input --multicore $cores
        echo "extract_methylation_controls complete for unmethyl control"
else
    echo "methylation control extraction for the unmethyl control found and already created"
fi
