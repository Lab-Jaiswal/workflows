#!/bin/bash

#script to extract methylation with error checking code wrapped around it. 

#TO DO: test all the paths in this script!

while getopts b:o: flag
do
    case "${flag}" in
        b) bam=${OPTARG};;
        o) output_directory=${OPTARG};;
    esac
done



module load bismark/0.22.3


if [ ! -f "$temp_path/$unmethyl_control/${trimmed_R1_file_name}_bismark_bt2_pe_splitting_report.txt" ]; then  
    cd $temp_path/$unmethyl_control
        echo "$temp_path/$unmethyl_control/${trimmed_R1_file_name}_bismark_bt2_pe_splitting_report.txt does not exist yet"
    bismark_methylation_extractor --gzip --cytosine_report --bedGraph --genome_folder "$temp_path/unmethyl_genome" $trimmed_R1_bismark_bam --multicore $cores
        echo "extract_methylation_controls complete for unmethyl control"
    rsync -vur --exclude "main_genome" --exclude "unmethyl_genome" --exclude "hydroxymethyl_genome" --exclude "methyl_genome" $temp_path/ $seq_path
else
    echo "methylation control extraction for the unmethyl control found and already created"
fi







