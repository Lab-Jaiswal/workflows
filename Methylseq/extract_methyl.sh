#!/bin/bash
echo "entering map_and_deduplicate script"

bismark_input=$1
bismark_output=$2
output_temp_directory=$3
output_directory=$4
genome_fasta_path=$5
cores=$6
parameter_file=$7

echo "map_and_deduplicate command used the following parameters:
$0 $1 $2 $3 $4 $5 $6 $7"

if [ $SLURM_ARRAY_TASK_ID -eq 1 ]; then
        echo "bismark_input=$1
        bismark_output=$2
        output_temp_directory=$3
        output_directory=$4
        genome_fasta_path=$5
        cores=$6
        parameter_file=$7
        " >> $parameter_file
fi

module load bismark/0.22.3

if [ ! -f $bismark_output ]; then  
        echo "$bismark_output does not exist yet"
    bismark_methylation_extractor --gzip --cytosine_report --bedGraph --genome_folder "$genome_fasta_path" $bismark_input -o $output_temp_directory --multicore $cores
        echo "extract_methylation_controls complete for unmethyl control"
else
    echo "methylation control extraction for the unmethyl control found and already created"
fi

rsync -vur $output_temp_directory/ $output_directory
