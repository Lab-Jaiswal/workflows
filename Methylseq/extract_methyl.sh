#!/bin/bash
echo "entering extract_methyl script"

bismark_extraction_input=$1
bismark_extraction_report=$2
output_temp_directory=$3
output_directory=$4
genome_fasta_path=$5
cores=$6
parameter_file=$7

echo "extract_methyl command used the following parameters:
$0 $1 $2 $3 $4 $5 $6 $7"

if [ $SLURM_ARRAY_TASK_ID -eq 1 ]; then
        echo "bismark_extraction_input=$1
        bismark_extraction_report=$2
        output_temp_directory=$3
        output_directory=$4
        genome_fasta_path=$5
        cores=$6
        parameter_file=$7
        " >> $parameter_file
fi

module load bismark/0.22.3

N_cores=$((cores/3))

if [ ! -f $bismark_extraction_report ]; then  #TO DO: should add requirement to also have the extracted methylation files too! Not just the report.
        echo "$bismark_extraction_report does not exist yet"
    bismark_methylation_extractor --gzip --cytosine_report --bedGraph --genome_folder "$genome_fasta_path" $bismark_extraction_input -o $output_temp_directory --multicore $N_cores
        echo "extract_methylation_controls complete for unmethyl control"
else
    echo "methylation control extraction for the unmethyl control found and already created"
fi

rsync -vur $output_temp_directory/ $output_directory
