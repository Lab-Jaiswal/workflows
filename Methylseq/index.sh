#!/bin/bash
echo "entering index script"

index_input=$1
index_output=$2
output_temp_directory=$3
output_directory=$4
parameter_file=$5

echo "index command used the following parameters:
$0 $1 $2 $3 $4 $5"

if [ $SLURM_ARRAY_TASK_ID -eq 1 ]; then
        echo "arguments used for the sort_and_index.sh script:
                index_input: $1
                index_output: $2
                output_temp_directory: $3
                output_directory: $4
                parameter_file: $5
                " >> $parameter_file
fi

module load samtools/1.9

if  [ ! -f "$index_output" ]; then
    echo "Starting to index $(basename "$index_input")"  
    samtools index $index_input
    echo "indexing of $(basename "$index_input") is complete"
else
    echo "indexing of $(basename "$index_input") already complete"
fi

rsync -vur $output_temp_directory/ $output_directory
