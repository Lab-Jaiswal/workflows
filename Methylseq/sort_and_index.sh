#!/bin/bash
echo "entering sort_and_index script"

sort_input=$1
index_input=$2
index_output=$3
output_temp_directory=$4
output_directory=$5
parameter_file=$6

echo "sort_and_index command used the following parameters:
$0 $1 $2 $3 $4 $5 $6"

if [ $SLURM_ARRAY_TASK_ID -eq 1 ]; then
        echo "arguments used for the sort_and_index.sh script:
                sort_input: $1
                index_input: $2
                index_output: $3
                output_temp_directory: $4
                output_directory: $5
                parameter_file: $6
                " >> $parameter_file
fi

module load samtools/1.9

if [ ! -f "$index_input" ]; then
    echo "Expected sorted output file is '$index_input'"
    echo "Starting to sort $(basename "$sort_input")"
    samtools sort $sort_input -o $index_input
    echo "sorting of $sort_input is complete"
else
    echo "sorting of $(basename "$sort_input") already complete"
fi

if  [ -f "$index_output" ]; then
    echo "Starting to index $(basename "$index_input")"  
    samtools index $index_input
    echo "indexing of $(basename "$index_input") is complete"
else
    echo "sorting of $(basename "$index_input") already complete"
fi

rsync -vur $output_temp_directory/ $output_directory

    	



