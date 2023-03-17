#!/bin/bash
echo "entering sort script"

sort_input=$1
sort_output=$2
output_temp_directory=$3
output_directory=$4
parameter_file=$5

echo "sort command used the following parameters:
$0 $1 $2 $3 $4 $5"

if [ $SLURM_ARRAY_TASK_ID -eq 1 ]; then
        echo "arguments used for the sort.sh script:
                sort_input: $1
                sort_output: $2
                output_temp_directory: $3
                output_directory: $4
                parameter_file: $5
                " >> $parameter_file
fi

module load samtools/1.9

if [ ! -f "$sort_output" ]; then
    echo "Expected sorted output file is '$sort_output'"
    echo "Starting to sort $(basename "$sort_input")"
    samtools sort $sort_input -o $sort_output
    echo "sorting of $sort_input is complete"
else
    echo "sorting of $(basename "$sort_input") already complete"
fi

rsync -vur $output_temp_directory/ $output_directory
