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

#copy over output files if they are already made
sort_output_name=$(basename "${sort_output}")
sort_input_name=$(basename "${sort_input}")

#if there is not the sort output in the output directory, then transfer files needed to make it and make it. else skip
if [ ! -f "$output_directory/$sort_output_name" ]; then
    #transfer files needed if not present yet in temp directory.
    if [ ! -f "$sort_input" ]; then
        rsync -vur --include="${sort_input_name}" --exclude="*" "$output_directory/" $output_temp_directory
    fi

    echo "Expected sorted output file is '$sort_output'"
    echo "Starting to sort $(basename "$sort_input")"
    samtools sort $sort_input -o $sort_output
    echo "sorting of $sort_input is complete"

    rsync -vur $output_temp_directory/ $output_directory
else
    echo "sorting of $(basename "$sort_input") already complete"
fi


