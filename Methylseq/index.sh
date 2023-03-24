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

#copy over output files if they are already made
index_output_name=$(basename "${index_output}")
index_input_name=$(basename "${index_input}")

#if there is not the sort output in the output directory, then transfer files needed to make it and make it. else skip
if [ ! -f "$output_directory/$index_output_name" ]; then
    #transfer files needed if not present yet in temp directory.
    if [ ! -f "$index_input" ]; then
        rsync -vur --include="${index_input_name}" --exclude="*" "$output_directory/" $output_temp_directory
    fi

    echo "Starting to index $(basename "$index_input")"  
    samtools index $index_input
    echo "indexing of $(basename "$index_input") is complete"

    rsync -vur $output_temp_directory/ $output_directory

else
    echo "indexing of $(basename "$index_input") already complete"
fi


