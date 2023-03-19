#!/bin/bash

#script to deduplicate the output from mapping.

echo "entering deduplicate script"

dedup_input=$1
dedup_output=$2
dedup_report=$3
output_temp_directory=$4
output_directory=$5
cores=$6
parameter_file=$7

echo "deduplicate command used the following parameters:
		dedup_input=$1
		dedup_output=$2
		dedup_report=$3
		output_temp_directory=$4
		output_directory=$5
		cores=$6
		parameter_file=$7"

if [ $SLURM_ARRAY_TASK_ID -eq 1 ]; then
	echo "command given to deduplicate.sh includes the following parameters:
		dedup_input=$1
		dedup_output=$2
		dedup_report=$3
		output_temp_directory=$4
		output_directory=$5
		cores=$6
		parameter_file=$7
		" >> $parameter_file
fi

module load bismark/0.22.3        

	echo "Deduplication requested for mapping output $dedup_input"
	echo "Expected deduplication output file is $dedup_output"

	if [ ! -f "$dedup_output" ] || [ ! -f "$dedup_report" ]; then
		echo "Begin deduplicating $dedup_input"
            	deduplicate_bismark -p --bam $dedup_input --output_dir $output_temp_directory -o $dedup_input #this uses $dedup_input as output file basename because bismark modifies it to add deduplicated.bam 
            	#output_dir needed, otherwise outputs get written to the working directory!
            	echo "Finished deduplicating $dedup_input"
    fi

#checkpoint
rsync -vur $output_temp_directory/ $output_directory

#TO DO: This script could use testing code that tests each 'then' and 'else' sections of each 'if' statement! So a total of 4 test cases.

#TO DO: halt code with any error to prevent errors from snowballing and not being detected by end user or by development team. http://web.archive.org/web/20110314180918/http://www.davidpashley.com/articles/writing-robust-shell-scripts.html
