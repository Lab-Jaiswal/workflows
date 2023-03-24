#!/bin/bash

#script to deduplicate the output from mapping.
echo ""
echo "entering deduplicate script"
echo ""

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

#copy over output files if they are already made
dedup_output_name=$(basename "${dedup_output}")
dedup_report_name=$(basename "${dedup_report}")

dedup_input_name=$(basename "${dedup_input}")



#if there is not the deduplicate output and deduplicate report in the output directory, then transfer files needed to make it and make it. else skip
	if [ ! -f "$output_directory/$dedup_output_name" ] || [ ! -f "$output_directory/$dedup_report_name" ]; then
		#transfer files needed if not present yet in temp directory.
		if [ ! -f "$dedup_input" ]; then
			rsync -vur --include="${dedup_input_name}" --exclude="*" "$output_directory/" $output_temp_directory
		fi

		echo "Begin deduplicating $dedup_input"
            	deduplicate_bismark -p --bam $dedup_input --output_dir $output_temp_directory -o $dedup_input 
            	#-o uses $dedup_input as output file basename because bismark modifies it to add deduplicated.bam ... -o just wants the basename to use for the output file.
            	#output_dir needed, otherwise outputs get written to the working directory!
        echo ""
        echo "Finished deduplicating $dedup_input"
        echo ""

        #checkpoint
		rsync -vur $output_temp_directory/ $output_directory

	else
	echo ""
	echo "$(basename "$dedup_input") has already been deduplicated"
	echo ""

    fi

#TO DO: This script could use testing code that tests each 'then' and 'else' sections of each 'if' statement! So a total of 4 test cases.

#TO DO: halt code with any error to prevent errors from snowballing and not being detected by end user or by development team. http://web.archive.org/web/20110314180918/http://www.davidpashley.com/articles/writing-robust-shell-scripts.html
