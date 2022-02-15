#!/bin/bash

#script to map to genome using Bismark and optionally deduplicate the output from mapping.

echo "entering map_and_deduplicate script"

read1_input=$1
read2_input=$2
read1_bismark=$3
dedup_input=$4
dedup_output=$5
genome_fasta_path=$6
output_temp_directory=$7
output_directory=$8
cores=$9
deduplicate=${10}
parameter_file=${11}

echo "trim command used the following parameters:
$0 $1 $2 $3 $4 $5 $6 $7 $8 $9 ${10} ${11}"

if [ $SLURM_ARRAY_TASK_ID -eq 1 ]; then
	echo "command given to map_and_deduplicate.sh includes the following parameters:
		read1_input=$1
		read2_input=$2
		read1_bismark=$3
		dedup_input=$4
		dedup_output=$5
		genome_fasta_path=$6
		output_temp_directory=$7
		output_directory=$8
		cores=$9
		deduplicate=${10}
		parameter_file=${11}
		" >> $parameter_file

module load bismark/0.22.3
        

if [ ! -f $read1_bismark ]; then
	echo "Expected mapping output file is $read1_bismark"
	echo "$(basename "$read1_input") and $(basename "$read2_input") read pair not yet mapped to "$genome_name" genome. Mapping now."
	echo "bismark command:
		bismark --bam --maxins 800 $genome_fasta_path -1 $read1_input -2 $read2_input -o $output_temp_directory --unmapped --nucleotide_coverage --multicore $cores"

	bismark --bam --maxins 800 $genome_fasta_path -1 $read1_input -2 $read2_input -o $output_temp_directory --unmapped --nucleotide_coverage --multicore $cores
			
	echo "finished mapping $(basename $read1_input) and $(basename $read2_input) read pair; map to control seqs complete for $genome_name"
else
	echo "$(basename "$read1_input") and $(basename "$read2_input") have already been Bismark mapped to $genome_fasta_path"
fi

	
if [ "$deduplicate" == "TRUE" ] || [ "$deduplicate" == TRUE ] || [ "$deduplicate" == "true" ] || [ "$deduplicate" == "True" ]; then
	echo "Deduplication requested for mapping output of $(basename "$read1_input") and $(basename "$read2_input") to genome $genome_fasta_path"
	echo "Expected deduplication output file is $dedup_output"

	if [ ! -f "$dedup_output" ]; then
		cd $output_temp_directory
                deduplication_input=$(basename $dedup_input)
		echo "Begin deduplicating $dedup_input"
            	deduplicate_bismark -p --bam $dedup_input -o $dedup_input
            	echo "Finished deduplicating $dedup_input"
        fi
    
    else        
    	echo "Deduplication NOT requested. ( If this is not desired then set '-d true' in the command $0 )"
fi

rsync -vur $output_temp_directory/ $output_directory

#TO DO: This script could use testing code that tests each 'then' and 'else' sections of each 'if' statement! So a total of 4 test cases.

#TO DO: halt code with any error to prevent errors from snowballing and not being detected by end user or by development team. http://web.archive.org/web/20110314180918/http://www.davidpashley.com/articles/writing-robust-shell-scripts.html

#TO DO: 
#should do test on relative and absolute paths and ensure that it is mapping the files. it's looking for the relative path . which exists always. It should map and output the mapped files to the current working directory. 
