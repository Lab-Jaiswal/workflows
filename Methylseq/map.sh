#!/bin/bash

#script to map to genome using Bismark

echo "entering map script"

read1_input=$1
read2_input=$2
bismark_map_output=$3
bismark_map_report=$4
genome_fasta_path=$5
output_temp_directory=$6
output_directory=$7
cores=$8
parameter_file=$9
previous_loop_output_directory=${10}

echo "map command used the following parameters:
		read1_input=$1
		read2_input=$2
		bismark_map_output=$3
		bismark_map_report=$4
		genome_fasta_path=$5
		output_temp_directory=$6
		output_directory=$7
		cores=$8
		parameter_file=$9
		previous_loop_output_directory=$10"

if [ $SLURM_ARRAY_TASK_ID -eq 1 ]; then
	echo "command given to map.sh includes the following parameters:
		read1_input=$1
		read2_input=$2
		bismark_map_output=$3
		bismark_map_report=$4
		genome_fasta_path=$5
		output_temp_directory=$6
		output_directory=$7
		cores=$8
		parameter_file=$9
		previous_loop_output_directory=$10
		" >> $parameter_file
fi

module load bismark/0.22.3        

#copy over output files if they are already made
bismark_map_output_name=$(basename "${bismark_map_output}")
bismark_map_report_name=$(basename "${bismark_map_report}")

read1_input_name=$(basename "${read1_input}")
read2_input_name=$(basename "${read2_input}")

previous_loop_temp_directory=$(dirname "${read1_input}")


# rsync -vur --include="${read1_name}" --include="${read2_name}" --include="${read1_trimmed_name}" --include="${read2_trimmed_name}" --exclude="*" "${data_path}/fastq/" $temp_path

#if there is the map output name and the map report, then skip without transferring (the next step will transfer the files just prior if the next step needs to be done) ... basically don't transfer file until the last second just before it will be needed
#else make the files
#rsync -vur --include="${bismark_map_output_name}" --include="${bismark_map_report_name}" --exclude="*" "$output_directory/" $output_temp_directory

#if there is not the map output and map report in the output directory, then transfer files needed to make it and make it. else skip

#if [ ! -f $bismark_map_output ] || [ ! -f $bismark_map_report ]; then
if [ ! -f $output_directory/$bismark_map_output_name ] || [ ! -f $output_directory/$bismark_map_report_name ]; then
	#transfer files needed if not present yet in temp directory. for this if it's the 1st loop then not needed, but for all other loops it is needed.
	if [ ! -f $read1_input ] || [ ! -f $read2_input ]; then
		rsync -vur --include="${read1_input_name}" --include="${read2_input_name}" --exclude="*" "$previous_loop_output_directory/" $previous_loop_temp_directory
	fi
	#rsync -vur --include="${bismark_map_output_name}" --include="${bismark_map_report_name}" --exclude="*" "$output_directory/" $output_temp_directory
	echo "Expected mapping output file is $bismark_map_report"
	echo "$(basename "$read1_input") and $(basename "$read2_input") read pair not yet mapped to '$genome_name' genome. Mapping now."
	echo "bismark command:
		bismark --bam --maxins 800 $genome_fasta_path -1 $read1_input -2 $read2_input -o $output_temp_directory --unmapped --nucleotide_coverage --multicore $cores"

#TO DO: if mapping is already done, then check if the next step after is done and if not then transfer the unmapped files and things needed for the next step to be done

#more concervative mapping
#	bismark --bam --maxins 800 $genome_fasta_path -1 $read1_input -2 $read2_input -o $output_temp_directory --unmapped --nucleotide_coverage --multicore $cores

N_cores=$((cores/3))

#more liberal mapping
	bismark --bam --maxins 800 $genome_fasta_path -1 $read1_input -2 $read2_input -o $output_temp_directory --unmapped --nucleotide_coverage --multicore $N_cores --score_min L,0,-0.6

	echo "finished mapping $(basename $read1_input) and $(basename $read2_input) read pair; map to control seqs complete for $genome_name"

	rsync -vur $output_temp_directory/ $output_directory

else
	echo "$(basename "$read1_input") and $(basename "$read2_input") have already been Bismark mapped to $genome_fasta_path"
fi

#checkpoint



#TO DO: This script could use testing code that tests each 'then' and 'else' sections of each 'if' statement! So a total of 4 test cases.

#TO DO: halt code with any error to prevent errors from snowballing and not being detected by end user or by development team. http://web.archive.org/web/20110314180918/http://www.davidpashley.com/articles/writing-robust-shell-scripts.html
