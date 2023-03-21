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

echo "map command used the following parameters:
		read1_input=$1
		read2_input=$2
		bismark_map_output=$3
		bismark_map_report=$4
		genome_fasta_path=$5
		output_temp_directory=$6
		output_directory=$7
		cores=$8
		parameter_file=$9"

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
		" >> $parameter_file
fi

module load bismark/0.22.3        

#copy over output files if they are already made
bismark_map_output_name=$(basename "${bismark_map_output}")
bismark_map_report_name=$(basename "${bismark_map_report}")
# rsync -vur --include="${read1_name}" --include="${read2_name}" --include="${read1_trimmed_name}" --include="${read2_trimmed_name}" --exclude="*" "${data_path}/fastq/" $temp_path
rsync -vur --include="${bismark_map_output_name}" --include="${bismark_map_report_name}" --exclude="*" "$output_directory/" $output_temp_directory

if [ ! -f $bismark_map_output ] || [ ! -f $bismark_map_report]; then
	echo "Expected mapping output file is $bismark_map_report"
	echo "$(basename "$read1_input") and $(basename "$read2_input") read pair not yet mapped to '$genome_name' genome. Mapping now."
	echo "bismark command:
		bismark --bam --maxins 800 $genome_fasta_path -1 $read1_input -2 $read2_input -o $output_temp_directory --unmapped --nucleotide_coverage --multicore $cores"

#more concervative mapping
#	bismark --bam --maxins 800 $genome_fasta_path -1 $read1_input -2 $read2_input -o $output_temp_directory --unmapped --nucleotide_coverage --multicore $cores

N_cores=$((cores/3))

#more liberal mapping
	bismark --bam --maxins 800 $genome_fasta_path -1 $read1_input -2 $read2_input -o $output_temp_directory --unmapped --nucleotide_coverage --multicore $N_cores --score_min L,0,-0.6

	echo "finished mapping $(basename $read1_input) and $(basename $read2_input) read pair; map to control seqs complete for $genome_name"
else
	echo "$(basename "$read1_input") and $(basename "$read2_input") have already been Bismark mapped to $genome_fasta_path"
fi

#checkpoint
rsync -vur $output_temp_directory/ $output_directory



#TO DO: This script could use testing code that tests each 'then' and 'else' sections of each 'if' statement! So a total of 4 test cases.

#TO DO: halt code with any error to prevent errors from snowballing and not being detected by end user or by development team. http://web.archive.org/web/20110314180918/http://www.davidpashley.com/articles/writing-robust-shell-scripts.html
