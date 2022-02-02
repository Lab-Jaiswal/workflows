#!/bin/bash

#script to map to genome using Bismark and optionally deduplicate the output from mapping.

while getopts t:T:g:o:c:d: flag
do
    case "${flag}" in
        t) read1_trimmed=${OPTARG};;
        T) read2_trimmed=${OPTARG};;
        g) genome_fasta_path=${OPTARG};;
        o) output_directory=${OPTARG};;
        c) cores=${OPTARG};;
        d) deduplicate=${OPTARG};;

    esac
done

echo ""
echo "mapping command used:"
echo "$0 -t $read1_trimmed -T $read2_trimmed -g $genome_fasta_path -o $output_directory -c $cores -d $deduplicate"
echo ""

module load bismark/0.22.3


#check if output directory exists yet, if it does do work there, else throw error message to fix upstream code.

if [ -d "$output_directory" ] 
	then
		echo "Output directory for mapping is "$output_directory""
		#if statement to stopcode from running this section if already done
		"$(dirname "$read1_trimmed")"
		expected_mapping_output="$(basename "$read1_trimmed" ".fastq.gz")_bismark_bt2_PE_report.txt"
		#expected_mapping_output=${output_directory}/${read1_trimmed##*/}_bismark_bt2_PE_report.txt
		echo "Expected mapping output file is $expected_mapping_output"
		if [ ! -f "$expected_mapping_output" ] 
			then
			#I don't think this part of the code needs to cd to the output directory, so it doesn't. 
			#map code
			echo "$(basename "$read1_trimmed") and $(basename "$read2_trimmed") read pair not yet mapped to "$genome" genome. Mapping now."
			bismark --bam --maxins 800 $genome_fasta_path -1 $read1_trimmed -2 $read2_trimmed -o $output_directory --unmapped --nucleotide_coverage --multicore $cores
			echo "finished mapping $(basename $read1_trimmed) and $(basename $read2_trimmed) read pair."

			#optional rsync to copy files back to seq_path directory
			#rsync -vur --exclude "main_genome" --exclude "unmethyl_genome" --exclude "hydroxymethyl_genome" --exclude "methyl_genome" $temp_path/ $seq_path

			else
		        echo "$(basename "$read1") and $(basename "$read2") have already been Bismark mapped to $genome_fasta_path"
		fi

	else
		echo "Error: Directory "$output_directory" does not exist. Directory "$output_directory" must be created before running "$0""
fi

#TO DO: This script could use testing code that tests each 'then' and 'else' sections of each 'if' statement! So a total of 4 test cases.

#TO DO: halt code with any error to prevent errors from snowballing and not being detected by end user or by development team. http://web.archive.org/web/20110314180918/http://www.davidpashley.com/articles/writing-robust-shell-scripts.html

#TO DO: 
#should do test on relative and absolute paths and ensure that it is mapping the files. it's looking for the relative path . which exists always. It should map and output the mapped files to the current working directory.


#deduplicate step

#check if deduplication is wanted
if [ "$deduplicate" == "True" ] || [ "$deduplicate" == "true" ] || [ "$deduplicate" == "TRUE" ]
	then
		echo ""
		echo "Deduplication requested for mapping output of $(basename "$read1") and $(basename "$read2") to genome $genome_fasta_path"
		#expected_deduplication_output="${output_directory}/${read1_trimmed##*/}_bismark_bt2_pe.deduplicated.bam"
		expected_deduplication_output="${output_directory}/$(basename "$read1_trimmed" ".fastq.gz")_bismark_bt2_pe.deduplication_report.txt.gz"
		#TO DO: this path needs to be checked that it is correct!
		echo "Expected deduplication output file is $expected_deduplication_output"
		if [ ! -f "$expected_deduplication_output" ] 
			then
				#deduplicate code
				bam="${output_directory}/$(basename "$read1_trimmed" ".fastq.gz")_bismark_bt2_pe.bam"
				#TO DO: check bam file is correct!
				echo "Begin deduplicating $bam"
            	deduplicate_bismark -p --bam $bam  
            	echo "Finished deduplicating $bam"
        	else
        		echo "$(basename "$bam") has already been deduplicated by Bismark."
        fi
    else
    	echo "Deduplication NOT requested. ( If this is not desired then set '-d true' in the command $0 )"
fi








