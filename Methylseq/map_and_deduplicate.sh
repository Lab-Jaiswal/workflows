#!/bin/bash

#script to map to genome using Bismark and optionally deduplicate the output from mapping.

echo "entered map_and_deduplicate script"

read1_input=$1
read2_input=$2
read1_bismark=$3
dedup_input=$4
dedup_output=$5
genome_fasta_path=$6
output_directory=$7
cores=$8
deduplicate=$9

echo "read1_input=$1
read2_input=$2
read1_bismark=$3
dedup_input=$4
dedup_output=$5
genome_fasta_path=$6
output_directory=$7
cores=$8
deduplicate=$9
"
#echo "mapping command used:"
#echo "$0 -t $read1_trimmed -T $read2_trimmed -g $genome_fasta_path -o $output_directory -c $cores -d $deduplicate -l $level" 
#echo ""

module load bismark/0.22.3

        if [ ! -d "$output_directory" ] 
        	then
                mkdir -p $output_directory
        		echo "Output directoryfor mapping is "$output_directory""
        		#if statement to stopcode from running this section if already done
        else
        		echo "Error: Directory "$output_directory" already exist."
        
        fi
        
        echo "bismark --bam --maxins 800 $genome_fasta_path -1 $read1_input -2 $read2_input -o $output_directory --unmapped --nucleotide_coverage --multicore $cores"
        if [ ! -f $read1_bismark ]; then
		    #expected_mapping_output=${output_directory}/${read1_trimmed##*/}_bismark_bt2_PE_report.txt
	        echo "Expected mapping output file is $read1_bismark"
            echo "$(basename "$read1_input") and $(basename "$read2_input") read pair not yet mapped to "$genome_name" genome. Mapping now."
			bismark --bam --maxins 800 $genome_fasta_path -1 $read1_input -2 $read2_input -o $output_directory --unmapped --nucleotide_coverage --multicore $cores
			echo "finished mapping $(basename $read1_input) and $(basename $read2_input) read pair."
            echo "map to control seqs complete for $genome_name"

        else
		    echo "$(basename "$read1_input") and $(basename "$read2_input") have already been Bismark mapped to $genome_fasta_path"
		fi

	
if [ "$deduplicate" == "TRUE" ] || [ "$deduplicate" == TRUE ] || [ "$deduplicate" == "true" ] || [ "$deduplicate" == "True" ]; then
		echo "Deduplication requested for mapping output of $(basename "$read1_input") and $(basename "$read2_input") to genome $genome_fasta_path"
		#TO DO: this path needs to be checked that it is correct!
		echo "Expected deduplication output file is $dedup_output"
	
        
		if [ ! -f "$dedup_output" ]; then
                cd $output_directory
                deduplication_input=$(basename $dedup_input)
				echo "Begin deduplicating $dedup_input"
            	deduplicate_bismark -p --bam $dedup_input -o $dedup_input
            	echo "Finished deduplicating $dedup_input"
        fi
    
    else        
    	echo "Deduplication NOT requested. ( If this is not desired then set '-d true' in the command $0 )"
fi

#TO DO: This script could use testing code that tests each 'then' and 'else' sections of each 'if' statement! So a total of 4 test cases.

#TO DO: halt code with any error to prevent errors from snowballing and not being detected by end user or by development team. http://web.archive.org/web/20110314180918/http://www.davidpashley.com/articles/writing-robust-shell-scripts.html

#TO DO: 
#should do test on relative and absolute paths and ensure that it is mapping the files. it's looking for the relative path . which exists always. It should map and output the mapped files to the current working directory. 
