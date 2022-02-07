#!/bin/bash

#script to sort and index bam file with error checking code wrapped around it. 

#TO DO: test all the paths in this script!

sort_input=$1
index_input=$2
index_output=$3
output_directory=$4

module load samtools/1.9


#TO DO: add code to check all inputs are present and if not to give error message with help information about all the script inputs

echo "sort_and_index command used:"
echo "$0 -b $bam -o $output_directory"
echo ""

echo "Expected sorted output file is '$index_input'"

if [ ! -f "$index_input" ]; then
    echo "Starting to sort $(basename "$sort_input")"
    samtools sort $sort_input -o $index_input
    echo "sorting of $sort_input is complete"
else
    echo "sorting of $(basename "$sort_input") already complete"
fi

if  [ -f "$index_output" ]; then
    echo "Starting to index $(basename "$index_input")"  
    samtools index $index_input
    echo "indexing of $(basename "$index_input") is complete"
else
    echo "sorting of $(basename "$index_input") already complete"
fi

    	
    	



