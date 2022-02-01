#!/bin/bash

#script to sort and index bam file with error checking code wrapped around it. 

#TO DO: test all the paths in this script!

while getopts b:o: flag
do
    case "${flag}" in
        b) bam=${OPTARG};;
        o) output_directory=${OPTARG};;
    esac
done


module load samtools/1.9

#TO DO: add code to check all inputs are present and if not to give error message with help information about all the script inputs

echo ""
echo "sort_and_index command used:"
echo "$0 -b $bam -o $output_directory"
echo ""

expected_sorted_output="${output_directory}/$(basename "$bam" ".bam").sorted.bam"
echo "Expected sorted output file is '$expected_sorted_output'"

if [ ! -f "${output_directory}/$(basename "$bam" ".bam").sorted.bam" ]
	then
    	echo "Starting to sort $(basename "$bam")"
    	#cd "$temp_path/$unmethyl_control/$hydroxymethyl_control/$methyl_control/genome_alignment"
        #echo "$temp_path/$unmethyl_control/$hydroxymethyl_control/$methyl_control/genome_alignment/${trimmed_R1_file_name}.fastq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1_bismark_bt2_pe.deduplicated.bam.sorted.bam does not yet exist"
    	samtools sort $bam -o $expected_sorted_output
    	echo "Finished sorting $(basename "$bam")"
    	echo "Starting to index $(basename "$expected_sorted_output")"
    	samtools index $expected_sorted_output
    	echo "Finished indexing $(basename "$expected_sorted_output")"

    	#rsync -vur --exclude "main_genome" --exclude "unmethyl_genome" --exclude "hydroxymethyl_genome" --exclude "methyl_genome" $temp_path/ $seq_path
    	echo "Sorting and indexing complete for the bam file $bam"
	else
    	echo "Sorted file already exists (and likely an indexed file too) for $bam"
fi





