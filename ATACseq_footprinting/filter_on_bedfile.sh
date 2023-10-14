#!/bin/bash
echo "entering filter_on_bed_file"

output_path=$1
output_temp_dir=$2
parameter_file=$3
PREFIX=$4
bed_file=$5

echo "filter_on_bed_file.sh used the following parameters:
$0 $1 $2 $3 $4 $5"

if [ $SLURM_ARRAY_TASK_ID -eq 1 ]; then
	echo "#####################peak_calling.sh file#####################
        command given to peak_calling.sh $0 $1 $2 $3 $4 $5
        The paramters given:            
            output_path=$1
            output_temp_dir=$2
            parameter_file=$3
            PREFIX=$4
            bed_file=$5
            " >> $parameter_file
fi

###########################---STEP 1: FILTERING WITH BEDTOOLS INTERSECT AND SAMTOOLS INDEX---########################################
if [ ! -f "$output_temp_dir/${PREFIX}.merged.sorted.filtered.bam" ]; then
    module load bedtools
    module load samtools/1.9
    echo "filtering ${PREFIX}.merged.sorted.bam using $bed_file"
    bedtools intersect -abam $output_temp_dir/${PREFIX}.merged.sorted.bam -b $bed_file > $output_temp_dir/${PREFIX}.merged.sorted.filtered.bam
    samtools index $output_temp_dir/${PREFIX}.merged.sorted.filtered.bam $output_temp_dir/${PREFIX}.merged.sorted.filtered.bai
    rsync -vur "$output_temp_dir/" "$output_path"
    echo "filtering of bam file complete"
else
    echo "filtering of bam file already completed"
fi 
