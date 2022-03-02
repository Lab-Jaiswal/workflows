#!/bin/bash

echo "entering coverage_file_creation.sh script"

output_path=$1
output_temp_dir=$2
parameter_file=$3
PREFIX=$4

echo "coverage_file_creation.sh used the following parameters:
$0 $1 $2 $3 $4"

if [ $SLURM_ARRAY_TASK_ID -eq 1 ]; then
	echo "#####################coverage_file_creation.sh file#####################
        coverage_file_creation.sh entered                
        command given to coverage_file_creation.sh $0 $1 $2 $3 $4
        The paramters given:
            output_path=$1
            output_temp_dir=$2
            parameter_file=$3
            PREFIX=$4
        " >> $parameter_file
fi

########################---STEP 5: CREATE COVERAGE FILES, THEN SORT---##################################
#Create: coverage.bg, coverage.sorted.bg, coverage.bw
#After creating the coverage files, sort them
if [ ! -f "$output_dir/coverage/${PREFIX}_coverage.bg" ]; then
    if [ ! -d "$output_dir/coverage" ]; then
        mkdir "$output_dir/coverage"
    fi
    
    module load samtools/1.9
    bedtools genomecov -ibam $output_dir/${PREFIX}.merged.sorted.bam -g "$genome_folder/chromsizes.txt" -bg  > "$output_dir/coverage/${PREFIX}_coverage.bg"
    echo "creation of coverage file complete"
else
    echo "coverage file already created"
fi

if [ ! -f "$output_dir/coverage/${PREFIX}_coverage.sorted.bg" ]; then
    if [ ! -d "$output_dir/coverage" ]; then
        mkdir "$output_dir/coverage"
    fi

    sort -k1,1 -k2,2n "$output_dir/coverage/${PREFIX}_coverage.bg" > "$output_dir/coverage/${PREFIX}_coverage.sorted.bg"
    echo "creation of coverage file complete"
else
    echo "coverage file already created"
fi

if [ ! -f "$output_dir/coverage/${PREFIX}_coverage.bw" ]; then
    if [ ! -d "$output_dir/coverage" ]; then
        mkdir "$output_dir/coverage"
    fi

    #chmod 775 $genome_folder/chromsizes.txt
    bedGraphToBigWig "$output_dir/coverage/${PREFIX}_coverage.sorted.bg" "$genome_folder/chromsizes.txt" "$output_dir/coverage/${PREFIX}_coverage.bw" 
    echo "begraph to BigWig complete"
else
    echo "bedgraph has already been converted to BigWig"
fi

rsync -vur "$output_dir/" "$output_path"

if [ $SLURM_ARRAY_TASK_ID -eq 1 ]; then
	echo "coverage_file_creation.sh script complete
        " >> $parameter_file
fi
