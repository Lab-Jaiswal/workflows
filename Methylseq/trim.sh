#!/bin/bash

#script to trim paired reads with error checking code wrapped around it. 
#Is set to autodetect cores available and use them.
module load bismark/0.22.3
module load samtools
while getopts r:R:t:T: flag
do
    case "${flag}" in
        r) read1=${OPTARG};;
        R) read2=${OPTARG};;
        t) read1_trimmed=${OPTARG};;
        T) read2_trimmed=${OPTARG};;
	o) data_path=${OPTARG};;
    esac
done

echo "entering trim.sh script"
echo "trim command used the following parameters:"
echo "$0 -r $read1 -R $read2 -t $read1_trimmed -T $read2_trimmed -o $data_path"

if [ $SLURM_ARRAY_TASK_ID -eq 1 ]; then
	echo "command given to trim.sh includes the following parameters:
		R1: $read1
		R2: $read2
		read1_trimmed: $read1_trimmed
		read2_trimmed: $read2_trimmed
		data_path: $data_path
		" >> $parameter_file 
fi

module load cutadapt/3.4

parent_directory="$(dirname "$read1_trimmed")"

echo "Parent directory for trimming is '$parent_directory'"

if [ ! -f "$read1_trimmed" ] && [ ! -f "$read2_trimmed" ]; then			
	echo "$(basename "$read1") and $(basename "$read2") read pair not yet trimmed. Trimming now."
	cutadapt -q 20 -m 15 -a AGATCGGAAGAGC -A AAATCAAAAAAAC -o $read1_trimmed -p $read2_trimmed $read1 $read2 -j 0
	echo "finished trimming $(basename $read1) and $(basename $read2) read pair."
	rsync -vur $parent_directory/ $data_path
else
	echo "$(basename "$read1") and $(basename "$read2") have already been trimmed"
fi
