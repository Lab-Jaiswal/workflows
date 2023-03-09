#!/bin/bash

#script to trim paired reads with error checking code wrapped around it. 
#Is set to autodetect cores available and use them.
module load bismark/0.22.3
module load samtools

read1=$1
read2=$2
read1_trimmed=$3
read2_trimmed=$4
data_path=$5
parameter_file=$6

echo "entering trim.sh script"
echo "trim command used the following parameters:"
echo "$0 $1 $2 $3 $4 $5 $6"

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

if [[ ! -f "$read1_trimmed" ]]; then			
	echo "$(basename "$read1") and $(basename "$read2") read pair not yet trimmed. Trimming now."
	cutadapt -q 20 -m 15 -a AGATCGGAAGAGC -A AAATCAAAAAAAC -o $read1_trimmed -p $read2_trimmed $read1 $read2 -j 0
	echo "finished trimming $(basename $read1) and $(basename $read2) read pair."
	rsync -vur $parent_directory/ "$data_path/fastq"
else
	echo "$(basename "$read1") and $(basename "$read2") have already been trimmed"
fi

#TO DO: 
#for EMseq, recommended to remove 10bp from ends to avoid bias and increase mapping rates. 
#http://felixkrueger.github.io/Bismark/bismark/library_types/#em-seq-neb
#https://github.com/FelixKrueger/Bismark/issues/509

