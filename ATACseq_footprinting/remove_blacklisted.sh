#!/bin/bash

echo "entering remove_blacklisted.sh"

output_path=$1
output_temp_dir=$2
blacklist=$3
whitelist=$4
parameter_file=$5
PREFIX=$6

echo "remove_blacklisted.sh used the following parameters:
$0 $1 $2 $3 $4 $5 $6"

if [ $SLURM_ARRAY_TASK_ID -eq 1 ]; then
	echo "#####################remove_blackisted.sh file#####################
        command given to remove_blacklisted.sh: $0 $1 $2 $3 $4 $5 $6
        The paramters given:
            output_path=$1
            output_temp_dir=$2
            blacklist=$3
            whitelist=$4
            parameter_file=$5
            PREFIX=$6
          " >> $parameter_file
fi
##########################---STEP 7: REMOVE BLACKLISTED REGIONS---######################################
if [ ! -f "$output_temp_dir/peak_calling/${PREFIX}/${PREFIX}_union_final.bed" ]; then
    module load samtools/1.9
    cat $output_temp_dir/peak_calling/${PREFIX}/${PREFIX}_raw.bed | cut -f1-3 | sort -k1,1 -k2,2n | bedtools merge -d 5 | \
        bedtools subtract -a - -b $blacklist -A | bedtools intersect -a - -b $whitelist -wa | awk '$1 !~ /[M]/' | \
        sed "s/$/ ${PREFIX}/"  > $output_temp_dir/peak_calling/${PREFIX}/${PREFIX}_union.bed
    echo "blacklisted region removal complete"
    rsync -vur "$output_temp_dir/" "$output_path"


else
    echo "blacklisted regions have already been removed"
fi

if [ $SLURM_ARRAY_TASK_ID -eq 1 ]; then
	echo "remove_blacklisted.sh script complete
        " >> $parameter_file
fi
