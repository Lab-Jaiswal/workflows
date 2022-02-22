#!/bin/bash
echo "entering sort_merge_index.sh"

output_path=$1
output_temp_dir=$2
reps=$3
parameter_file=$4
PREFIX=$5

echo "sort_merge_index.sh used the following parameters:
$0 $1 $2 $3 $4 $5"

if [ $SLURM_ARRAY_TASK_ID -eq 1 ]; then
	echo "#####################sort_merge_index.sh file#####################
            sort_merge_index.sh entered
            command given to sort_merge_index.sh $0 $1 $2 $3 $4 $5
            The paramters given:
                output_path=$1
                output_temp_dir=$2
                reps=$3
                parameter_file=$4
                PREFIX=$5
            " >> $parameter_file
fi

######################---STEP 2: SORT BAMS (STILL SPLIT UP BY REPLICATE)---#############################
if [ ! -f "$output_temp_dir/${PREFIX}_Rep1_treat.sorted.bam" ]; then
      module load samtools/1.9
      replicates=$(seq $reps)
          for i in ${replicates[@]}
          do
                  rep="${PREFIX}_Rep${i}_treat.bam"
                  bam_split_files="${output_temp_dir}/BAMs_rep_${PREFIX}"
                  rep_sorted="${output_temp_dir}/${PREFIX}_Rep${i}_treat.sorted.bam" 
                  echo $rep_sorted >> $bam_split_files

                  if [ ! -f "$rep_sorted" ]; then
                      samtools sort $rep -o $rep_sorted
                  fi
              done
      echo "sorting bams complete"
else
      echo "bams were already sorted"
fi

####################---STEP 3: MERGE REPLICATE BAMS OF THE SAME CONDITION---############################
#Only complete this step for conditions that have more than 1 replicate

if [ ! -f "$output_temp_dir/${PREFIX}.merged.bam" ] && [ $reps -gt 1 ]; then
      samtools merge -b "$output_temp_dir/BAMs_rep_${PREFIX}" "$output_temp_dir/${PREFIX}.merged.bam"
      echo "merging of sorted bams complete"
  else
      echo "sorted bams were already merged"
      cp "$output_temp_dir/${PREFIX}_Rep1_treat.sorted.bam" "$output_temp_dir/${PREFIX}.merged.bam"
fi

##############################---STEP 4: INDEX MERGED BAMS---###########################################
#Result should be 1 final sorted bam for each experimental condition
#In order to index the merged bams, you must sort them first

if [ ! -f "$output_temp_dir/${PREFIX}.merged.sorted.bam" ]; then
      samtools sort "$output_temp_dir/${PREFIX}.merged.bam" -o "$output_temp_dir/${PREFIX}.merged.sorted.bam"
      echo "sorting of merged bams complete"
else
      echo "merged bams already sorted"
fi     

if [ ! -f "$output_temp_dir/${PREFIX}.merged.sorted.bai" ]; then
      samtools index "$output_temp_dir/${PREFIX}.merged.sorted.bam" "$output_temp_dir/${PREFIX}.merged.sorted.bai"
      echo "indexing of sorted merged bams complete"
else
      echo "sorted merged bams already indexed"
fi      

rsync -vur "$output_temp_dir/" "$output_path"

if [ $SLURM_ARRAY_TASK_ID -eq 1 ]; then
	echo "sort_merge_index.sh script complete
        " >> $parameter_file
fi
