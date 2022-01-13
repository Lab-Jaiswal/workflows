#!/bin/bash

#SBATCH --job-name=sort
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --partition=batch
#SBATCH --mem=500G
#SBATCH --time=3:00:00
#SBATCH --account=sjaiswal
#cd /oak/stanford/groups/sjaiswal/kameronr/JG97/methylseq/miseq/213566405/FASTQ/

bam_path=$1

#cd $temp_path

module load samtools/1.9

line_number=$SLURM_ARRAY_TASK_ID #get index of which file to process from $SLURM_ARRAY_TASK_ID provided by SLURM
bam_file="${bam_path}/BAMs" #provide path to file containing list of fastq files
bam_prefix="$(sed "${line_number}q; d" "${bam_file}")" #extract only the line number corresponding to $SLURM_ARRAY_TASK_ID

FILENAME=$(basename "${bam_prefix}")
PREFIX=$FILENAME
echo "filename: $PREFIX"

cd $bam_path
number_replicates=$(find . -name "${PREFIX}_Rep[1-9]*_treat_rep1.bam" |
            sed "s/${PREFIX}_Rep\([0-9][0-9]*\)\_treat_rep1.bam/\1/" |
            sort -n |
            tail -n 1)
reps=$(basename "${number_replicates}")

echo "number of replicates for ${PREFIX}: $reps"

if [ ! -f "$bam_path/${PREFIX}_Rep1_treat_rep1.sorted.bam" ]; then
      replicates=$(seq $reps)
          for i in ${replicates[@]}
              do
                  rep="${PREFIX}_Rep${i}_treat_rep1.bam"
                  bam_split_files="${bam_path}/BAMs_rep_${PREFIX}"
                  rep_sorted="${PREFIX}_Rep${i}_treat_rep1.sorted.bam" 
                  echo $rep_sorted >> $bam_split_files

                  if [ ! -f "$bam_path/$rep_sorted" ]; then
                      samtools sort $rep -o $rep_sorted
                  fi
              done
      echo "sorting bams complete"
else
      echo "bams were already sorted"
fi

#1b. merge bam files with the same condition (for any conditions that have more than 1 replicate)

if [ ! -f "$bam_path/${PREFIX}.merged.bam" ] && [ $reps -gt 1 ]; then
      samtools merge -b "${bam_path}/BAMs_rep_${PREFIX}" "${PREFIX}.merged.bam"
      echo "merging of sorted bams complete"
else
      echo "sorted bams were already merged"
      cp "$bam_path/${PREFIX}_Rep1_treat_rep1.sorted.bam" "$bam_path/${PREFIX}.merged.bam"
fi

#1c. index all resulting files (should be 1 final sorted bam for each experimental condition)

if [ ! -f "$bam_path/${PREFIX}.merged.sorted.bam" ]; then
      samtools sort "${PREFIX}.merged.bam" -o "${PREFIX}.merged.sorted.bam"
      echo "sorting of merged bams complete"
else
      echo "merged bams already sorted"
fi     

if [ ! -f "$bam_path/${PREFIX}.merged.sorted.bai" ]; then
      samtools index "${PREFIX}.merged.sorted.bam" "${PREFIX}.merged.sorted.bai"
      echo "indexing of sorted merged bams complete"
else
      echo "sorted merged bams already indexed"
fi      

