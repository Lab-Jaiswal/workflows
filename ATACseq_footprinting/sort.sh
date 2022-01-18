#!/bin/bash

#SBATCH --job-name=sort
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --partition=batch
#SBATCH --mem=500G
#SBATCH --time=3:00:00
#SBATCH --account=sjaiswal

#####################################---STEP 1: SET UP---############################################### 
bam_path=$1 #set arguments
module load samtools/1.9 #load necessary modules

line_number=$SLURM_ARRAY_TASK_ID #get index of which file to process from $SLURM_ARRAY_TASK_ID provided by SLURM
bam_file="${bam_path}/BAMs" #provide path to file containing list of fastq files
bam_prefix="$(sed "${line_number}q; d" "${bam_file}")" #extract only the line number corresponding to $SLURM_ARRAY_TASK_ID

#Collect the sample name from the BAMs file
#Ex: If $bam_prefix is "/atac_seq/data/ATAC_tet2_KO_LDL", then the PREFIX is "ATAC_tet2_KO_LDL"
PREFIX=$(basename "${bam_prefix}")
echo "filename: $PREFIX"

cd $bam_path

#Find the number of replicates associated with each prefix
number_replicates=$(find . -name "${PREFIX}_Rep[1-9]*_treat.bam" |
            sed "s/${PREFIX}_Rep\([0-9][0-9]*\)\_treat.bam/\1/" |
            sort -n |
            tail -n 1)
reps=$(basename "${number_replicates}")

echo "number of replicates for ${PREFIX}: $reps"

######################---STEP 2: SORT BAMS (STILL SPLIT UP BY REPLICATE)---#############################
if [ ! -f "$bam_path/${PREFIX}_Rep1_treat.sorted.bam" ]; then
      replicates=$(seq $reps)
          for i in ${replicates[@]}
              do
                  rep="${PREFIX}_Rep${i}_treat.bam"
                  bam_split_files="${bam_path}/BAMs_rep_${PREFIX}"
                  rep_sorted="${PREFIX}_Rep${i}_treat.sorted.bam" 
                  echo $rep_sorted >> $bam_split_files

                  if [ ! -f "$bam_path/$rep_sorted" ]; then
                      samtools sort $rep -o $rep_sorted
                  fi
              done
      echo "sorting bams complete"
else
      echo "bams were already sorted"
fi

####################---STEP 3: MERGE REPLICATE BAMS OF THE SAME CONDITION---############################
#Only complete this step for conditions that have more than 1 replicate

if [ ! -f "$bam_path/${PREFIX}.merged.bam" ] && [ $reps -gt 1 ]; then
      samtools merge -b "${bam_path}/BAMs_rep_${PREFIX}" "${PREFIX}.merged.bam"
      echo "merging of sorted bams complete"
else
      echo "sorted bams were already merged"
      cp "$bam_path/${PREFIX}_Rep1_treat.sorted.bam" "$bam_path/${PREFIX}.merged.bam"
fi

##############################---STEP 4: INDEX MERGED BAMS---###########################################
#Result should be 1 final sorted bam for each experimental condition
#In order to index the merged bams, you must sort them first

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

##############################---STEP 6: INDEX MERGED BAMS---###########################################

gsize=2620345972 #is this correct?!
macs="$bam_path/peak_calling/${PREFIX}/${PREFIX}_peaks.broadPeak"
raw="$bam_path/peak_calling/${PREFIX}/${PREFIX}_raw.bed"
log="$bam_path/Logs/${PREFIX}_peak_calling.log"

echo "Running macs2 with .bam file: ${PREFIX}.merged.sorted.bam"

#Arguments:
#-t: the IP data file
#--outdir: the output directory

mkdir "$bam_path/peak_calling/${PREFIX}"

macs2 callpeak -t "$bam_path/${PREFIX}.merged.sorted.bam" --name ${PREFIX} --outdir "$bam_path/peak_calling/${PREFIX}" --gsize $gsize --nomodel --shift -100 --extsize 200 --broad &> "$bam_path/Logs/{PREFIX}_peak_calling.log"
cp "$bam_path/peak_calling/${PREFIX}/{PREFIX}_peaks.broadPeak" "$bam_path/peak_calling/${PREFIX}/{PREFIX}_raw.bed"


