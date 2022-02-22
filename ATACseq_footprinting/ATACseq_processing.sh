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
output_path=$2
gsize=$3
extsize=$4
shifts=$5
broad=$6
nomodel=$7
blacklist=$8
whitelist=$9
genome_folder=${10}
parameter_file=${11}

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

temp_path=$(mktemp -d /tmp/tmp.XXXXXXXXXX)
echo "temp_path is: " $temp_path
echo "copying bams from the data path..."
rsync -vur "$bam_path/" $temp_path
rsync -vur "$output_path/" "$temp_path/output_path"

output_dir="$temp_path/output_path"

cd $temp_path

######################---STEP 2: SORT BAMS (STILL SPLIT UP BY REPLICATE)---#############################
if [ ! -f "$output_dir/${PREFIX}_Rep1_treat.sorted.bam" ]; then
      replicates=$(seq $reps)
          for i in ${replicates[@]}
          do
                  rep="${PREFIX}_Rep${i}_treat.bam"
                  bam_split_files="${output_dir}/BAMs_rep_${PREFIX}"
                  rep_sorted="${output_dir}/${PREFIX}_Rep${i}_treat.sorted.bam" 
                  echo $rep_sorted >> $bam_split_files

                  if [ ! -f "$rep_sorted" ]; then
                      samtools sort $rep -o $rep_sorted
                  fi
              done
      echo "sorting bams complete"
else
      echo "bams were already sorted"
fi

rsync -vur "$output_dir/" "$output_path"

####################---STEP 3: MERGE REPLICATE BAMS OF THE SAME CONDITION---############################
#Only complete this step for conditions that have more than 1 replicate

if [ ! -f "$output_dir/${PREFIX}.merged.bam" ] && [ $reps -gt 1 ]; then
      samtools merge -b "${output_dir}/BAMs_rep_${PREFIX}" "${output_dir}/${PREFIX}.merged.bam"
      echo "merging of sorted bams complete"
  else
      echo "sorted bams were already merged"
      cp "$output_dir/${PREFIX}_Rep1_treat.sorted.bam" "$output_dir/${PREFIX}.merged.bam"
fi

rsync -vur "$output_dir/" "$output_path"

##############################---STEP 4: INDEX MERGED BAMS---###########################################
#Result should be 1 final sorted bam for each experimental condition
#In order to index the merged bams, you must sort them first

if [ ! -f "$output_dir/${PREFIX}.merged.sorted.bam" ]; then
      samtools sort "$output_dir/${PREFIX}.merged.bam" -o "$output_dir/${PREFIX}.merged.sorted.bam"
      echo "sorting of merged bams complete"
else
      echo "merged bams already sorted"
fi     

if [ ! -f "$output_dir/${PREFIX}.merged.sorted.bai" ]; then
      samtools index "$output_dir/${PREFIX}.merged.sorted.bam" "$output_dir/${PREFIX}.merged.sorted.bai"
      echo "indexing of sorted merged bams complete"
else
      echo "sorted merged bams already indexed"
fi      

rsync -vur "$output_dir/" "$output_path"

########################---STEP 5: CREATE COVERAGE FILES, THEN SORT---##################################
#Create: coverage.bg, coverage.sorted.bg, coverage.bw
#After creating the coverage files, sort them
if [ ! -f "$output_dir/coverage/${PREFIX}_coverage.bg" ]; then
    if [ ! -d "$output_dir/coverage" ]; then
        mkdir "$output_dir/coverage"
    fi

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

###########################---STEP 6: PEAK CALLING WITH MACS2---########################################
if [ ! -f "$output_dir/peak_calling/${PREFIX}/${PREFIX}_peaks.broadPeak" ]; then

      if [ ! -f "$output_dir/peak_calling" ]; then
            mkdir "$output_dir/peak_calling/"
      fi

      if [ ! -f "$output_dir/peak_calling/${PREFIX}" ]; then
            mkdir "$output_dir/peak_calling/${PREFIX}"
      fi
      
      module load macs2
      echo "Running macs2 with .bam file: $output_dir/${PREFIX}.merged.sorted.bam"

      if [ $broad == true ] && [ $nomodel == true ]; then
            macs2 callpeak -t "$output_dir/${PREFIX}.merged.sorted.bam" --name ${PREFIX} --outdir "$output_dir/peak_calling/${PREFIX}" --gsize $gsize --nomodel --shift -$shifts --extsize $extsize --broad 
      fi

      if [ $broad == true ] && [ $nomodel == false ]; then
            macs2 callpeak -t "$output_dir/${PREFIX}.merged.sorted.bam" --name ${PREFIX} --outdir "$output_dir/peak_calling/${PREFIX}" --gsize $gsize --shift -$shifts --extsize $extsize --broad 
      fi

      if [ $broad == false ] && [ $nomodel == true ]; then
            macs2 callpeak -t "$output_dir/${PREFIX}.merged.sorted.bam" --name ${PREFIX} --outdir "$output_dir/peak_calling/${PREFIX}" --gsize $gsize --nomodel --shift -$shifts --extsize $extsize
      fi
      
      echo "peak calling with macs2 complete"
else
      echo "peaking calling with macs2 already done"
fi
      
if [ ! -f "$output_dir/peak_calling/${PREFIX}/${PREFIX}_raw.bed" ]; then
      
    cp "$output_dir/peak_calling/${PREFIX}/${PREFIX}_peaks.broadPeak" "$output_dir/peak_calling/${PREFIX}/${PREFIX}_raw.bed"

fi

rsync -vur "$output_dir/" "$output_path"

##########################---STEP 7: REMOVE BLACKLISTED REGIONS---######################################
if [ ! -f "$output_dir/peak_calling/${PREFIX}/${PREFIX}_union_final.bed" ]; then
    cat $output_dir/peak_calling/${PREFIX}/${PREFIX}_raw.bed | cut -f1-3 | sort -k1,1 -k2,2n | bedtools merge -d 5 | \
        bedtools subtract -a - -b $blacklist -A | bedtools intersect -a - -b $whitelist -wa | awk '$1 !~ /[M]/' | \
        sed "s/$/ ${PREFIX}/"  > $output_dir/peak_calling/${PREFIX}/${PREFIX}_union.bed
    #cat $bam_path/peak_calling/${PREFIX}/${PREFIX}_union.bed | tr ' ' '\t' > $bam_path/peak_calling/${PREFIX}/${PREFIX}_union_final.bed
        #excludes mitochondria chromosome (M)
        #adds condition name to each peak
    echo "blacklisted region removal complete"
else
    echo "blacklisted regions have already been removed"
fi

rsync -vur "$output_dir/" "$output_path"
