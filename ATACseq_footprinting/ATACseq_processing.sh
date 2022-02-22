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
gsize=$2
extsize=$3
shifts=$4
broad=$5
nomodel=$6
blacklist=$7
whitelist=$8
genome_folder=$9

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

########################---STEP 5: CREATE COVERAGE FILES, THEN SORT---##################################
#Create: coverage.bg, coverage.sorted.bg, coverage.bw
#After creating the coverage files, sort them
if [ ! -f "$bam_path/coverage/${PREFIX}_coverage.bg" ]; then
    if [ ! -d "$bam_path/coverage" ]; then
        mkdir "$bam_path/coverage"
    fi

    bedtools genomecov -ibam ${PREFIX}.merged.sorted.bam -g "$genome_folder/chromsizes.txt" -bg  > "$bam_path/coverage/${PREFIX}_coverage.bg"
    echo "creation of coverage file complete"
else
    echo "coverage file already created"
fi

if [ ! -f "$bam_path/coverage/${PREFIX}_coverage.sorted.bg" ]; then
    if [ ! -d "$bam_path/coverage" ]; then
        mkdir "$bam_path/coverage"
    fi

    sort -k1,1 -k2,2n "$bam_path/coverage/${PREFIX}_coverage.bg" > "$bam_path/coverage/${PREFIX}_coverage.sorted.bg"
    echo "creation of coverage file complete"
else
    echo "coverage file already created"
fi

if [ ! -f "$bam_path/coverage/${PREFIX}_coverage.bw" ]; then
    if [ ! -d "$bam_path/coverage" ]; then
        mkdir "$bam_path/coverage"
    fi

    chmod 775 $genome_folder/chromsizes.txt
    bedGraphToBigWig "$bam_path/coverage/${PREFIX}_coverage.sorted.bg" "$genome_folder/chromsizes.txt" "$bam_path/coverage/${PREFIX}_coverage.bw" 
    echo "begraph to BigWig complete"
else
    echo "bedgraph has already been converted to BigWig"
fi

###########################---STEP 6: PEAK CALLING WITH MACS2---########################################
if [ ! -f "$bam_path/peak_calling/${PREFIX}/${PREFIX}_peaks.broadPeak" ]; then

      if [ ! -f "$bam_path/peak_calling" ]; then
            mkdir "$bam_path/peak_calling/"
      fi

      if [ ! -f "$bam_path/peak_calling/${PREFIX}" ]; then
            mkdir "$bam_path/peak_calling/${PREFIX}"
      fi
      
      module load macs2
      echo "Running macs2 with .bam file: ${PREFIX}.merged.sorted.bam"

      if [ $broad == true ] && [ $nomodel == true ]; then
            macs2 callpeak -t "$bam_path/${PREFIX}.merged.sorted.bam" --name ${PREFIX} --outdir "$bam_path/peak_calling/${PREFIX}" --gsize $gsize --nomodel --shift -$shifts --extsize $extsize --broad 
      fi

      if [ $broad == true ] && [ $nomodel == false ]; then
            macs2 callpeak -t "$bam_path/${PREFIX}.merged.sorted.bam" --name ${PREFIX} --outdir "$bam_path/peak_calling/${PREFIX}" --gsize $gsize --shift -$shifts --extsize $extsize --broad 
      fi

      if [ $broad == false ] && [ $nomodel == true ]; then
            macs2 callpeak -t "$bam_path/${PREFIX}.merged.sorted.bam" --name ${PREFIX} --outdir "$bam_path/peak_calling/${PREFIX}" --gsize $gsize --nomodel --shift -$shifts --extsize $extsize
      fi
      
      echo "peak calling with macs2 complete"
else
      echo "peaking calling with macs2 already done"
fi
      
if [ ! -f "$bam_path/peak_calling/${PREFIX}/${PREFIX}_raw.bed" ]; then
      
    cp "$bam_path/peak_calling/${PREFIX}/${PREFIX}_peaks.broadPeak" "$bam_path/peak_calling/${PREFIX}/${PREFIX}_raw.bed"

fi

##########################---STEP 7: REMOVE BLACKLISTED REGIONS---######################################
if [ ! -f "$bam_path/peak_calling/${PREFIX}/${PREFIX}_union_final.bed" ]; then
    cat $bam_path/peak_calling/${PREFIX}/${PREFIX}_raw.bed | cut -f1-3 | sort -k1,1 -k2,2n | bedtools merge -d 5 | \
        bedtools subtract -a - -b $blacklist -A | bedtools intersect -a - -b $whitelist -wa | awk '$1 !~ /[M]/' | \
        sed "s/$/ ${PREFIX}/"  > $bam_path/peak_calling/${PREFIX}/${PREFIX}_union.bed
    #cat $bam_path/peak_calling/${PREFIX}/${PREFIX}_union.bed | tr ' ' '\t' > $bam_path/peak_calling/${PREFIX}/${PREFIX}_union_final.bed
        #excludes mitochondria chromosome (M)
        #adds condition name to each peak
    echo "blacklisted region removal complete"
else
    echo "blacklisted regions have already been removed"
fi
