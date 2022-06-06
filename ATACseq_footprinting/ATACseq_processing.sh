#!/bin/bash

#SBATCH --job-name=sort
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=batch
#SBATCH --mem=128G
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
code_directory=${12}
bed_file=${13}
filter=${14}

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

temp_path=$(mktemp -d)
echo "temp_path is: " $temp_path
echo "copying bams from the data path..."
rsync -vur "$bam_path/" $temp_path
rsync -vur "$output_path/" "$temp_path/output_path"

output_temp_dir="$temp_path/output_path"

cd $temp_path

######################---STEP 2: SORT BAMS (STILL SPLIT UP BY REPLICATE)---#############################
####################---STEP 3: MERGE REPLICATE BAMS OF THE SAME CONDITION---############################
##############################---STEP 4: INDEX MERGED BAMS---###########################################

$code_directory/sort_merge_index.sh  $output_path $output_temp_dir $reps $parameter_file $PREFIX 

########################---STEP 5: CREATE COVERAGE FILES, THEN SORT---##################################

 $code_directory/coverage_file_creation.sh $output_path $output_temp_dir $genome_folder $parameter_file $PREFIX 

###########################---STEP 6: PEAK CALLING WITH MACS2---########################################

$code_directory/peak_calling.sh $output_path $output_temp_dir $gsize $extsize $shifts $broad $nomodel $whitelist $parameter_file $PREFIX 

##########################---STEP 7: REMOVE BLACKLISTED REGIONS---######################################

$code_directory/remove_blacklisted.sh $output_path $output_temp_dir $blacklist $whitelist $parameter_file $PREFIX 

#################################---STEP 7: FILTER BEDFILE---###########################################
if [ $filter -ne 0 ]; then
            echo "bed_file: $bed_file"
            $code_directory/filter_on_bedfile.sh $output_path $output_temp_dir $parameter_file $PREFIX $bed_file
else
            echo "filtering of bam files via bed file was not requested. 
            Please use argument "--filter bed_file_location" if this is not correct.
            " >> $parameter_file
fi
