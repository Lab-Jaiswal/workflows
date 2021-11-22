#!/bin/bash

#data_path=/oak/stanford/groups/sjaiswal/kameronr/JG97/methylseq/miseq/214460250
data_path=$1
unmethyl_control=$2
unmethyl_control_fasta=$3
hydroxymethyl_control=$4
hydroxymethyl_control_fasta=$5
genome_path=$6
phix_path=$7
temp_path=$8

code_directory=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )

#cores=$(expr ${SLURM_CPUS_PER_TASK}/2)
cores=24
echo core = $cores

echo "$data_path"
fastq_file="${data_path}/fastq/FASTQs" #give a path to a file to store the paths to the fastq files in $fastq_directory
echo "$fastq_file"
find "${data_path}/" -maxdepth 2 -type f | grep fastq | grep -v "Undetermined" | sed -e 's/_R1.*$//g' | sed -e 's/_R2.*$//g' | sort -u > "${fastq_file}" #generate list of full paths to fastq files and save to the file in $fastq_list
tail -n +2 "$fastq_file" > "$fastq_file.tmp" && mv "$fastq_file.tmp" "$fastq_file"
number_fastq=$(wc -l < "${fastq_file}") #get the number of files
echo "$number_fastq"
array_length_trim=$(( $number_fastq * 2 ))
array_length="$number_fastq"

if ! [ -d "$data_path/fastq" ]; then
    cd $data_path
    module load bcl2fastq2
    bcl2fastq -o ./fastq -p 8
    cd $code_directory
else
    echo "bcls already transformed into fastqs"
fi

if ! [ -d "$data_path/Logs" ]; then
    mkdir -p "$data_path/Logs"
fi

trimmed=$(find "$data_path/fastq" -type f | grep "trimmed" | sort -u | wc -l)

if [ $trimmed -le 1 ]; then
         sbatch -o "${data_path}/Logs/%A_%a.log" `#put into log` \
        -a "1-${array_length_trim}" `#initiate job array equal to the number of fastq files` \
        -W `#indicates to the script not to move on until the sbatch operation is complete` \
            "${code_directory}/trim.sh" \
            $data_path $temp_path
        
        wait
    else
        echo "trimmed files found"
fi

sbatch -o "${data_path}/Logs/%A_%a.log" `#put into log` \
        -a "1-${array_length}" `#initiate job array equal to the number of fastq files` \
        -W `#indicates to the script not to move on until the sbatch operation is complete` \
        "${code_directory}/methylseq.sh" \
        $data_path $unmethyl_control_fasta $unmethyl_control $hydroxymethyl_control_fasta $hydroxymethyl_control $genome_path $phix_path $cores
        
    wait


