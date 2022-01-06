#!/bin/bash

if [ -z $1 ] || [ -z $2 ] || [ -z $3 ] || [ -z $4 ] || [ -z $5 ] || [ -z $6 ] || [ -z $7 ] || [ -z $8 ] || [ -z $9 ]; then
    echo "Format: ./submit_methylseq.sh [data_directory] [unmethyl_control] [unmethyl_control_fasta] [hydroxymethyl_control] [hydroxymethyl_control_fasta] [methyl_control] [methyl_control_fasta] [genome_path] [phix_path]"
    exit 1
else

    TEMP=`getopt -o vdm: --long force  -n 'submit_methylseq.sh' -- "$@"`
        eval set -- "$TEMP"
        
        force=false
                        
    while true; do
        case "$1" in
            -f | --force ) force=true; shift 2 ;;
            -- ) shift; break ;;
            * ) break ;;
        esac
    done

    echo "force: $force"
    #data_path=/oak/stanford/groups/sjaiswal/kameronr/JG97/methylseq/miseq/214460250
    data_path=$1
    unmethyl_control=$2
    unmethyl_control_fasta=$3
    hydroxymethyl_control=$4
    hydroxymethyl_control_fasta=$5
    methyl_control=$6
    methyl_control_fasta=$7
    genome_path=$8
    phix_path=$9
    temp_path=$10

    code_directory=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )

    #cores=$(expr ${SLURM_CPUS_PER_TASK}/2)
    cores=24
    echo core = $cores

    echo "$data_path"
    fastq_file="${data_path}/fastq/FASTQs" #give a path to a file to store the paths to the fastq files in $fastq_directory
    echo "$fastq_file"
    find "$data_path/fastq" -type f | grep ".*\.fastq.gz$" | grep -v ".*\.trimmed.fastq.gz$" | sed -e 's/_R1.*$//g' | sed -e 's/_R2.*$//g' | sort -u > "${fastq_file}" #generate list of full paths to fastq files and save to the file in $fastq_list
    #tail -n +2 "$fastq_file" > "$fastq_file.tmp" && mv "$fastq_file.tmp" "$fastq_file"
    number_fastq=$(wc -l < "${fastq_file}") #get the number of files
    echo "$number_fastq"
    array_length_trim=$(( $number_fastq * 2 ))
    array_length="$number_fastq"

    if [ -d "$data_path/fastq/$unmethyl_control" ];then
        mkdir "$data_path/fastq/$unmethyl_control"
    fi

    if [ -d "$data_path/fastq/$unmethyl_control/$hydroxymethyl_control" ]; then
        mkdir "$data_path/fastq/$unmethyl_control/$hydroxymethyl_control"
    fi

    if [ -d "$data_path/fastq/$unmethyl_control/$hydroxymethyl_control/$hydroxymethyl_control" ]; then
        mkdir "$data_path/fastq/$unmethyl_control/$hydroxymethyl_control/$methyl_control"
    fi

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

    picard=$(find "$data_path/fastq/$unmethyl_control/$hydroxymethyl_control/$methyl_control/genome_alignment" -type f | grep ".*\.bam_picard_insert_size_plot.pdf$" | sort -u | wc -l)
    echo "picard: $picard"

    if [[ $picard -lt 1 ]]; then
        #took out force here - add it back in
            sbatch -o "${data_path}/Logs/%A_%a.log" `#put into log` \
                    -a "1-${array_length}" `#initiate job array equal to the number of fastq files` \
                   -W `#indicates to the script not to move on until the sbatch operation is complete` \
                    "${code_directory}/methylseq.sh" \
                    $data_path $unmethyl_control_fasta $unmethyl_control $hydroxymethyl_control_fasta $hydroxymethyl_control $methyl_control_fasta $methyl_control $genome_path $phix_path $cores $temp_path
            wait
        else
            echo "picard files found"
    fi


    bam_file="$data_path/fastq/$unmethyl_control/$hydroxymethyl_control/$methyl_control/genome_alignment/bam_files" 
    find "$data_path/fastq/$unmethyl_control/$hydroxymethyl_control/$methyl_control/genome_alignment/split_bams" -type f | grep ".*\.bam$" | sort -u > "${bam_file}" #generate list of full paths to fastq files and save to the file in $fastq_list
    #tail -n +2 "$fastq_file" > "$fastq_file.tmp" && mv "$fastq_file.tmp" "$fastq_file"
    number_bam=$(wc -l < "${bam_file}") #get the number of files
    echo "number_bam: $number_bam"

    bedgraph=$(find "$data_path/fastq/$unmethyl_control/$hydroxymethyl_control/$methyl_control/genome_alignment/split_bams" -type f | grep ".*\.bedGraph.gz$" | sort -u | wc -l)
    echo "bedgraph: $bedgraph"

    if [[ $bedgraph -lt 2 ]] || [[ $force = true ]]; then
            sbatch -o "${data_path}/Logs/%A_%a.log" `#put into log` \
                    -a "1-${number_bam}" `#initiate job array equal to the number of fastq files` \
                    -W `#indicates to the script not to move on until the sbatch operation is complete` \
                    "${code_directory}/extract_methylation.sh" \
                    $data_path $unmethyl_control_fasta $unmethyl_control $hydroxymethyl_control_fasta $hydroxymethyl_control $methyl_control_fasta $methyl_control $genome_path $phix_path $cores $temp_path
                wait
            else
                echo "split coverage files and bed graphs have already been created for the genome alignemnt"
    fi

    bedgraph=$(find "$data_path/fastq/$unmethyl_control/$hydroxymethyl_control/$methyl_control/genome_alignment" -maxdepth 1 -type f  | grep ".*\.bedGraph.gz$" | sort -u | wc -l)
    echo "bedgraph: $bedgraph"

    if [[ $bedgraph -lt 1 ]] || [[ $force = true ]] ; then
            echo "combining bed graphs and .cov files currently split by chromosome"
            sbatch -o "${data_path}/Logs/%A_%a.log" `#put into log` \
                    -a "1-${array_length}" `#initiate job array equal to the number of fastq files` \
                    -W `#indicates to the script not to move on until the sbatch operation is complete` \
                    "${code_directory}/join_coverage.sh" \
                    $data_path  $unmethyl_control  $hydroxymethyl_control $methyl_control  $cores $temp_path
                wait
        else
            echo "combined bed graphs and coverage files have already been created for the genome alignemnt"
    fi

fi
