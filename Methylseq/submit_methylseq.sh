#!/bin/bash
##################################################################################################################################
#############################################---STEP 1: SET UP PARAMETERS---###################################################### 
##################################################################################################################################
if [ -z $1 ] || [ -z $2 ] ; then
    echo "Format: ./submit_methylseq.sh [data_directory] [genetic_locations_file]"
    echo "genetic_locations_file is a txt file"
    echo "user can use argument --force to force the code to disregard the checks preventing methylseq.sh, extract_methylation.sh & join_coverage.sh from running"
    echo "user can use argument --cores to select a number of cores different from the default (24)"
    exit 1
else

    TEMP=`getopt -o vdm: --long force,cores  -n 'submit_methylseq.sh' -- "$@"`
        eval set -- "$TEMP"
        
        force=false
        cores=24
                        
    while true; do
        case "$1" in
            -f | --force ) force=true; shift 2 ;;
            --cores ) cores="$2"; shift 2;;
            -- ) shift; break ;;
            * ) break ;;
        esac
    done

    echo "force: $force"
    data_path=$1
    genetic_locations=$2
    unmethyl_control=$(sed -n '1p' $genetic_locations)
    unmethyl_control_fasta=$(sed -n '2p' $genetic_locations)
    hydroxymethyl_control=$(sed -n '3p' $genetic_locations)
    hydroxymethyl_control_fasta=$(sed -n '4p' $genetic_locations)
    methyl_control=$(sed -n '5p' $genetic_locations)
    methyl_control_fasta=$(sed -n '6p' $genetic_locations)
    genome_path=$(sed -n '8p' $genetic_locations)
    phix_path=$(sed -n '10p' $genetic_locations)

    echo "parameters:
    unmethyl_control: $unmethyl_control
    unmethyl_control_fasta: $unmethyl_control_fasta
    hydroxymethyl_control: $hydroxymethyl_control
    hydroxymethyl_control_fasta: $hydroxymethyl_control_fasta
    methyl_control: $methyl_control
    methyl_control_fasta: $methyl_control_fasta
    genome_path: $genome_path
    phix_path: $phix_path"

    code_directory=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )

##################################################################################################################################
#######################################---STEP 2: CREATE NECESSARY FOLDERS---#####################################################
##################################################################################################################################
    
    if [ ! -d "$data_path/fastq/$unmethyl_control" ];then
        mkdir "$data_path/fastq/$unmethyl_control"
    fi

    if [ ! -d "$data_path/fastq/$unmethyl_control/$hydroxymethyl_control" ]; then
        mkdir "$data_path/fastq/$unmethyl_control/$hydroxymethyl_control"
    fi

    if [ ! -d "$data_path/fastq/$unmethyl_control/$hydroxymethyl_control/$methyl_control" ]; then
        mkdir "$data_path/fastq/$unmethyl_control/$hydroxymethyl_control/$methyl_control"
    fi

    if [ ! -d "$data_path/fastq/$unmethyl_control/$hydroxymethyl_control/$methyl_control/genome_alignment" ]; then
        mkdir "$data_path/fastq/$unmethyl_control/$hydroxymethyl_control/$methyl_control/genome_alignment"
    fi

    if ! [ ! -d "$data_path/Logs" ]; then
        mkdir -p "$data_path/Logs"
    fi

##################################################################################################################################
#######################################---STEP 3: CREATE PARAMETER LOG---#####################################################
##################################################################################################################################
    parameter_file="${data_path}/Logs/methylSeqParameters.txt" #give a path to a file to store the paths to the fastq files in $fastq_directory
    
    if [ $force == true ] && [ $cores -eq 24 ]; then 
        echo "call made to execute code: ./submit_methylseq $data_path $genetic_locations --force
        " > $parameter_file
    fi
    
    if [ $force != true ] && [ $cores -ne 24 ]; then 
        echo "call made to execute code: ./submit_methylseq $data_path $genetic_locations --cores $cores
        " > $parameter_file
    fi
    
     if [ $force == true ] && [ $cores -ne 24 ]; then 
        echo "call made to execute code: ./submit_methylseq $data_path $genetic_locations --force --cores $cores
        " > $parameter_file
    fi
    
     if [ $force != true ] && [ $cores -eq 24 ]; then 
        echo "call made to execute code: ./submit_methylseq $data_path $genetic_locations
        " > $parameter_file
    fi
        
    echo "location of file with genome directions: $genetic_locations
    " >> $parameter_file
    echo "contents of file with genome directions:
    $unmethyl_control
    $unmethyl_control_fasta
    $hydroxymethyl_control
    $hydroxymethyl_control_fasta
    $methyl_control
    $methyl_control_fasta
    $main_genome_name
    $genome_path
    $phix_genome_name
    $phix_path
    " >> $parameter_file
    echo "parameters:
    unmethyl_control: $unmethyl_control
    unmethyl_control_fasta: $unmethyl_control_fasta
    hydroxymethyl_control: $hydroxymethyl_control
    hydroxymethyl_control_fasta: $hydroxymethyl_control_fasta
    methyl_control: $methyl_control
    methyl_control_fasta: $methyl_control_fasta
    genome_path: $genome_path
    phix_path: $phix_path
    " >> $parameter_file 

##################################################################################################################################
##############################################---STEP 4: BCL TO FASTQ---######################################################### 
##################################################################################################################################
    fastq=$(find "$data_path/fastq" -type f | grep ".*\.fastq.gz$" | sort -u | wc -l)

    if [ $fastq -lt 1 ]; then
        echo "converting bcls to fastqs"
        cd $data_path
        module load bcl2fastq2
        bcl2fastq -o ./fastq -p 8
        cd $code_directory
        echo "conversion of bcls to fastqs complete"
    else
        echo "bcls already transformed into fastqs"
    fi

##################################################################################################################################
##############################################---STEP 5: RUN methylseq.sh---###################################################### 
##################################################################################################################################    
    fastq_file="${data_path}/fastq/FASTQs" #give a path to a file to store the paths to the fastq files in $fastq_directory
    echo "location of fastq_file: $fastq_file"
    find "$data_path/fastq" -type f | grep ".*\.fastq.gz$" | grep -v ".*\.trimmed.fastq.gz$" | sed -e 's/_R1.*$//g' | sed -e 's/_R2.*$//g' | sort -u > "${fastq_file}" #generate list of full paths to fastq files and save to the file in $fastq_list
    array_length=$(wc -l < "${fastq_file}") #get the number of files
    echo "array length: $array_length"
    echo "array length: $array_length
    " >> $parameter_file
   
    picard=$(find "$data_path/fastq/$unmethyl_control/$hydroxymethyl_control/$methyl_control/genome_alignment" -type f | grep ".*\.bam_picard_insert_size_plot.pdf$" | sort -u | wc -l)
    echo "picard: $picard"

    if [[ $picard -lt 1 ]] || [[ $force = true ]]; then
            echo "methylseq.sh running"
            sbatch -o "${data_path}/Logs/%A_%a.log" `#put into log` \
                    -a "1-${array_length}" `#initiate job array equal to the number of fastq files` \
                   -W `#indicates to the script not to move on until the sbatch operation is complete` \
                    "${code_directory}/methylseq.sh" \
                    $data_path $unmethyl_control_fasta $unmethyl_control $hydroxymethyl_control_fasta $hydroxymethyl_control $methyl_control_fasta $methyl_control $genome_path $phix_path $cores
            wait
            echo "methylseq.sh complete"
        else
            echo "picard files already created, methylseq.sh skipped"
    fi

##################################################################################################################################
#######################################---STEP 6: RUN extract_methylation.sh---################################################### 
##################################################################################################################################
    bam_file="$data_path/fastq/$unmethyl_control/$hydroxymethyl_control/$methyl_control/genome_alignment/bam_files" 
    find "$data_path/fastq/$unmethyl_control/$hydroxymethyl_control/$methyl_control/genome_alignment/split_bams" -type f | grep ".*\.bam$" | sort -u > "${bam_file}" #generate list of full paths to fastq files and save to the file in $fastq_list
    number_bam=$(wc -l < "${bam_file}") #get the number of files
    echo "number_bam: $number_bam"

    bedgraph=$(find "$data_path/fastq/$unmethyl_control/$hydroxymethyl_control/$methyl_control/genome_alignment/split_bams" -type f | grep ".*\.bedGraph.gz$" | sort -u | wc -l)
    echo "bedgraph: $bedgraph"

    if [[ $bedgraph -lt 2 ]] || [[ $force = true ]]; then
            echo "extract_methylation.sh running"
            sbatch -o "${data_path}/Logs/%A_%a.log" `#put into log` \
                    -a "1-${number_bam}" `#initiate job array equal to the number of fastq files` \
                    -W `#indicates to the script not to move on until the sbatch operation is complete` \
                    "${code_directory}/extract_methylation.sh" \
                    $data_path $unmethyl_control_fasta $unmethyl_control $hydroxymethyl_control_fasta $hydroxymethyl_control $methyl_control_fasta $methyl_control $genome_path $phix_path $cores
                wait
            echo "extract_methylation.sh complete"
            else
                echo "split coverage files and bed graphs have already been created for the genome alignemnt, extract_methylation.sh skipped"
    fi

##################################################################################################################################
##########################################---STEP 7: RUN join_coverage.sh---###################################################### 
##################################################################################################################################
    bedgraph=$(find "$data_path/fastq/$unmethyl_control/$hydroxymethyl_control/$methyl_control/genome_alignment" -maxdepth 1 -type f  | grep ".*\.bedGraph.gz$" | sort -u | wc -l)
    echo "bedgraph: $bedgraph"

    if [[ $bedgraph -lt 1 ]] || [[ $force = true ]] ; then
            echo "combining bed graphs and .cov files currently split by chromosome via join_coverage.sh"
            sbatch -o "${data_path}/Logs/%A_%a.log" `#put into log` \
                    -a "1-${array_length}" `#initiate job array equal to the number of fastq files` \
                    -W `#indicates to the script not to move on until the sbatch operation is complete` \
                    "${code_directory}/join_coverage.sh" \
                    $data_path  $unmethyl_control  $hydroxymethyl_control $methyl_control  $cores $temp_path
                wait
            echo "join_coverage.sh complete"
        else
            echo "combined bed graphs and coverage files have already been created for the genome alignemnt, join_coverage.sh skipped"
    fi

fi
