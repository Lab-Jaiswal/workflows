#!/bin/bash
##################################################################################################################################
#############################################---STEP 1: SET UP PARAMETERS---###################################################### 
##################################################################################################################################
if [ -z $1 ] || [ -z $2 ] || [ -z $3 ]; then
    echo "Format: ./submit_methylseq.sh [data_directory] [output_path] [genetic_locations_file]"
    echo "genetic_locations_file is a txt file"
    echo "user can use argument --force to force the code to disregard the checks preventing methylseq.sh, extract_methylation.sh & join_coverage.sh from running"
    echo "user can use argument --cores to select a number of cores different from the default (24)"
    echo "user can use argument --log_name to specify a name for the log files"
    exit 1
else

    TEMP=`getopt -o vdm: --long cores:,log_name:,force  -n 'submit_methylseq.sh' -- "$@"`
        eval set -- "$TEMP"
        
        force=false
        cores=24
        log_name="log_"
                        
    while true; do
        case "$1" in
            -f | --force ) force=true; shift 2 ;;
            --cores ) cores="$2"; shift 2;;
            --log_name ) log_name="$2"; shift 2;;
            -- ) shift; break ;;
            * ) break ;;
        esac
    done

    if [ $log_name != "log_" ]; then
        log_name=$(echo "${log_name}_")
    fi

    echo "force: $force"
    data_path=$1
    output_path=$2
    genetic_locations=$3
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
   
    if [ ! -d $output_path ]; then
        mkdir $output_path
    fi

    if [ ! -d "$output_path/$unmethyl_control" ];then
        mkdir "$output_path/$unmethyl_control"
    fi

    if [ ! -d "$output_path/$unmethyl_control/$hydroxymethyl_control" ]; then
        mkdir "$output_path/$unmethyl_control/$hydroxymethyl_control"
    fi

    if [ ! -d "$output_path/$unmethyl_control/$hydroxymethyl_control/$methyl_control" ]; then
        mkdir "$output_path/$unmethyl_control/$hydroxymethyl_control/$methyl_control"
    fi

    if [ ! -d "$output_path/$unmethyl_control/$hydroxymethyl_control/$methyl_control/genome_alignment" ]; then
        mkdir "$output_path/$unmethyl_control/$hydroxymethyl_control/$methyl_control/genome_alignment"
    fi

    if [ ! -d "$output_path/Logs" ]; then
        mkdir "$output_path/Logs"
    fi

    if [ ! -d "$output_path/Parameters" ]; then
        mkdir "$output_path/Parameters"
    fi

##################################################################################################################################
#######################################---STEP 3: CREATE PARAMETER LOG---#####################################################
##################################################################################################################################
    now=$(date +%m_%d_%H_%M)
    parameter_file="${output_path}/Parameters/${now}_parameters.txt" #give a path to a file to store the paths to the fastq files in $fastq_directory
    
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
   
    picard=$(find "$output_path/$unmethyl_control/$hydroxymethyl_control/$methyl_control/genome_alignment" -type f | grep ".*\.bam_picard_insert_size_plot.pdf$" | sort -u | wc -l)
    echo "picard: $picard"

    if [[ $picard -lt 1 ]] || [[ $force = true ]]; then
            echo "methylseq.sh running"
            sbatch -o "${output_path}/Logs/${log_name}%A_%a.log" `#put into log` \
                    -a "1-${array_length}" `#initiate job array equal to the number of fastq files` \
                   -W `#indicates to the script not to move on until the sbatch operation is complete` \
                    "${code_directory}/methylseq.sh" \
                    $data_path $output_path $unmethyl_control_fasta $unmethyl_control $hydroxymethyl_control_fasta $hydroxymethyl_control $methyl_control_fasta $methyl_control $genome_path $phix_path $cores $log_name $parameter_file
            wait
            echo "methylseq.sh complete"
        else
            echo "picard files already created, methylseq.sh skipped"
    fi

#####################previously report_controls.sh################################
################################################################################
module load bismark
    if [ ! -f "$output_path/$unmethyl_control/bismark_summary_report.txt" ]; then 
        cd $output_path/$unmethyl_control
            echo "$output_path/$unmethyl_control/bismark_summary_report.txt does not exist yet"
        #http://felixkrueger.github.io/Bismark/Docs/
        #bismark report options:
        #--alignment_report FILE
        #--dedup_report FILE
        #--splitting_report FILE
        #--mbias_report FILE
        #--nucleotide_report FILE
        
        #https://rawgit.com/FelixKrueger/Bismark/master/Docs/Bismark_User_Guide.html
        bismark2report
        bismark2summary
        #need to specify a nucleotide coverage report file in the above command!
            echo "created report for unmethyl control"
    else 
            echo "bismark summary already completed for unmethyl control"
    fi
        
    if [ ! -f "$output_path/$unmethyl_control/$hydroxymethyl_control/bismark_summary_report.txt" ]; then
        cd $output_path/$unmethyl_control/$hydroxymethyl_control
            echo "$output_path/$unmethyl_control/$hydroxymethyl_control/bismark_summary_report.txt does not exist yet"
        bismark2report
        bismark2summary
        #need to specify a nucleotide coverage report file in the above command! ^
            echo "created report for hydroxymethyl control"
    else
        echo "bismark summary already completed for hydroxymethyl control"
    fi

    if [ ! -f "$output_path/$unmethyl_control/$hydroxymethyl_control/$methyl_control/bismark_summary_report.txt" ]; then
        cd $output_path/$unmethyl_control/$hydroxymethyl_control/$methyl_control
            echo "$output_path/$unmethyl_control/$hydroxymethyl_control/$methyl_control/bismark_summary_report.txt does not exist yet"
        bismark2report
        bismark2summary
        #need to specify a nucleotide coverage report file in the above command! ^
            echo "created report for methyl control"
    else
        echo "bismark summary already completed for methyl control"
    fi

    echo "report_controls complete for unmethyl, hydroxymethyl, and methyl control sequences"

    if [ ! -f "$output_path/$unmethyl_control/$hydroxymethyl_control/$methyl_control/genome_alignment/bismark_summary_report.txt" ]; then 
        #https://rawgit.com/FelixKrueger/Bismark/master/Docs/Bismark_User_Guide.html
        cd $output_path/$unmethyl_control/$hydroxymethyl_control/$methyl_control/genome_alignment
            echo "$output_path/$unmethyl_control/$hydroxymethyl_control/$methyl_control/genome_alignment/bismark_summary_report.txt does not exist yet"
        bismark2report
        bismark2summary
            echo "report complete"
    else
        echo "bismark summary found and already created for genome alignment"
    fi
fi
