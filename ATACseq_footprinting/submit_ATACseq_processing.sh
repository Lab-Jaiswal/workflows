#!/bin/bash

#SBATCH --time=72:00:00
#SBATCH --account=sjaiswal
#SBATCH --cpus-per-task=8
#SBATCH --mem=256GB
#SBATCH --job-name=submit_methylseq

##################################################################################################################################
#############################################---STEP 1: SET UP PARAMETERS---###################################################### 
##################################################################################################################################
if [ -z $1 ] || [ -z $2 ] || [ -z $3 ]; then
    echo "Format: ./submit_ATACseq_processing.sh [data_directory] [genome folder] [output path]"
    echo "user can use argument --gsize # to set gsize to a value other than 2620345972"
    echo "user can use argument --extsize # to set extsize to a vlaue other than 200"
    echo "user can use argument --shifts # to set shifts to a value other than 100"
    echo "user can use argument --broad false to change nomodel to false"
    echo "user can use argument --nomodel false to change nomodel to false"
    echo "user can use argument --blacklist to set a black list other than /oak/stanford/groups/sjaiswal/kameronr/ATACseq/blacklist/mm9-blacklist.bed.gz"
    echo "user can use argument --force to force the code to disregard the checks preventing methylseq.sh, extract_methylation.sh & join_coverage.sh from running"
    echo "user can use argument --log_name to give the log file a unique name"
    exit 1
else
#pipeline in shell for ATACseq footprinting
    TEMP=`getopt -o vdm: --long gsize:,extsize:,shifts:,broad:,nomodel:,blacklist:,log_name \
        -n './submit_atacseq' -- "$@"`

       if [ $? != 0 ]; then
           echo "Unrecognized argument. Possible arguments: --gize, --extsize, shifts, --broad true, --broad false, --nomodel true, --nomodel false, --blacklist path_to_blacklist --log_name chosen_log_name." >&2 ; exit 1 ; 
       fi
           eval set -- "$TEMP"

            gsize=2620345972
            extsize=200
            shifts=100
            broad=true
            nomodel=true
            blacklist="/oak/stanford/groups/sjaiswal/kameronr/ATACseq/blacklist/mm9-blacklist.bed.gz"
            log_name="log_"
            
        while true; do
            case "$1" in
                --gsize ) gsize="$2"; shift 2 ;;
                --extsize ) extsize="$2"; shift 2 ;;
                --shifts ) shifts="$2"; shift 2 ;;
                --broad ) broad="$2"; shift 2 ;;
                --nomodel ) nomodel="$2"; shift 2 ;;
                --blacklist ) blacklist="$2"; shift 2 ;;
                --log_name ) log_name="$2"; shift 2;;
                -- ) shift; break ;;
                * ) break ;;
            esac
        done

        echo "macs2 parameters: gsize: $gsize, extsize: $extsize, shifts: $shifts, broad: $broad, nomodel: $nomodel"

        if ( [[ $broad != true ]] && [[ $broad != false ]] ); then
            echo "broad can only be set to true or false"
            echo "example: --broad true"
            echo "example: --broad false"
            exit 1
        fi

        if ( [[ $nomodel != true ]] && [[ $nomodel != false ]] ); then
            echo "nomodel can only be set to true or false"
            echo "example: --nomodel true"
            echo "example: --nomodel false"
            exit 1
        fi
        
    bam_path=$1
    genome_path=$2
    genome_folder="$(dirname "${genome_path}")"
    output_path=$3
    echo "$genome_folder"

    code_directory=$( realpath . )

    module load samtools/1.9

##################################################################################################################################
#######################################---STEP 2: CREATE NECESSARY FOLDERS---#####################################################
##################################################################################################################################
    if [ ! -d "$output_path/Logs" ]; then 
            mkdir -p "$output_path/Logs"
    fi
            
    if [ ! -d "$output_path/Parameters" ]; then
            mkdir -p "$output_path/Parameters"
    fi

    if [ ! -d "$output_path/peak_calling" ]; then
            mkdir -p "$output_path/peak_calling"
    fi

    if [ ! -d "$output_path/coverage" ]; then
            mkdir -p "$output_path/coverage"
    fi

    Logs="${output_path}/Logs"
    Parameters="${initial_path}/Parameters"

##################################################################################################################################
#######################################---STEP 3: CREATE PARAMETER LOG---#########################################################
##################################################################################################################################
    now=$(date +%m_%d_%H_%M)
    if [ $log_name == "log_" ]; then                                                             #give a path to a file to store the parameter files (so they are unique)
        parameter_file="$Parameters/${now}_parameters.txt"                                 #add date stamp to parameter files and, if provided, the log name
    else
        parameter_file="$Parameters/${log_name}${now}_parameters.txt"
    fi
    
    touch $parameter_file
        
    if [ $gzise -ne 2620345972 ]; then                                                            
        set_gsize="--gsize $gsize"                                                                     
    fi
    if [ $extsize -ne 200 ]; then                                                             
        set_extsize="--extsize ${extsize}"                                                            
    fi
    if [ $shifts -ne 100 ]; then                                                      
        set_shifts="--shifts ${shifts}"                                                         
    fi
    if [ $broad != true ]; then                                                            
        set_broad="--broad $broad"                                                                     
    fi
    if [ $nomodel != true ]; then                                                             
        set_nomodel="--nomodel ${nomodel}"                                                            
    fi
    if [ $blacklist != "/oak/stanford/groups/sjaiswal/kameronr/ATACseq/blacklist/mm9-blacklist.bed.gz" ]; then                                                             
        set_blacklist="--blacklist ${blacklist}"                                                            
    fi
    if [ $log_name != "log_" ]; then                                                      
        set_log="--log_name ${log_name}"                                                         
    fi

    echo "location of scripts used to run code : $code_directory
        " > $parameter_file 
        
    echo "call made to execute code: $0 $1 $2 $3 $set_gsize $set_extsize $set_shifts $set_broad $set_nomodel $set_blacklist $set_log
    " >> $parameter_file
    
    if [ $gzise -ne 2620345972 ] || [ $extsize -ne 200 ] || [ $shifts -ne 100 ] || [ $broad != true ] || [ $nomodel != true ] || [ $blacklist != "/oak/stanford/groups/sjaiswal/kameronr/ATACseq/blacklist/mm9-blacklist.bed.gz" ] || [ $log_name != "log_" ]; then
        echo "you selected the following optional arguments: $set_gsize, $set_extsize, $set_shifts, $set_broad, $set_nomodel, $set_blacklist, $set_log 
            gsize is now equal to $gsize
            extsize is now equal to $extsize
            shifts is now equal to $shifts
            broad is now equal to $broad
            nomodel is now equal to $nomodel
            blacklist is now equal to $blacklist
            and your log files will begin with $log_name
        " >> $parameter_file
    fi
        
##################################################################################################################################
########################################---STEP 1: CREATE COVERAGE BED FILES---################################################### 
##################################################################################################################################
    if [ ! -f $genome_folder/chromsizes.bed ]; then
        cd $genome_folder
        samtools faidx mm9_bgzip.fa.gz
        cut -f1,2 mm9_bgzip.fa.gz.fai > sizes.genome

        cat chromsizes.txt | awk '{{ print $1, 0, $2 }}' | tr ' ' '\t' > chromsizes.bed

        echo "creation of chromosome-based coverage bed files complete"
    else
        echo "chromosome-based coverage bed files have already been created"
    fi

    whitelist=$genome_folder/chromsizes.bed

##################################################################################################################################
#####################################---STEP 2: SORT, MERGE, AND INDEX BAM FILES---############################################### 
##########################################---CREATE COVERAGE BIGWIG TRACK---######################################################
############################################---PEAK CALLING WITH MACS2---#########################################################
###########################################---REMOVE BLACKLISTED REGIONS---#######################################################
##################################################################################################################################
    cd $bam_path 
    bam_file="$bam_path/BAMs" #give a path to a file to store the paths to the bams files in $bam_directory

    find "${bam_path}/" -type f `#list all files in ${fastq_directory}` | \
            grep ".*\.bam$" `#only keep files with FASTQ in name (case insensitive)` | \
            grep -v ".sorted.bam" `#remove sorted BAMs` | \
            grep -v ".merged.bam" `#remove merged BAMs` | \
            sed -e 's/\_Rep.*//g' `#remove _Rep1/2_treat_rep1.bam extension` | \
            sort -u  `#sort and remove duplicate names` > ${bam_file}

        bam_array_length=$(wc -l < ${bam_file}) #get the number of FASTQs 
        echo "$bam_array_length"

    number_bams=$(wc -l < "${bam_file}") #get the number of files
    array_length="$number_bams"

    bed_file=$(find "$output_path/peak_calling" -maxdepth 2 -type f | grep "raw.bed" | sort -u | wc -l)

    echo "entering the script" >> $parameter_file

    if [ $bed_file -lt 1 ]; then
             sbatch -o "$Logs/%A_%a.log" `#put into log` \
            -a "1-${array_length}" `#initiate job array equal to the number of bam files` \
            -W `#indicates to the script not to move on until the sbatch operation is complete` \
                "${code_directory}/ATACseq_processing.sh" \
                $bam_path $output_path $gsize $extsize $shifts $broad $nomodel $blacklist $whitelist $genome_folder $parameter_file
            wait
        else
            echo "sorting, merging, and indexing of files already completed"
    fi

##################################################################################################################################
#########################---STEP 4: PEAK PROCESSING: REDUCE GENOMIC LOCATION COLUMNS AND SORT---################################## 
##################################################################################################################################
    if [ ! -f "$output_path/peak_calling/all_merged.bed" ]; then
        temp_path=$(mktemp -d /tmp/tmp.XXXXXXXXXX)
        echo "temp_path is: " $temp_path
        mkdir $temp_path

        cat $bam_file | while read line; do 
            PREFIX=$(basename $line)
            chmod 775 $output_path/peak_calling/$PREFIX/${PREFIX}_union.bed
            cp $output_path/peak_calling/$PREFIX/${PREFIX}_union.bed $temp_path
        done

        cd $temp_path
        touch $temp_path/merged.txt
        cat *.bed > $temp_path/merged.txt
        cat merged.txt | tr ' ' '\t' | sort -k1,1 -k2,2n | bedtools merge -d 5 -c 4 -o distinct > all_merged.bed
        cp $temp_path/all_merged.bed $output_path/peak_calling/all_merged.bed
        echo "merged union file complete"
    else
        echo "creation of merged union file already finished"
    fi

fi
