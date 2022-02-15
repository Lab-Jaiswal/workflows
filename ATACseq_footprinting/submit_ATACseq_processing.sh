#!/bin/bash

#pipeline in shell for ATACseq footprinting
TEMP=`getopt -o vdm: --long gsize:,extsize:,shifts:,broad:,nomodel:,blacklist \
    -n './submit_atacseq' -- "$@"`

   if [ $? != 0 ]; then
       echo "Unrecognized argument. Possible arguments: --gize, --extsize, shifts, --broad true, --broad false, --nomodel true, --nomodel false, --blacklist path_to_blacklist." >&2 ; exit 1 ; 
   fi
       eval set -- "$TEMP"

        gsize=2620345972
        extsize=200
        shifts=100
        broad=true
        nomodel=true
        blacklist="/oak/stanford/groups/sjaiswal/kameronr/ATACseq/blacklist/mm9-blacklist.bed.gz"
        
    while true; do
        case "$1" in
            --gsize ) gsize="$2"; shift 2 ;;
            --extsize ) extsize="$2"; shift 2 ;;
            --shifts ) shifts="$2"; shift 2 ;;
            --broad ) broad="$2"; shift 2 ;;
            --nomodel ) nomodel="$2"; shift 2 ;;
            --blacklist ) blacklist="$2"; shift 2 ;;
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
echo "$genome_folder"

code_directory=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )

#organism=$1
#fasta=$2
#blacklist=$3
#gtf=$4
#motifs=$5
#output=$6
#macs=$7
#mac="--nomodel --shift -100 --extsize 200 --broad"

#FASTA_PATH=/oak/stanford/groups/sjaiswal/kameronr/sjaiswal_old/genomes/mm9/mm9_bgzip.fa.gz
#BLACKLIST_PATH=/oak/stanford/groups/sjaiswal/kameronr/ATACseq/blacklist/mm9-blacklist.bed.gz
#GTF_PATH=/oak/stanford/groups/sjaiswal/kameronr/sjaiswal_old/genomes/mm9/GTF/gencode.vM1.annotation.gtf.gz
#MOTIFS_DIR=/oak/stanford/groups/sjaiswal/kameronr/ATACseq/TOBIAS_snakemake/data/individual_motifs
#OUTPUT_DIR=/oak/stanford/groups/smontgom/kameronr/ATACseq/output

module load samtools/1.9

##################################################################################################################################
########################################---STEP 1: CREATE COVERAGE BED FILES---################################################### 
##################################################################################################################################
#get chromosomes available in fasta (fasta chroms) -- because it's needed for making the bigwig tracks.
#OUTPUT_DIR=/oak/stanford/groups/smontgom/kameronr/ATACseq/test
#this next line might not work? seems like a python command? https://www.biostars.org/p/173963/

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
            #| \
        #head -n -1 > ${bam_list} `#remove the last line and generate a list of unique FASTQs`
    bam_array_length=$(wc -l < ${bam_file}) #get the number of FASTQs 
    echo "$bam_array_length"

number_bams=$(wc -l < "${bam_file}") #get the number of files
array_length="$number_bams"

if ! [ -d "$bam_path/Logs" ]; then
    mkdir -p "$bam_path/Logs"
fi

bed_file=$(find "$bam_path/peak_calling" -maxdepth 2 -type f | grep "raw.bed" | sort -u | wc -l)

if [ $bed_file -lt 1 ]; then
         sbatch -o "${bam_path}/Logs/%A_%a.log" `#put into log` \
        -a "1-${array_length}" `#initiate job array equal to the number of bam files` \
        -W `#indicates to the script not to move on until the sbatch operation is complete` \
            "${code_directory}/ATACseq_processing.sh" \
            $bam_path $gsize $extsize $shifts $broad $nomodel $blacklist $whitelist $genome_folder
        
        wait
    else
        echo "sorting, merging, and indexing of files already completed"
fi

##################################################################################################################################
#########################---STEP 4: PEAK PROCESSING: REDUCE GENOMIC LOCATION COLUMNS AND SORT---################################## 
##################################################################################################################################
if [ ! -f "$bam_path/peak_calling/all_merged.bed" ]; then
    temp_path=$(mktemp -d /tmp/tmp.XXXXXXXXXX)
    echo "temp_path is: " $temp_path
    mkdir $temp_path

    cat $bam_file | while read line; do 
        PREFIX=$(basename $line)
        chmod 775 $bam_path/peak_calling/$PREFIX/${PREFIX}_union.bed
        cp $bam_path/peak_calling/$PREFIX/${PREFIX}_union.bed $temp_path
    done

    cd $temp_path
    touch $temp_path/merged.txt
    cat *.bed > $temp_path/merged.txt
    cat merged.txt | tr ' ' '\t' | sort -k1,1 -k2,2n | bedtools merge -d 5 -c 4 -o distinct > all_merged.bed
    cp $temp_path/all_merged.bed $bam_path/peak_calling/all_merged.bed
    cd $bam_path
    echo "merged union file complete"
else
    echo "creation of merged union file already finished"
fi

#4d. peak annotation. peaks per condition or across conditions, dependent on run info output
#need to make .config file for uropa to use.
#could later add expression information to each peaks if needed?
