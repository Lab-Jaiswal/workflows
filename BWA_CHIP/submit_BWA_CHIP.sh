#!/bin/bash
if [ -z "$1" ] || [ -z "$2" ]; then
    echo "run_BWA_mutect [fastq_directory] [output_directory]"
    echo "fastq_directory: path to raw .fastq or .fastq.gz files"
    echo "output_directory: path for BWA and mutect output"
else
    TEMP=`getopt -o vdm: --long mutect,varscan,haplotypecaller,all,min_coverage:,min_var_freq:,p_value: \
    -n 'javawrap' -- "$@"`

    if [ $? != 0 ] ; then echo "Unrecognized argument. Possible arguments: mutect, varscan, haplotypecaller, all, min_coverage, min_var_freq, and p_value." >&2 ; exit 1 ; fi
        eval set -- "$TEMP"

        get_mutect=false
        get_varscan=false
        get_haplotype=false
        all=false
        min_coverage="10"
        min_var_freq="0.001"
        p_value="0.1"

    while true; do
        case "$1" in
        -m | --mutect ) get_mutect=true; shift ;;
        -v | --varscan ) get_varscan=true; shift ;;
        -h | --haplotypecaller ) get_haplotype=true; shift ;;
        --all ) get_mutect=true; get_varscan=true; get_haplotype=true; shift ;;
        --min_coverage ) min_coverage="$2"; shift 2 ;;
        --min_var_freq ) min_var_freq="$2"; shift 2 ;;
        --p_value ) p_value="$2"; shift 2;;
        -- ) shift; break ;;
        * ) break ;;
        esac
    done

    fastq_directory=$1 #get directory path from second argument (first argument $0 is the path of this script)
    output_directory=$2
    parent_directory=$(dirname $fastq_directory) #get parent directory of $fastq_directory
    code_directory=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P ) #specify location of star_align_and_qc.sh
    fastq_list="${parent_directory}/fastq_files" #give a path to a file to store the paths to the fastq files in $fastq_directory

    mkdir -p $output_directory
    mkdir -p "${output_directory}/Logs"

    find "${fastq_directory}/" -type f `#list all files in ${fastq_directory}` | \
    grep "fastq" `#only keep files with FASTQ in name (case insensitive)` | \
    grep -v "Undetermined" `#remove Undetermined FASTQs` | \
    sed -e 's/_R1.*$//g' | sed -e 's/_R2.*$//g' `#remove _R1/2_fastq.gz file extension`| \
    sort -u > ${fastq_list} #sort and remove duplicate names, to generate list of unique FASTQs
    array_length=$(wc -l < ${fastq_list}) #get the number of FASTQs 

    sbatch -o "${output_directory}/Logs/%A_%a.log" `#put into log` \
            -a "1-${array_length}" `#initiate job array equal to the number of fastq files` \
            -W `#indicates to the script not to move on until the sbatch operation is complete` \
            "${code_directory}/BWA_CHIP.sh" \
            ${parent_directory} \
            ${output_directory} \
            ${min_coverage} \
            ${min_var_freq} \
            ${p_value} \
            ${get_mutect}\
            ${get_varscan}\
            ${get_haplotype}
fi

wait

if [ $get_mutect = true ]; then
    Rscript aggregate_variants_mutect.R /labs/sjaiswal/chip_submitted_targets_Twist.xls "$output_directory" > "$output_directory/Logs/mutectOutFile.Rout" 2>&1
else 
    echo "no mutect analysis requested"
fi 

if [ $get_varscan = true ]; then
    Rscript aggregate_variants_varscan.R /labs/sjaiswal/chip_submitted_targets_Twist.xls "$output_directory" > "$output_directory/Logs/mutectOutFile.Rout" 2>&1
else 
    echo "no varscan analysis requested"
fi 

if [ $get_haplotype = true ]; then
    Rscript aggregate_variants_haplotypecaller.R /labs/sjaiswal/chip_submitted_targets_Twist.xls "$output_directory" > "$output_directory/Logs/haplotypeOutFile.Rout" 2>&1
else
    echo "no haplotypecaller analysis requested"
fi
