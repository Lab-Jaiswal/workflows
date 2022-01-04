#!/bin/bash
if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ]; then
    echo "Format: ./submit_BWA_CHIP.sh [fastq_directory] [output_directory] --argument"
    echo "fastq_directory: path to raw .fastq or .fastq.gz files"
    echo "output_directory: path for BWA and mutect output"
    echo "argument: indicates requested analysis technique(s) (--mutect, --varscan, --haplotypecaller, or --all)."
    echo "If --varscan is selected, you may use the optional arguments --p_value, --min_var_freq, and --min_coverage (if so desired)." 
    echo "If --mutect is selected, you may use --twist to indicate you would like your results filtered by the Twist panel"
    exit 1
else
    TEMP=`getopt -o vdm: --long min_coverage:,min_var_freq:,p_value:,twist,mutect,varscan,haplotypecaller,all \
    -n 'submit_BWA_CHIP.sh' -- "$@"`

   if [ $? != 0 ]; then
       echo "Unrecognized argument. Possible arguments: mutect, varscan, haplotypecaller, all, min_coverage, min_var_freq, p_value, and twist." >&2 ; exit 1 ; 
   fi
       eval set -- "$TEMP"

        get_mutect=false
        get_varscan=false
        get_haplotype=false
        all=false
        twist="0"
        min_coverage="10"
        min_var_freq="0.001"
        p_value="0.1"
        
    while true; do
        case "$1" in
            --min_coverage ) min_coverage="$2"; shift 2 ;;
            --min_var_freq ) min_var_freq="$2"; shift 2 ;;
            --p_value ) p_value="$2"; shift 2;;
            --twist ) twist="1"; shift ;;
            -m | --mutect ) get_mutect=true; shift ;;
            -v | --varscan ) get_varscan=true; shift ;;
            -h | --haplotypecaller ) get_haplotype=true; shift ;;
            --all ) get_mutect=true; get_varscan=true; get_haplotype=true; shift ;;
            -- ) shift; break ;;
            * ) break ;;
        esac
    done

    if ( [[ $min_coverage != "10" ]] || [[ $min_var_freq != "0.001" ]] || [[ $p_value != "0.1" ]] ) && \
        ( [[ $varscan = false ]] || [[ $all = false ]] ); then
        echo "The p_value, min_coverage, and min_var_freq arguments are not used in the mutect and haplotypecaller workflows"; exit 1
    fi

    if ( [[ twist = true ]] ) && \
        ( [[ $mutect = false ]] || [[ $all = false ]] ); then
        echo "The --twist argument is not applicable to the varscan and haplotypecaller workflows"; exit 1
    fi


    data_directory=$1 #get directory path from second argument (first argument $0 is the path of this script)
    #find "${data_directory}/" -type f | grep "bam" | grep -v ".bam.bai" | sed -e 's/\.bam$//g'
    output_directory=$2
    parent_directory=$(dirname $data_directory) #get parent directory of $fastq_directory
    code_directory=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P ) #specify location of star_align_and_qc.sh
    fastq_list="${parent_directory}/fastq_files" #give a path to a file to store the paths to the fastq files in $fastq_directory
    bam_list="${parent_directory}/bam_files" #give a path to a file to store the paths to the fastq files in $fastq_directory


    mkdir -p $output_directory
    mkdir -p "${output_directory}/Logs"


    find "${data_directory}/" -type f `#list all files in ${fastq_directory}` | \
        grep ".*\.bam$" `#only keep files with FASTQ in name (case insensitive)` | \
        grep -v ".bam.bai" `#remove Undetermined FASTQs` | \
        sed -e 's/\.bam$//g' `#remove _R1/2_fastq.gz file extension` | \
        sort -u  `#sort and remove duplicate names` > ${bam_list}
            #| \
        #head -n -1 > ${bam_list} `#remove the last line and generate a list of unique FASTQs`
    bam_array_length=$(wc -l < ${bam_list}) #get the number of FASTQs 
    echo "$bam_array_length"


    find "${data_directory}/" -type f `#list all files in ${fastq_directory}` | \
        grep ".*\.fastq.gz$" `#only keep files with FASTQ in name (case insensitive)` | \
        grep -v "Undetermined" `#remove Undetermined FASTQs` | \
        sed -e 's/_R1.*$//g' | sed -e 's/_R2.*$//g' `#remove _R1/2_fastq.gz file extension`| \
        sort -u  `#sort and remove duplicate names` > ${fastq_list} 
            #| \
        #head -n -1 > ${fastq_list} `#remove the last line and generate a list of unique FASTQs`
    fastq_array_length=$(wc -l < ${fastq_list}) #get the number of FASTQs 
    echo "$fastq_array_length"

    if [ $bam_array_length -ge 1 ] && [ $fastq_array_length -eq 0 ]; then
        sbatch -o "${output_directory}/Logs/%A_%a.log" `#put into log` \
            -a "1-${bam_array_length}" `#initiate job array equal to the number of fastq files` \
            -W `#indicates to the script not to move on until the sbatch operation is complete` \
            "${code_directory}/bam_to_fastq.sh" \
            "$data_directory" \
            "$code_directory"

        cp $bam_list $fastq_list
        fastq_array_length=$(wc -l < ${fastq_list}) #get the number of FASTQs 
        echo "$fastq_array_length"
            wait
        
    fi


    sbatch -o "${output_directory}/Logs/%A_%a.log" `#put into log` \
        -a "1-${fastq_array_length}" `#initiate job array equal to the number of fastq files` \
        -W `#indicates to the script not to move on until the sbatch operation is complete` \
        "${code_directory}/BWA_CHIP.sh" \
        ${parent_directory} \
        ${output_directory} \
        ${min_coverage} \
        ${min_var_freq} \
        ${p_value} \
        ${get_mutect} \
        ${get_varscan} \
        ${get_haplotype}
fi

wait

module load R/4.0

if [ $get_mutect = true ]; then
    Rscript aggregate_variants_mutect.R /labs/sjaiswal/chip_submitted_targets_Twist.xls \
        "$output_directory" "$twist" > "$output_directory/Logs/mutectOutFile.Rout" 2>&1
        
    #Rscript WhiteList/whitelist_mutect_join.R "/labs/sjaiswal/variant_whitelist.xlsx" \
      #  "$output_directory/mutect_aggregated_simple.tsv" "$output_directory" > "$output_directory/Logs/annotationOutFile.Rout" 2>&1

fi

if [ $get_varscan = true ]; then
    Rscript aggregate_variants_varscan.R /labs/sjaiswal/chip_submitted_targets_Twist.xls \
        "$output_directory" > "$output_directory/Logs/mutectOutFile.Rout" 2>&1
    Rscript WhiteList/whitelist_varscan_join.R "/labs/sjaiswal/variant_whitelist.xlsx" \
        "$output_directory/varscan_aggregated.tsv" "$output_directory" > "$output_directory/Logs/annotationOutFile.Rout" 2>&1
fi 

if [ $get_haplotype = true ]; then
    Rscript aggregate_variants_haplotypecaller.R /labs/sjaiswal/chip_submitted_targets_Twist.xls \
        "$output_directory" > "$output_directory/Logs/haplotypeOutFile.Rout" 2>&1
fi

mutectRout="$output_directory/Logs/mutectOutFile.Rout"
if grep -q "Can't combine" "$mutectRout"; then
    echo "There is an issue with the list containing maf dataframes. The column type varies from data frame to dataframe within list_of_mafs. Check mutectOutFile.Rout to determine which column(s) are causing this error. All columns sharing the same name must share the same type in order for the bind_rows function line 108 to work."
fi
