#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --account=sjaiswal
#SBATCH --cpus-per-task=8
#SBATCH --mem=256GB
#SBATCH --job-name=submit_methylseq

##################################################################################################################################
#############################################---STEP 1: SET UP PARAMETERS---###################################################### 
##################################################################################################################################
if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ]; then
    echo "Format: ./submit_BWA_CHIP.sh [fastq_directory] [output_directory] --argument"
    echo "fastq_directory: path to raw .fastq or .fastq.gz files"
    echo "output_directory: path for BWA and mutect output"
    echo "argument: indicates requested analysis technique(s) (--mutect, --varscan, --haplotypecaller, or --all)."
    echo "If --varscan is selected, you may use the optional arguments --p_value, --min_var_freq, and --min_coverage (if so desired)." 
    echo "If --mutect is selected, you may use --twist to indicate you would like your results filtered by the Twist panel"
    exit 1
else
    TEMP=`getopt -o vdm: --long min_coverage:,min_var_freq:,p_value:,intervals:,paired:,log_name:,bam,twist,mutect,varscan,haplotypecaller,all \
    -n 'submit_BWA_CHIP.sh' -- "$@"`

   if [ $? != 0 ]; then
       echo "Unrecognized argument. Possible arguments: mutect, varscan, haplotypecaller, all, min_coverage, min_var_freq, p_value, and twist." >&2 ; exit 1 ; 
   fi
       eval set -- "$TEMP"

        get_mutect=false
        get_varscan=false
        get_haplotype=false
        all=false
        intervals=false
        paired=false
        normal_sample=false
        log_name="log_"
        use_bam=false
        twist="0"
        min_coverage="10"
        min_var_freq="0.001"
        p_value="0.1"
        
    while true; do
        case "$1" in
            --min_coverage ) min_coverage="$2"; shift 2 ;;
            --min_var_freq ) min_var_freq="$2"; shift 2 ;;
            --p_value ) p_value="$2"; shift 2 ;;
            --intervals ) intervals="$2"; shift 2 ;;
            -p | --paired ) paired=true; normal_sample="$2"; shift 2 ;; 
            --log_name ) log_name="$2"; shift 2;;
            --bam ) use_bam=true; shift ;;
            --twist ) twist="1"; shift ;;      
            -m | --mutect ) get_mutect=true; shift ;;
            -v | --varscan ) get_varscan=true; shift ;;
            -h | --haplotypecaller ) get_haplotype=true; shift ;;
            --all ) get_mutect=true; get_varscan=true; get_haplotype=true; shift ;;
            -- ) shift; break ;;
            * ) break ;;
        esac
    done

    echo "intervals: $intervals"
    if ( [[ $min_coverage != "10" ]] || [[ $min_var_freq != "0.001" ]] || [[ $p_value != "0.1" ]] ) && \
        ( [[ $get_varscan = false ]] || [[ $all = false ]] ); then
        echo "The p_value, min_coverage, and min_var_freq arguments are not used in the mutect and haplotypecaller workflows"; exit 1
    fi

    if ( [[ twist = true ]] ) && \
        ( [[ $get_mutect = false ]] || [[ $all = false ]] ); then
        echo "The --twist argument is not applicable to the varscan and haplotypecaller workflows"; exit 1
    fi


    data_directory=$1 #get directory path from second argument (first argument $0 is the path of this script)
    #find "${data_directory}/" -type f | grep "bam" | grep -v ".bam.bai" | sed -e 's/\.bam$//g'
    output_directory=$2
    parent_directory=$(dirname $data_directory) #get parent directory of $fastq_directory
    code_directory=$( realpath . ) #specify location of star_align_and_qc.sh
    fastq_list="${parent_directory}/fastq_files" #give a path to a file to store the paths to the fastq files in $fastq_directory
    bam_list="${parent_directory}/bam_files" #give a path to a file to store the paths to the fastq files in $fastq_directory

##################################################################################################################################
#######################################---STEP 2: CREATE NECESSARY FOLDERS---#####################################################
##################################################################################################################################
    if [ ! -d "$output_directory/Logs" ]; then
        
        mkdir -p "$output_directory/Logs"
    fi
        

    if [ ! -d "$output_directory/Parameters" ]; then
        mkdir -p "$output_directory/Parameters"
    fi
    
    Logs="${output_directory}/Logs"
    Parameters="${output_directory}/Parameters"

##################################################################################################################################
#######################################---STEP 3: CREATE PARAMETER LOG---#####################################################
##################################################################################################################################
    now=$(date +%m_%d_%H_%M)
    if [ $log_name == "log_" ]; then                                                             #give a path to a file to store the parameter files (so they are unique)
        parameter_file="$Parameters/${now}_parameters.txt"                                 #add date stamp to parameter files and, if provided, the log name
    else
        parameter_file="$Parameters/${log_name}${now}_parameters.txt"
    fi
    touch $parameter_file
        
    if [ $get_mutect != false ]; then                                           #Create a variable called set_log which is empty if the
        set_mutect="--mutect"                                                   #user did not selecte --log_name and "--log_name value" if they did
    fi                                                                          #We have to create these variable because these argument are NOT recognized
    if [ $get_varscan != false ]; then                                          #using traditional methods ($4, $5, etc.) 
        set_varscan="--varscan"
    fi
    if [ $get_haplotype != false ]; then
        set_haplotype="--haplotypecaller"
    fi
    if [ $all != false ]; then
        set_all="--mutect --varscan --haplotypecaller"
    fi
    if [ $intervals != false ]; then
        set_intervals="--intervals $intervals"
    fi
    if [ $paired != false ]; then
        set_paired="--paired $normal_sample"
    fi
    if [ $log_name != "log_" ]; then
        set_log="--log_name ${log_name}"
    fi
    if [ $use_bam != false ]; then
        set_bam="--bam"
    fi
    if [ $twist -ne 0 ]; then                                                      
        set_twist="--twist"                                                         
    fi
    if [ $min_coverage != "10" ]; then
        set_min_coverage="--min_coverage $min_coverage"
    fi
    if [ $min_var_freq != "0.001" ]; then
        set_min_var_freq="--min_var_freq $min_var_freq"
    fi
    if [ $p_value != "0.1" ]; then
        set_p_value="--p_value $p_value"
    fi
    
    
    echo "location of scripts used to run code : $code_directory
        " > $parameter_file 
        
    echo "call made to execute code: $0 $1 $2 $set_mutect $set_varscan $set_haplotype $set_all $set_intervals $set_paired $set_log $set_bam $set_twist $set_min_coverage $set_min_var_freq $set_p_value
    " >> $parameter_file
    
    echo "you chose the following optional arguments: 
        $set_mutect $set_varscan $set_haplotype $set_all $set_intervals $set_paired $set_log $set_bam $set_twist $set_min_coverage $set_min_var_freq $set_p_value
    " >> $parameter_file
        
 
##################################################################################################################################
#############################################--STEP 4: GET ARRAY LENGTHS---####################################################### 
##################################################################################################################################
    find -L "${data_directory}/" -type f `#list all files in ${fastq_directory}` | \
        grep ".*\.bam$" `#only keep files with FASTQ in name (case insensitive)` | \
        grep -v ".bai" `#remove Undetermined FASTQs` | \
        sed -e 's/\.bam$//g' `#remove _R1/2_fastq.gz file extension` | \
        sort -u  `#sort and remove duplicate names` > ${bam_list}
            #| \
        #head -n -1 > ${bam_list} `#remove the last line and generate a list of unique FASTQs`
    bam_array_length=$(wc -l < ${bam_list}) #get the number of FASTQs 
    echo "$bam_array_length"


    find -L "${data_directory}/" -type f `#list all files in ${fastq_directory}` | \
        grep ".*\.fastq.gz$" `#only keep files with FASTQ in name (case insensitive)` | \
        grep -v "Undetermined" `#remove Undetermined FASTQs` | \
        sed -e 's/_R1.*$//g' | sed -e 's/_R2.*$//g' `#remove _R1/2_fastq.gz file extension`| \
        sort -u  `#sort and remove duplicate names` > ${fastq_list} 
            #| \
        #head -n -1 > ${fastq_list} `#remove the last line and generate a list of unique FASTQs`
    fastq_array_length=$(wc -l < ${fastq_list}) #get the number of FASTQs 
    echo "$fastq_array_length"

##################################################################################################################################
################################################--STEP 5: BAM TO FASTQ---######################################################### 
##################################################################################################################################
    if [[ $bam_array_length -ge 1 ]] && [[ $use_bam == false ]]; then
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

        array_length=$fastq_array_lenth
    else
        array_length=$bam_array_length
    fi

##################################################################################################################################
##################################################--STEP 4: BWA_CHIP.sh---######################################################## 
##################################################################################################################################
    
    sbatch -o "${output_directory}/Logs/%A_%a.log" `#put into log` \
    -a "1-${array_length}%50" `#initiate job array equal to the number of fastq files` \
    -W `#indicates to the script not to move on until the sbatch operation is complete` \
    "${code_directory}/BWA_CHIP.sh" \
    ${parent_directory} \
    ${output_directory} \
    ${min_coverage} \
    ${min_var_freq} \
    ${p_value} \
    ${get_mutect} \
    ${get_varscan} \
    ${get_haplotype} \
    ${use_bam} \
    ${intervals} \
    ${paired} \
    ${normal_sample} \
    ${code_directory}

    wait
fi
    
##################################################################################################################################
################################################--STEP 4: RUN R SCRIPTS---######################################################## 
##################################################################################################################################
module load R/4.0

if [ $get_mutect = true ]; then
    Rscript aggregate_variants_mutect.R /labs/sjaiswal/chip_submitted_targets_Twist.xls \
        "$output_directory" "$twist" > "$output_directory/Logs/mutectOutFile.Rout" 2>&1
        
    Rscript WhiteList/whitelist_mutect_join.R "/labs/sjaiswal/variant_whitelist.xlsx" \
        "$output_directory/mutect_aggregated_simple.tsv" "$output_directory" > "$output_directory/Logs/annotationOutFile.Rout" 2>&1

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
