#!/bin/bash
##################################################################################################################################
#############################################---STEP 1: SET UP PARAMETERS---###################################################### 
##################################################################################################################################
if [ -z $1 ] || [ -z $2 ] || [ -z $3 ]; then
    echo "Format: ./submit_methylseq.sh [data_directory] [output_path] [genetic_locations_file]"
    echo "genetic_locations_file is a txt file"
    echo "user can use argument --cores to select a number of cores different from the default (24)"
    echo "user can use argument --log_name to specify a name for the log files"
    echo "user can use argument --force to force the code to disregard the checks preventing methylseq.sh, extract_methylation.sh & join_coverage.sh from running"
    exit 1
else

    TEMP=`getopt -o vdm: --long cores:,log_name:,force  -n 'submit_methylseq.sh' -- "$@"` #create optional arguments --cores, --force and --log_name
        eval set -- "$TEMP"
        
        cores=0                                                                           #if --cores # is not called, the default number is 0
        log_name="log"                                                                    #if --log_name is not called, log files will being with "log"
        force=false                                                                       #if --force is not called, force is equal to false
                        
    while true; do
        case "$1" in
            --cores ) cores="$2"; shift 2;;
            --log_name ) log_name="$2"; shift 2;;
            -f | --force ) force=true; shift ;;
            -- ) shift; break ;;
            * ) break ;;
        esac
    done                                                                 
    
    data_path=$1                                                                             
    output_path=$2
    genetic_locations=$3
    unmethyl_control=$(sed -n '1p' $genetic_locations)                                    #Collect the folder names and genome locations from the provided .txt file 
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
    if [ log == "log_" ]; then                                                             #give a path to a file to store the parameter files (so they are unique)
        parameter_file="${output_path}/Parameters/${now}_parameters.txt"                   #add date stamp to parameter files and, if provided, the log name
    else
        parameter_file="${log_name}${now}_parameters.txt"
    fi
        
    if [ $force == true ]; then                                                           #Create a variable called set_force which is empty if the 
        set_force="--force"                                                                     #user did not selecte --force and "--force" if they did
    fi
    if [ $cores -ne 0 ]; then                                                             #Create a variable called set_cores which is empty if the
        set_cores="--cores ${cores}"                                                             #user did not selecte --cores and "--cores #" if they did
    fi
    if [ $log_name != "log_" ]; then                                                      #Create a variable called set_log which is empty if the
        set_log="--log_name ${log_name}"                                                         #user did not selecte --log_name and "--log_name value" if they did
    fi
                                                                                          #We have to create these variable because these argument are NOT recognized
                                                                                                 #using traditional methods ($4, $5, etc.)     
        
    echo "call made to execute code: $0 $1 $2 $3 $set_force $set_cores $set_long
    " > $parameter_file
    
    if [ $force == true ] || [ $cores -ne 0 ] || [ $log_name != "log_" ]; then
        echo "you selected the following optional arguments: ${set_force}, ${set_cores}, ${set_log} 
            force is now equal to $force
            cores is now equal to $cores
            and your log files will begin with $log_name
        " >> $parameter_file
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
    fastq_file="${data_path}/fastq/FASTQs"                                              #give a path to a file to store the fastq file paths in $fastq_directory
    echo "location of fastq_file: $fastq_file"
    find "$data_path/fastq" -type f | grep ".*\.fastq.gz$" | \                          #generate list of full paths to fastq files and save to the file in $fastq_list
           grep -v ".*\.trimmed.fastq.gz$" | sed -e 's/_R1.*$//g' | \
           sed -e 's/_R2.*$//g' | sort -u > "${fastq_file}" 
    array_length=$(wc -l < "${fastq_file}")                                             #get the number of files and, thus, array length
    echo "array length: $array_length
    " >> $parameter_file
   
    picard=$(find "$output_path/$unmethyl_control/$hydroxymethyl_control/$methyl_control/genome_alignment" -type f | \
            grep ".*\.bam_picard_insert_size_plot.pdf$" | sort -u | wc -l)
    
    if [[ $picard -lt 1 ]] || [[ $force = true ]]; then
            echo "methylseq.sh running"
            sbatch -o "${output_path}/Logs/${log_name}_%A_%a.log" `#put into log` \
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
