#!/bin/bash

#SBATCH --time=167:00:00
#SBATCH --partition=nih_s10
#SBATCH --account=sjaiswal
#SBATCH --cpus-per-task=1
#SBATCH --mem=4GB
#SBATCH --job-name=submit_methylseq

##################################################################################################################################
#############################################---STEP 1: SET UP PARAMETERS---###################################################### 
##################################################################################################################################
if [ -z $1 ] || [ -z $2 ] || [ -z $3 ]; then
    echo "Format: ./submit_methylseq.sh [data_directory] [output_path] [genetic_locations_file]"
    echo "genetic_locations_file is a txt file"
    echo "user can use argument --cores to select a number of cores different from the default (24)"
    echo "user can use argument --log_name to specify a name forthe log files"
    echo "user can use argument --force to force the code to disregard the checks preventing methylseq.sh, extract_methylation.sh & join_coverage.sh from running"
    exit 1
else

    TEMP=`getopt -o vdm: --long cores:,log_name:,force  -n 'submit_methylseq.sh' -- "$@"` #create optional arguments --cores, --force and --log_name
        eval set -- "$TEMP"
        
        cores=24                                                                           #if --cores # is not called, the default number is 0
        log_name="log_"                                                                    #if --log_name is not called, log files will being with "log"
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
    initial_path=$output_path
    genetic_locations=$3
    line_count=$( wc -l < "${genetic_locations}" )
    genome_count=$(bc -l <<< "scale=1; ($line_count / 3)")
    code_directory=$( realpath . )
    
##################################################################################################################################
#######################################---STEP 2: CREATE NECESSARY FOLDERS---#####################################################
##################################################################################################################################
    if [ ! -d "$initial_path/Logs" ]; then
        
        mkdir -p "$initial_path/Logs"
    fi
        

    if [ ! -d "$initial_path/Parameters" ]; then
        mkdir -p "$initial_path/Parameters"
    fi
    
    Logs="${initial_path}/Logs"
    Parameters="${initial_path}/Parameters"

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
        
    if [ $force == true ]; then                                                           #Create a variable called set_force which is empty if the 
        set_force="--force"                                                                     #user did not selecte --force and "--force" if they did
    fi
    if [ $cores -ne 24 ]; then                                                             #Create a variable called set_cores which is empty if the
        set_cores="--cores ${cores}"                                                             #user did not selecte --cores and "--cores #" if they did
    fi
    if [ $log_name != "log_" ]; then                                                      #Create a variable called set_log which is empty if the
        set_log="--log_name ${log_name}"                                                         #user did not selecte --log_name and "--log_name value" if they did
    fi
                                                                                          #We have to create these variable because these argument are NOT recognized
                                                                                                 #using traditional methods ($4, $5, etc.)     

    echo "location of scripts used to run code : $code_directory
        " > $parameter_file 
        
    echo "call made to execute code: $0 $1 $2 $3 $set_force $set_cores $set_long
    " >> $parameter_file
    
    if [ $force == true ] || [ $cores -ne 24 ] || [ $log_name != "log_" ]; then
        echo "you selected the following optional arguments: ${set_force}, ${set_cores}, ${set_log} 
            force is now equal to $force
            cores is now equal to $cores
            and your log files will begin with $log_name
        " >> $parameter_file
    fi
        
    echo "location of file with genome directions: $genetic_locations
    " >> $parameter_file

    total_genomes=$(bc -l <<< "scale=0; (($line_count / 3) - 1)")
    echo "number of control genomes provided: $total_genomes
    " >> $parameter_file

    for i in $(seq 0 $total_genomes); do
        number1=$(bc -l <<< "scale=0; (($i * 3) +1)")
        number2=$(bc -l <<< "scale=0; (($i * 3) +2)")
        number3=$(bc -l <<< "scale=0; (($i * 3) +3)") 
        genome_name=$(sed -n ${number1}'p' $genetic_locations)
        genome_fasta_path=$(sed -n ${number2}'p' $genetic_locations)
        deduplicate=$(sed -n ${number3}'p' $genetic_locations)
        output_path=${output_path}/${genome_name}
        
        echo "parameters:$genome_name
                            ${genome_fasta_path}
                            deduplicate: $deduplicate"   
        echo "parameters:$genome_name
                            ${genome_fasta_path}
                            deduplicate: $deduplicate
                            " >> $parameter_file             
    done
        mkdir -p $output_path 
        echo "output path: $output_path"

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
        echo ""
        echo "conversion of bcls to fastqs complete"
        echo ""
    else
        echo ""
        echo "bcls already transformed into fastqs"
        echo ""
    fi

##################################################################################################################################
##############################################---STEP 5: RUN methylseq.sh---###################################################### 
##################################################################################################################################    
    fastq_file="${data_path}/fastq/FASTQs"                                              #give a path to a file to store the fastq file paths in $fastq_directory
    echo "location of fastq_file: $fastq_file"
    find "$data_path/fastq" -type f | grep ".*\.fastq.gz$" | grep -v ".*\.trimmed.fastq.gz$" | sed -e 's/_R1.*$//g' | sed -e 's/_R2.*$//g' | sort -u > "${fastq_file}"                        
                                                                                        #generate list of full paths to fastq files and save to the file in $fastq_list
    array_length=$(wc -l < "${fastq_file}")                                             #get the number of files and, thus, array length
    echo "array length: $array_length
    " >> $parameter_file
   
    picard=$(find "$output_path" -type f | grep ".*\.bam_picard_insert_size_plot.pdf$" | sort -u | wc -l)
    #this should actually check for coverage output file from bismark, as that's the real important output file, not picard.

    if [[ $picard -lt 1 ]] || [[ $force = true ]]; then
            echo ""
            echo "entering methylseq.sh"
            echo ""
            sbatch -o "$Logs/${log_name}_%A_%a.log" `#put into log` \
                    -a "1-${array_length}" `#initiate job array equal to the number of fastq files` \
                   -W `#indicates to the script not to move on until the sbatch operation is complete` \
                   --cpus-per-task "${cores}" \
                   --mem "$((cores * 16 / 3))GB" `#the ratio of memory needed for bismark to work` \
                    "${code_directory}/methylseq.sh" \
                    $data_path $output_path $genetic_locations $cores $log_name $parameter_file $code_directory $Logs $parameter_file $initial_path

            wait
            echo ""
            echo "methylseq.sh complete"
            echo ""
        else
            echo ""
            echo "picard files already created, methylseq.sh skipped"
            echo ""
    fi

#####################previously report_controls.sh################################
################################################################################
    previous_loop_output_directory="${initial_path}"

    for i in $(seq 0 $total_genomes); do
        number1=$(bc -l <<< "scale=0; (($i * 3) +1)")
        genome_name=$(sed -n ${number1}'p' $genetic_locations)

        output_directory="${previous_loop_output_directory}/${genome_name}"
        bismark_summary="${output_directory}/bismark_summary_report.txt" 

        if [ ! -f $bismark_summary ]; then
            cd $output_directory
                echo "$bismark_summary does not exist yet"
            module load bismark
            bismark2report
            bismark2summary
                echo ""
                echo "bismark summary now complete"
                echo ""
        else
            echo ""
            echo "bismark summary already complete, skipping bismark2report and bismark2summary for $genome_name"
            echo ""
        fi

        #set previous_loop_output_directory variable for next possible loop iteration
        previous_loop_output_directory="${output_directory}"
    done

fi
