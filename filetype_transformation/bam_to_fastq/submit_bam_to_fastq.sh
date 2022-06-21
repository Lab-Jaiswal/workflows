##################################################################################################################################
#############################################---STEP 1: SET UP PARAMETERS---###################################################### 
##################################################################################################################################
if [ -z $1 ] || [ -z $2 ]; then
    echo "Format: ./submit_bam_to_fastq.sh [data_directory] [output_path]"
    exit 1
else
    data_directory=$1
    output_path=$2

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

    Logs="${output_path}/Logs"
    Parameters="${output_path}/Parameters"

##################################################################################################################################
#######################################---STEP 3: CREATE PARAMETER LOG---#########################################################
##################################################################################################################################
    now=$(date +%m_%d_%H_%M)
    
    parameter_file="$Parameters/${now}_parameters.txt"                                 #add date stamp to parameter files and, if provided, the log name
    
    touch $parameter_file
        
    echo "location of scripts used to run code : $code_directory
        " > $parameter_file 
        
    echo "call made to execute code: $0 $1 $2
    " >> $parameter_file
            
##################################################################################################################################
##########################################---STEP 5: CONVERT BAMS TO FASTQs---#################################################### 
##################################################################################################################################
    cd $data_directory 
    
    bam_list="$data_directory/bam_files" #give a path to a file to store the paths to the fastq files in $fastq_directory

   find -L "${data_directory}/" -type f `#list all files in ${fastq_directory}` | \
        grep ".*\.bam$" `#only keep files with FASTQ in name (case insensitive)` | \
        grep -v ".bai" `#remove Undetermined FASTQs` | \
        sed -e 's/\.bam$//g' `#remove _R1/2_fastq.gz file extension` | \
        sort -u  `#sort and remove duplicate names` > ${bam_list}
    bam_array_length=$(wc -l < ${bam_list}) #get the number of FASTQs 
    echo "bam array length: $bam_array_length"
     
    sbatch -o "$Logs/%A_%a.log" `#put into log` \
            -a "1-${bam_array_length}" `#initiate job array equal to the number of fastq files` \
            -W `#indicates to the script not to move on until the sbatch operation is complete` \
            "${code_directory}/bam_to_fastq.sh" \
            "$data_directory" \
            "$output_path"

    echo "all bams have been converted to fastqs"

fi
