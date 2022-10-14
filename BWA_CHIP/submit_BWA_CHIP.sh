#!/bin/bash
##################################################################################################################################
#############################################---STEP 1: SET UP PARAMETERS---###################################################### 
##################################################################################################################################
 if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ]; then
    echo "Format: ./submit_BWA_CHIP.sh [fastq_directory] [output_directory] --argument"
    echo "fastq_directory: path to raw .fastq or .fastq.gz files"
    echo "output_directory: path for BWA and mutect output"
    echo "argument: indicates requested analysis technique(s) (--mutect, --varscan, --haplotypecaller, or --all)."
    echo "If --varscan is selected, you may use the optional arguments --p_value, --min_var_freq, and --min_coverage (if so desired)." 
    echo "If --mutect is selected, you may use --twist to indicate you would like your results remove_silent by the Twist panel"
    #exit 1
else

    TEMP=`getopt -o vdm: --long min_coverage:,input:,output:,working_dir:,min_var_freq:,p_value:,intervals:,normal_sample:,log_name:,reference_genome:,panel:,assembly:,funcotator_sources:,transcript_list:,mode:,docker_image:,container_engine:,sequence_dictionary:,chr_intervals:,normal_pileups:,n_jobs:,gnomad_genomes:,bam,remove_silent,mutect,varscan,haplotypecaller,all,skip_funcotator,no_bam_out,bam,fastq,realign,normal_pileups,split_by_chr \
    -n 'submit_BWA_CHIP.sh' -- "$@"`

    code_directory=$(realpath .)  #specify location of star_align_and_qc.sh
    echo "CODE_DIRECTORY: $code_directory"
    
    #input, output, and working directory using -- 
    #input and output for slurm
    #working for cloud
    #add in failsafes and checks
    #all refrence data must be within the reference folder, inputs in inputs, and outputs will be made into a folder called utputs
        
   if [ $? != 0 ]; then
       echo "Unrecognized argument. Possible arguments: mutect, varscan, haplotypecaller, all, min_coverage, min_var_freq, p_value, and twist." >&2 ; exit 1 ; 
   fi
       eval set -- "$TEMP"

        get_mutect=false
        get_varscan=false
        get_haplotype=false
        all=false
        normal_sample=false
        run_funcotator=true
        bam_out=true
        realign=false
        remove_silent="0"
        log_name="log_"
        #twist_snps="/labs/sjaiswal/workflows/BWA_mutect_twist/twist_snps.bed"
        assembly="GRCh38"
        calculate_contamination=false
        panel=false
        min_coverage="10"
        min_var_freq="0.001"
        p_value="0.1"
        bam=false
        fastq=false
        normal_pileups=false
        mode="slurm"
        docker_image="path"
        container_engine="path"
        split_by_chr=false
        n_jobs=1
        working_directory=false
        data_directory=false
        output_directory=false
        
    while true; do
        case "$1" in
            --input | --input_dir | --input_directory ) data_directory="$2"; shift 2 ;;
            --output | --output_dir | --output_directory ) output_directory="$2"; shift 2 ;; 
            --working_dir ) working_directory="$2"; shift 2;;
            --min_coverage ) min_coverage="$2"; shift 2 ;;
            --min_var_freq ) min_var_freq="$2"; shift 2 ;;
            --p_value ) p_value="$2"; shift 2 ;;
            --normal_sample ) normal_sample="$2"; shift 2 ;; 
            --log_name ) log_name="$2"; shift 2;;
            --panel ) panel="$2"; shift 2;;
            --assembly ) assembly="$2"; shift 2;;
            --normal_pileups ) normal_pileups="$2"; shift 2 ;;
            --mode ) mode="$2"; shift 2 ;;
            --docker_image ) docker_image="$2"; shift 2 ;;
            --container_engine ) container_engine="$2"; shift 2 ;;
            --n_jobs ) n_jobs="$2"; shift 2 ;;
            --remove_silent ) remove_silent="1"; shift ;; 
            -m | --mutect ) get_mutect=true; shift ;;
            -v | --varscan ) get_varscan=true; shift ;;
            -h | --haplotypecaller ) get_haplotype=true; shift ;;
            --skip_funcotator ) run_funcotator=false; shift ;;
            --no_bam_out ) bam_out=false; shift ;;
            --all ) get_mutect=true; get_varscan=true; get_haplotype=true; shift ;;
            --realign ) realign=true; shift ;;
            --calculate_contamination ) calculate_contamination=true; shift ;;
            --bam ) bam=true; shift ;; 
            --fastq ) fastq=true; shift ;;
            --whole_genome | -WG | --Whole_Genome | --split_by_chr | --split ) split_by_chr=true; shift ;;
            -- ) shift; break ;;
            * ) break ;;
        esac
    done
    
    if [[ ! -f "/oak/stanford/groups/smontgom/maurertm/ADRC/Cloud_Testing_Folder/variant_calling/workflows/BWA_CHIP/Params/GRCh38.p12.genome.u2af1l5_mask.fa" ]]; then
        echo "THE FILE ISN'T THERE"
    fi

    echo "realign: $realign"
    echo "###################MODE##################: $mode"

    if [[ $calculate_contamination = true ]] && [[ $normal_pileups = false ]] && [[ $normal_sample != false ]]; then
        echo "if contamination calculation is requested while using tumor normal mode, an additonal argument containing --normal_pileups /path/to/pileups_file.pileups must be provided"
    fi

    if [[ $bam = false ]] && [[ $fastq = false ]]; then
        echo "please use either --bam or --fastq to indicate the file type of your data"
        echo "if your bam data is not aligned to hg38, please use the following arguments: --bam --realign"
        exit 1
    fi

    if [[ $bam = false ]] && [[ $realign = true ]]; then
        echo "all fastq files are aligned to hg38"
        echo "please remove --relign from your argument list. It cannot be used with --fastq"
        exit 1
    fi

    echo "intervals: $intervals"
    if ( [[ $min_coverage != "10" ]] || [[ $min_var_freq != "0.001" ]] || [[ $p_value != "0.1" ]] ) && \
        ( [[ $get_varscan = false ]] || [[ $all = false ]] ); then
        echo "The p_value, min_coverage, and min_var_freq arguments are not used in the mutect and haplotypecaller workflows"; 
        exit 1
    fi

    if ( [[ $panel != false ]] ) && \
        ( [[ $get_mutect = false ]] || [[ $all = false ]] ); then
        echo "The --panel argument is not applicable to the varscan and haplotypecaller workflows"; 
        exit 1
    fi

    if [[ $mode == "slurm" ]]; then
        echo "You have selected to run the BWA_CHIP pipeline on a slurm job management system."
        echo "If you did not mean to make this selection, please cancel this job and run with any of the following options:"
        echo "--mode "docker" to run on the cloud from a docker file OR"
        echo "--mode "cloud" to run on cloud from a docker image"
    fi

    if [[ $mode != "slurm" ]] && [[ $mode != "cloud" ]]; then
        echo "You did not select a recognizable option for mode."
        echo "Please select --mode cloud or --mode slurm"
    fi

    #data_directory=$1 #get directory path from second argument (first argument $0 is the path of this script)
    #find "${data_directory}/" -type f | grep "bam" | grep -v ".bam.bai" | sed -e 's/\.bam$//g'
    #output_directory=$2
    if [[ $data_directory == false ]]; then
        data_directory="${working_directory}/Inputs" 
    fi
    if [[ $output_directory == false ]]; then
        output_directory="${working_directory}/Outputs"
    fi
    

    echo "################### input - $data_directory ####### output - $output_directory ##### working - $working_directory"



    parent_directory=$(dirname $data_directory) #get parent directory of $fastq_directory
    fastq_list="${parent_directory}/fastq_files" #give a path to a file to store the paths to the fastq files in $fastq_directory
    bam_list="${parent_directory}/bam_files"
    data_count=${#data_directory}
    output_count=${#output_directory}
    sum=$(($data_count + $output_count + 9))
    TEMPS=${TEMP[0]} 
    TEMP_ARGUMENTS=${TEMPS::-$sum}
    echo "$TEMP_ARGUMENTS"  

    source "${code_directory}/config.sh"
    echo "$reference_genome"

    if [[ -z "$reference_genome" ]]; then
        echo "You must provide a path to the reference_genome in config.sh"
        exit 1
    elif [[ -z "$gnomad_genomes" ]]; then
        echo "You must provide a path to the gnomad_genomes in config.sh"
        exit 1
    fi

    if [[ $run_funcotator == true ]] && [[ -z "$funcotator_sources" ]]; then
        echo "You did not use the flag --skip_funcotator, which means you would like to run funcotator"
        echo "In order to run funcotator, you must include a path to the folder containing the funcotator sources in config.sh"
        exit 1
    fi

    if [[ $split_by_chr == true ]] && [[ -z "$sequence_dictionary" ]]; then
        echo "You used one of the following flags: --whole_genome, -WG, --Whole_Genome, --split_by_chr, --split"
        echo "These flags indicate that you would like you data to be split up by its chromosomes, which will  increase speed"
        echo "In order to split the data into chromosomal segments, you must provide a path to the sequence dictionary in config.sh"
        exit 1
    fi

    if [[ -z "$chr_intervals" ]]; then
        chr_intervals="${code_directory}/whole_genome_intervals.interval_list"
    fi
    
    if [[ -z "$intervals" ]]; then
        intervals="${code_directory}/CHIP_exons.interval_list"
    fi


    if [[ ! -d $data_directory ]]; then
        echo "data directory does not exist"
    elif [[ ! -f $reference_genome ]]; then
        echo "refrence genome path does not exist"
    elif [[ ! -f $funcotator_sources ]]; then
        echo "funcotator sources path does not exist"
    elif [[ ! -f $transcript_list ]]; then
        echo "transcript_list path dos not exist"
    elif [[ ! -f $chr_intervals ]]; then
        echo "chr_intervals path does not exist"
    elif [[ ! -f $intervals ]]; then
        echo "intervals path does not exist"
    fi


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
    if [ $log_name = "log_" ]; then                                                             #give a path to a file to store the parameter files (so they are unique)
        parameter_file="$Parameters/${now}_parameters.txt"                                 #add date stamp to parameter files and, if provided, the log name
    else
        parameter_file="$Parameters/${log_name}${now}_parameters.txt"
    fi
    touch $parameter_file
        
    echo "location of scripts used to run code : $code_directory
        " > $parameter_file 
        
    
    echo "call made to execute code: $0 $1 $2 $TEMP_ARGUMENTS 
    " >> $parameter_file
    
    echo "you chose the following optional arguments: 
        $TEMP_ARGUMENTS    
        " >> $parameter_file
        
##################################################################################################################################
#############################################--STEP 4: GET NORMAL PILEUPS---####################################################### 
##################################################################################################################################
    if [[ $normal_sample != false ]]; then
        #singularity
                
        NORMAL_SAMPLE_BASENAME=$(basename "${normal_sample}")
        NORMAL_SAMPLE_DIRNAME=$(dirname "${normal_sample}")
        NORMAL_SAMPLE_NAME=${NORMAL_SAMPLE_BASENAME%.*}

        array_length=1
        line_number=1
        run_mutect=false
        normal_sample_path="${parent_directory}/normal_sample_path"

        find -L "${NORMAL_SAMPLE_DIRNAME}" -type f `#list all files in ${fastq_directory}` | \
                grep "$NORMAL_SAMPLE_BASENAME" `#only keep files with FASTQ in name (case insensitive)` > $normal_sample_path
         
        
        if [[ $mode == "slurm" ]]; then
            sbatch -o "${output_directory}/Logs/%A_%a.log" `#put into log` \
                -a "1-${array_length}" `#initiate job array equal to the number of fastq files` \
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
                ${bam} \
                ${intervals} \
                ${normal_sample} \
                ${code_directory} \
                ${parameter_file} \
                ${reference_genome} \
                ${panel} \
                ${assembly} \
                ${funcotator_sources} \
                ${transcript_list} \
                ${remove_silent} \
                ${run_funcotator} \
                ${bam_out} \
                ${normal_pileups} \
                ${mode} \
                ${docker_image} \
                ${container_engine} \
                ${split_by_chr} \
                ${sequence_dictionary} \
                ${chr_intervals} \
                ${gnomad_genomes} \
                ${run_mutect} \
                ${working_dir}
                wait
        else
            echo "CODE_DIRECTORY: $code_directory"
            TASK_ID=1 ${code_directory}/BWA_CHIP.sh \
                ${parent_directory} \
                ${output_directory} \
                ${min_coverage} \
                ${min_var_freq} \
                ${p_value} \
                ${get_mutect} \
                ${get_varscan} \
                ${get_haplotype} \
                ${bam} \
                ${intervals} \
                ${normal_sample} \
                ${code_directory} \
                ${parameter_file} \
                ${reference_genome} \
                ${panel} \
                ${assembly} \
                ${funcotator_sources} \
                ${transcript_list} \
                ${remove_silent} \
                ${run_funcotator} \
                ${bam_out} \
                ${normal_pileups} \
                ${mode} \
                ${docker_image} \
                ${container_engine} \
                ${split_by_chr} \
                ${sequence_dictionary} \
                ${chr_intervals} \
                ${gnomad_genomes} \
                ${run_mutect} \
                ${working_dir}
        fi

            normal_pileups="${output_directory}/NORMAL_PILEUPS/${NORMAL_SAMPLE_NAME}_${assembly}_pileups.table" 
            echo "$normal_pileups"
            
            if [ ! -f ${normal_pileups} ]; then
                echo "${normal_pileups} does not exist"
            fi
    fi

##################################################################################################################################
#############################################--STEP 4: GET ARRAY LENGTHS---####################################################### 
##################################################################################################################################
    run_mutect=true
    if [ $bam = true ]; then
        find -L "${data_directory}" -type f `#list all files in ${fastq_directory}` | \
                    grep ".*\.bam$" `#only keep files with FASTQ in name (case insensitive)` | \
                    grep -v ".*\.bai$" `#remove Undetermined FASTQs` | \
                    sed -e 's/\.bam$//g' `#remove _R1/2_fastq.gz file extension` | \
                    sort -u  `#sort and remove duplicate names` > ${bam_list}
                        #| \
                    #head -n -1 > ${bam_list} `#remove the last line and generate a list of unique FASTQs`
            bam_array_length=$(wc -l < ${bam_list}) #get the number of FASTQs
            echo "bam array length: $bam_array_length"

        if [ $realign = true ]; then 
            echo "entering bam_to_fastq.sh script"

            if [[ $mode == "slurm" ]]; then
                sbatch -o "${output_directory}/Logs/%A_%a.log" `#put into log` \
                    -a "1-${bam_array_length}" `#initiate job array equal to the number of fastq files` \
                    -W `#indicates to the script not to move on until the sbatch operation is complete` \
                    "${code_directory}/bam_to_fastq.sh" \
                    "$data_directory" \
                    "$code_directory"
            fi

            cp $bam_list $fastq_list
            array_length=$(wc -l < ${fastq_list}) #get the number of FASTQs 
            echo "array length: $array_length" 
            $bam=false

        else 
            array_length=$bam_array_length
        fi
    else
       find -L "${data_directory}" -type f `#list all files in ${fastq_directory}` | \
            grep ".*\.fastq.gz$" `#only keep files with FASTQ in name (case insensitive)` | \
            grep -v "Undetermined" `#remove Undetermined FASTQs` | \
            sed -e 's/_R1.*$//g' | sed -e 's/_R2.*$//g' `#remove _R1/2_fastq.gz file extension`| \
            sort -u  `#sort and remove duplicate names` > ${fastq_list} 
                #| \
            #head -n -1 > ${fastq_list} `#remove the last line and generate a list of unique FASTQs`
        array_length=$(wc -l < ${fastq_list}) #get the number of FASTQs 
        echo "fastq array length: $array_length" 
    fi

##################################################################################################################################
##################################################--STEP 4: BWA_CHIP.sh---######################################################## 
##################################################################################################################################
    echo "sbatch -o "${output_directory}/Logs/%A_%a.log" `#put into log` \
    -a "1-${array_length}" `#initiate job array equal to the number of fastq files` \
    -W `#indicates to the script not to move on until the sbatch operation is complete` \
    "${code_directory}/BWA_CHIP.sh" \
    parent_directory ${parent_directory} \
    output_directory ${output_directory} \
    min_coverage ${min_coverage} \
    min_var_freq ${min_var_freq} \
    p_value ${p_value} \
    get_mutect${get_mutect} \
    get_varscan ${get_varscan} \
    get_haplotype ${get_haplotype} \
    bam ${bam} \
    intervals ${intervals} \
    normal_sample ${normal_sample} \
    code_directory ${code_directory} \
    parameter_file ${parameter_file} \
    reference_genome ${reference_genome} \
    panel  ${panel} \
    assembly ${assembly} \
    funcotato sources ${funcotator_sources} \
    transcript lists ${transcript_list} \
    remove_silent ${remove_silent} \
    run funcotator ${run_funcotator} \
    bam out ${bam_out} \
    normal_pileups ${normal_pileups} \
    calculate_contamination ${calculate_contamination} \
    mode ${mode} \
    docker_image ${docker_image} \
    container_image ${container_engine} \
    split_by_chr ${split_by_chr} \
    sequence_dictionary ${sequence_dictionary} \
    chr_intervals ${chr_intervals} \
    gnomad_genomes ${gnomad_genomes} \
    run_mutect ${run_mutect}"         

    if [[ $mode == "slurm" ]]; then 
        sbatch -o "${output_directory}/Logs/%A_%a.log" `#put into log` \
        -a "1-${array_length}" `#initiate job array equal to the number of fastq files` \
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
        ${bam} \
        ${intervals} \
        ${normal_sample} \
        ${code_directory} \
        ${parameter_file} \
        ${reference_genome} \
        ${panel} \
        ${assembly} \
        ${funcotator_sources} \
        ${transcript_list} \
        ${remove_silent} \
        ${run_funcotator} \
        ${bam_out} \
        ${normal_pileups} \
        ${mode} \
        ${docker_image} \
        ${container_engine} \
        ${split_by_chr} \
        ${sequence_dictionary} \
        ${chr_intervals} \
        ${gnomad_genomes} \
        ${run_mutect} \
        ${working_dir}
    else
        echo "CODE DIRECTORY: $code_directory"
        #. `which env_parallel.bash`
        seq 1 ${array_length} | parallel --ungroup --progress -j ${n_jobs} TASK_ID={} ${code_directory}/BWA_CHIP.sh \
        ${parent_directory} \
        ${output_directory} \
        ${min_coverage} \
        ${min_var_freq} \
        ${p_value} \
        ${get_mutect} \
        ${get_varscan} \
        ${get_haplotype} \
        ${bam} \
        ${intervals} \
        ${normal_sample} \
        ${code_directory} \
        ${parameter_file} \
        ${reference_genome} \
        ${panel} \
        ${assembly} \
        ${funcotator_sources} \
        ${transcript_list} \
        ${remove_silent} \
        ${run_funcotator} \
        ${bam_out} \
        ${normal_pileups} \
        ${mode} \
        ${docker_image} \
        ${container_engine} \
        ${split_by_chr} \
        ${sequence_dictionary} \
        ${chr_intervals} \
        ${gnomad_genomes} \
        ${run_mutect} \
        ${working_dir} 
        echo "Test"
    fi
    wait
##################################################################################################################################
################################################--STEP 4: RUN R SCRIPTS---######################################################## 
##################################################################################################################################
    if [[ $mode = "slurm" ]]; then
        module load R/4.0

        if [ $get_mutect = true ]; then
            Rscript aggregate_variants_mutect.R /labs/sjaiswal/chip_submitted_targets_Twist.xls \
                "$output_directory" "$panel" "$remove_silent" > "$output_directory/Logs/mutectOutFile.Rout" 2>&1
                
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
    fi
fi
