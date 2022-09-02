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
    exit 1
else
    TEMP=`getopt -o vdm: --long min_coverage:,min_var_freq:,p_value:,intervals:,normal_sample:,log_name:,reference_genome:,panel:,assembly:,funcotator_sources:,transcript_list:,mode:,docker_image:,container_engine:,sequence_dictionary:,chr_intervals:,bam,remove_silent,mutect,varscan,haplotypecaller,all,skip_funcotator,no_bam_out,bam,fastq,realign,normal_pileups,split_by_chr \
    -n 'submit_BWA_CHIP.sh' -- "$@"`

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
        reference_genome="/oak/stanford/groups/sjaiswal/Herra/CHIP_Panel_AmpliSeq/GRCh38.p12.genome.u2af1l5_mask.fa"
        #twist_snps="/labs/sjaiswal/workflows/BWA_mutect_twist/twist_snps.bed"
        assembly="GRCh38"
        funcotator_sources="/labs/sjaiswal/tools/funcotator/funcotator_dataSources.v1.6.20190124s"
        #funcotator_sources="/home/maurertm/smontgom/maurertm/funcotator_dataSources.v1.6.20190124s"
        transcript_list="/oak/stanford/groups/sjaiswal/Herra/CHIP_TWIST-PANEL_ATHEROMA/chip_transcript_list.txt"
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
        split_by_chr=false"/oak/stanford/groups/sjaiswal/Herra/CHIP_Panel_AmpliSeq/GRCh38.p12.genome.u2af1l5_mask.fa.dict"
        sequence_dictionary="/oak/stanford/groups/sjaiswal/Herra/CHIP_Panel_AmpliSeq/GRCh38.p12.genome.u2af1l5_mask.fa.dict"
        code_directory=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )  #specify location of star_align_and_qc.sh
        chr_intervals="${code_directory}/whole_genome_intervals.interval_list"
        intervals="${code_directory}/CHIP_exons.interval_list"
        
    while true; do
        case "$1" in
            --min_coverage ) min_coverage="$2"; shift 2 ;;
            --min_var_freq ) min_var_freq="$2"; shift 2 ;;
            --p_value ) p_value="$2"; shift 2 ;;
            --intervals ) intervals="$2"; shift 2 ;;
            --normal_sample ) normal_sample="$2"; shift 2 ;; 
            --log_name ) log_name="$2"; shift 2;;
            --reference_genome ) reference_genome="$2"; shift 2 ;;
            --panel ) panel="$2"; shift 2;;
            --assembly ) assembly="$2"; shift 2;;
            --funcotator_sources ) funcotator_sources="$2"; shift 2 ;;
            --transcript_list ) transcript_list="$2"; shift 2 ;;
            --normal_pileups ) normal_pileups="$2"; shift 2 ;;
            --mode ) mode ="$2"; shift 2 ;;
            --docker_image ) docker_image="$2"; shift 2 ;;
            --container_engine ) container_engine="$2"; shift 2 ;;
            --chr_intervals | --whole_genome_intervals | --WG_intervals ) chr_intervals="$2"; shift 2 ;;
            --sequence_dictionary | --sequence_dict ) sequence_dictionary="$2"; shift 2 ;;
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

    echo "realign: $realign"

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

    if ( [[ $mode != "slurm" ]] || [[ $mode != "cloud" ]] ); then
        echo "You did not select a recognizable option for mode."
        echo "Please select --mode cloud or --mode slurm"
    fi

    data_directory=$1 #get directory path from second argument (first argument $0 is the path of this script)
    #find "${data_directory}/" -type f | grep "bam" | grep -v ".bam.bai" | sed -e 's/\.bam$//g'
    output_directory=$2
    parent_directory=$(dirname $data_directory) #get parent directory of $fastq_directory
    fastq_list="${parent_directory}/fastq_files" #give a path to a file to store the paths to the fastq files in $fastq_directory
    bam_list="${parent_directory}/bam_files"
    data_count=${#data_directory}
    output_count=${#output_directory}
    sum=$(($data_count + $output_count + 9))
    TEMPS=${TEMP[0]} 
    TEMP_ARGUMENTS=${TEMPS::-$sum}
    echo "$TEMP_ARGUMENTS"

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
if ( [[ $normal_sample != false ]] && [[ $normal_pileups != false ]] ); then
            
    NORMAL_SAMPLE_BASENAME=$(basename $NORMAL_SAMPLE)
    NORMAL_SAMPLE_NAME=${SAMPLE_BASENAME%.*}

    if [[ $split_by_chr = true ]]; then
            
            #CHR_INTERVALS="/oak/stanford/groups/smontgom/maurertm/ADRC/workflows/Interval_filtering/whole_genome_intervals.interval_list"
            num_intervals=$(grep -v "@" < $chr_intervals | wc -l)
            
            if [[ $mode = "slurm" ]]; then
                module load bcftools
                module load gatk4 
     
                seq 1 ${num_intervals} | parallel -j8 --progress --ungroup "${CODE_DIRECTORY}/calculate_pileups.sh $normal_sample $intervals $parameter_file $code_directory $output_directory $split_by_chr {}"
            
            else
                echo "some other command"
            fi

            pileup_tables=$(find $output_directory/Normal_Pileups/pileups_by_chr -type f | sed -e 's/^/-I /g' | tr '\n' ' ')
            
            gatk GatherPileupSummaries --sequence-dictionary $sequence_dictionary ${pileup_tables} -O ${output_directory}/Normal_Pileups/${NORMAL_SAMPLE_NAME}_pileups.table

        else
            if [[ $mode = "slurm" ]]; then
                array_length=1

                sbatch -o "${output_directory}/Logs/%A_%a.log" `#put into log` \
                -a "1-${array_length}" `#initiate job array equal to the number of fastq files` \
                -W `#indicates to the script not to move on until the sbatch operation is complete` \
                "${code_directory}/calculate_pileups.sh" \
                ${normal_sample} \
                ${intervals} \
                ${parameter_file} \
                ${code_directory} \
                ${output_directory} \
                ${split_by_chr}

                wait
            else
                echo "some other command"
            fi

        fi

            normal_pileups="${output_directory}/Normal_Pileups/${NORMAL_SAMPLE_NAME}_pileups.table"    
fi

##################################################################################################################################
#############################################--STEP 4: GET ARRAY LENGTHS---####################################################### 
##################################################################################################################################

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
    if [[ $mode == "slurm" ]]; then 
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
        chr_intervals ${chr_intervals}"         

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
        ${calculate_contamination} \
        ${mode} \
        ${docker_image} \
        ${container_engine} \
        ${split_by_chr} \
        ${sequence_dictionary} \
        ${chr_intervals}

        wait

    else
        #gnu parallele command
        echo "Test"
    fi
fi
    
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
