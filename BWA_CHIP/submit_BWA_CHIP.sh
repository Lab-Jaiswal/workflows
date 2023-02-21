#!/usr/bin/env bash

# Run scripts to enable activating conda environments
. "${HOME}/micromamba/etc/profile.d/conda.sh"
. "${HOME}/micromamba/etc/profile.d/mamba.sh"
mamba activate base

# Set bash options for verbose output and to fail immediately on errors or if variables are undefined.
set -o xtrace -o nounset -o pipefail -o errexit

check_for_file() {
    argument_name="${1}"
    file_path="${2}"
    if [[ ${file_path} != "none" ]] && [[ ! -f ${file_path} ]]; then
        echo "Error: file ${file_path} passed with ${argument_name} does not exist."
        exit 1
    fi
}

check_for_directory() {
    argument_name="${1}"
    directory_path="${2}"
    if [[ ${directory_path} != "none" ]] && [[ ! -d ${directory_path} ]]; then
        echo "Error: directory ${directory_path} passed with ${argument_name} does not exist."
        exit 1
    fi
}

##################################################################################################################################
#############################################---Step 1: set up parameters---######################################################
##################################################################################################################################

options_array=(
    config_file
    input_directory
    input_file_list
    output_directory
    bam_extension
    fastq_extension
    assembly
    reference_genome
    slurm_runtime
    slurm_account
    slurm_ncpus
    slurm_memory
    slurm_jobname
    interval_list
    sequence_dictionary
    gnomad_reference
    normal_bam
    normal_pileups_table
    funcotator_sources
    transcript_list
    mpileup_interval_bed
    varscan_min_coverage
    varscan_min_var_freq
    varscan_max_pvalue
    annovarroot
    germline_snps
    n_jobs
)

flags_array=(
    slurm_mode
    run_mutect
    run_varscan
    run_haplotypecaller
    all
    mutect_bam_output
    run_annovar
    run_funcotator
    split_intervals
    realign
)

eval "$(printf "%s\n" "${options_array[@]}" | xargs --replace=% echo "declare %=none;")"
eval "$(printf "%s\n" "${flags_array[@]}" | xargs --replace=% echo "declare %=false;")"

options_joined=$(echo "${options_array[@]}" | sed -e 's/ /:,/g')
flags_joined=$(echo "${flags_array[@]}" | tr ' ' ',')
longoptions="${options_joined}:,${flags_joined}"

arguments=$(getopt --options a --longoptions "${longoptions}" --name 'submit_BWA_CHIP.sh' -- "$@")
eval set -- "${arguments}"

while true; do
    case "$1" in
        --config_file )
            config_file="$2"; check_for_file "${1}" "${2}"; shift 2 ;;
        -- )
            shift; break ;;
        * )
            shift; 
    esac
done

if [[ ${config_file} != "none" ]]; then
    source "${config_file}"
fi

eval set -- "${arguments}"
while true; do
    case "$1" in
        --input_directory )
            input_directory="$2"; check_for_directory "${1}" "${2}"; shift 2 ;;
        --input_file_list )
            input_file_list="$2"; check_for_file "${1}" "${2}"; shift 2 ;;
        --output_directory )
            output_directory="$2"; shift 2 ;;
        --bam_extension )
            bam_extension="$2"; shift 2 ;;
        --fastq_extension )
            fastq_extension="$2"; shift 2 ;;
        --assembly )
            assembly="$2"; shift 2 ;;
        --reference_genome )
            reference_genome="$2"; check_for_file "${1}" "${2}"; shift 2 ;;
        --slurm_mode )
            slurm_mode=true; shift ;;
        --slurm_runtime )
            slurm_runtime="$2"; shift 2 ;;
        --slurm_account )
            slurm_account="$2"; shift 2 ;;
        --slurm_ncpus )
            slurm_ncpus="$2"; shift 2 ;;
        --slurm_memory )
            slurm_memory="$2"; shift 2 ;;
        --slurm_jobname )
            slurm_jobname="$2"; shift 2 ;;
        --mutect )
            run_mutect=true; shift ;;
        --interval_list )
            interval_list="$2"; check_for_file "${1}" "${2}"; shift 2 ;;
        --split_intervals )
            split_intervals=true; shift ;;
        --sequence_dictionary )
            sequence_dictionary="$2"; check_for_file "${1}" "${2}"; shift 2 ;;
        --mutect_bam_output )
            mutect_bam_output=true; shift ;;
        --gnomad_reference )
            gnomad_reference="$2"; check_for_file "${1}" "${2}"; shift 2 ;;
        --normal_bam )
            normal_bam="$2"; check_for_file "${1}" "${2}"; shift 2 ;;
        --normal_pileups_table )
            normal_pileups_table="$2"; check_for_file "${1}" "${2}"; shift 2 ;;
        --skip_funcotator )
            run_funcotator=false; shift ;;
        --funcotator_sources )
            funcotator_sources="$2"; check_for_directory "${1}" "${2}"; shift 2 ;;
        --transcript_list )
            transcript_list="$2"; check_for_file "${1}" "${2}"; shift 2 ;;
        --varscan )
            run_varscan=true; shift ;;
        --mpileup_interval_bed )
            mpileup_interval_bed="$2"; check_for_file "${1}" "${2}"; shift 2 ;;
        --varscan_min_coverage )
            varscan_min_coverage="$2"; shift 2 ;;
        --varscan_min_var_freq )
            varscan_min_var_freq="$2"; shift 2 ;;
        --varscan_max_pvalue )
            varscan_max_pvalue="$2"; shift 2 ;;
        --skip_annovar )
            run_annovar=false; shift ;;
        --annnovarroot )
            annovarroot="$2"; check_for_directory "${1}" "${2}"; shift 2 ;;
        --haplotypecaller )
            run_haplotypecaller=true; shift ;;
        --germline_snps )
            germline_snps="$2"; check_for_file "${1}" "${2}"; shift 2 ;;
        --all )
            run_mutect=true; run_haplotypecaller=true; run_varscan=true; shift ;;
        --n_jobs )
            n_jobs="$2"; shift 2 ;;
        --realign )
            realign=true; shift ;;
        -- )
            shift; break ;;
        * )
            shift; break ;;
            #echo "Invalid argument ${1}" >&2
            #exit 1
    esac
done

#if [[ -z ${input_directory}  ]] && [[ -z ${input_file_list} ]]; then
    #echo "Error: an input directory or input file list must be provided as arguments (--input_directory [input directory] or --input_file_list [input file list]) or defined in a config file (input_directory=[input directory] or input_file_list=[input file list])"
    #exit 1
#fi

#if [[ -z ${output_directory} ]] && [[ ${slurm_mode} == true ]]; then
    #echo "Warning: running in SLURM mode (--slurm_mode was provided as an argument or slurm_mode was set to true in a config file) without a dedicated output directory (--output_directory [output directory] was not provided as an argument and output_directory=[output directory] was not defined in a config file).  All output will be sent to $HOME, which usually has a low storage quota on HPC systems."
#fi

#if [[ -z ${bam_extension} ]] && [[ -z ${fastq_extension} ]]; then
    #echo "Error: a file extension must be provided as an argument (--bam_extension [BAM file extension] or --fastq_extension [FASTQ file extension]) or defined in a config file (bam_extension=[BAM file extension] or fastq_extension=[FASTQ file extension])"
    #exit 1
#fi

#if [[ -n ${bam_extension} ]] && [[ -n ${fastq_extension} ]]; then
    #echo "Error: only a file extension for BAM/CRAM or FASTQ files may be provided as an argument (--bam_extension [BAM file extension] or --fastq_extension [FASTQ file extension]) or defined in a config file (bam_extension=[BAM file extension] or fastq_extension=[FASTQ file extension])"
    #exit 1
#fi

#if [[ ${bam_extension} != "bam" ]] && [[ ${bam_extension} != "cram" ]] && [[ -n ${bam_extension} ]]; then
    #echo "Error: Valid values for the file extension for a BAM/CRAM provided as an argument (--bam_extension [BAM file extension]) or defined in a config file (bam_extension=[BAM extension]) are \"bam\" and \"cram\"."
    #exit 1
#fi

#if [[ -z ${assembly} ]] && [[ -n ${fastq_extension} ]]; then
    #echo "Error: a reference genome assembly (e.g. GRCh38 or hg38) name must be provided as an argument (--assembly [reference genome assembly name]) or defined in a config file (assembly=[reference genome assembly name]) if FASTQs are being aligned (--fastq_extension [FASTQ file extension] was provided as an argument or fastq_extension=[FASTQ file extension] was defined in a config file."
    #exit 1
#fi

#if [[ -z ${reference_genome} ]]; then
    #echo "Error: a reference genome FASTA must be provided as an argument (--reference_genome [path to reference genome FASTA]) or defined in a config file (reference_genome=[path to reference genome FASTA])."
    #exit 1
#fi

#if [[ ${slurm_mode} == true ]] && [[ -z ${slurm_runtime} ]]; then
    #echo "Error: SLURM runtime must be provided as an argument (--slurm_runtime [SLURM runtime]) or defined in a config file (slurm_runtime=[SLURM runtime]) when running in SLURM mode (--slurm_mode provided as an argument or slurm_mode set as true in config file (slurm_mode=true))."
    #exit 1
#fi

#if [[ ${slurm_mode} == true ]] && [[ -z ${slurm_account} ]]; then
    #echo "Error: SLURM account must be provided as an argument (--slurm_account [SLURM account]) or defined in a config file (slurm_account=[SLURM account]) when running in SLURM mode (--slurm_mode provided as an argument or slurm_mode set as true in config file (slurm_mode=true))."
    #exit 1
#fi

#if [[ ${slurm_mode} == true ]] && [[ -z ${slurm_ncpus} ]]; then
    #echo "Error: number of CPU cores for SLURM job must be provided as an argument (--slurm_ncpus [number of CPU cores for SLURM job] or in a config file (slurm_ncpus=[number of CPU cores for SLURM job]) when running in SLURM mode (--slurm_mode provided as an argument or slurm_mode set as true in config file (slurm_mode=true))."
    #exit 1
#fi

#if [[ ${slurm_mode} == true ]] && [[ -z ${slurm_memory} ]]; then
    #echo "Error: memory usage for SLURM job must be provided as an argument (--slurm_memory [memory usage for SLURM job]) or defined in a config (slurm_memory=[memory usage for SLURM job]) file when running in SLURM mode (--slurm_mode provided as an argument or slurm_mode set as true in config file (slurm_mode=true))."
    #exit 1
#fi

#if [[ ${slurm_mode} == true ]] && [[ -z ${slurm_jobname} ]]; then
    #echo "Warning: a job name for the SLURM job was not provided as an argument (--slurm_jobname [name for SLURM job] or defined in a config file (slurm_jobname=[name for SLURM job]) when running in SLURM mode (--slurm_mode provided as an argument or slurm_mode set as true in config file).  Default job names will be used by SLURM."
#fi

#if [[ ${run_mutect} == false ]] && [[ ${split_intervals} == true ]]; then
    #echo "Error: splitting Mutect analysis across intervals can only be provided as an argument (--split_intervals) or set to true in a config file (split_intervals=true) when Mutect is enabled (--mutect is provided as an argument or run_mutect is set to true in a config file (run_mutect=true))."
    #exit 1
#fi

#if [[ -n ${sequence_dictionary} ]] && ([[ ${run_mutect} == false ]] || [[ ${split_intervals} == false ]]); then
    #echo "Error: a sequence dictionary can only be provided as an argument (--sequence_dictionary [path to sequence dictionary]) or defined in a config file (sequence_dictionary=[path to sequence dictionary]) when Mutect and splitting across intervals are enabled (--mutect and --split_intervals provided as arguments or run_mutect and split_intervals set as true in a config file (run_mutect=true and split_intervals=true))."
    #exit 1
#fi

#if [[ -z ${sequence_dictionary} ]] && [[ ${split_intervals} == true ]]; then
    #echo "Error: a sequence dictionary must be provided as an argument (--sequence_dictionary [path to sequence dictionary]) or defined in a config file (sequence_dictionary=[path to sequence dictionary]) when splitting across intervals is enabled (--split_intervals provided as an argument or split_intervals set to true in a config file (split_intervals=true))."
    #exit 1
#fi

#if [[ ${run_mutect} == false ]] && [[ ${mutect_bam_output} == true ]]; then
    #echo "Error: enabling BAM output from Mutect can only be provided as an argument (--mutect_bam_out) or set to true in a config file (mutect_bam_out=true) when Mutect is enabled (--mutect is provided as an argument or run_mutect is set to true in a config file (run_mutect=true))."
    #exit 1
#fi

#if [[ ${run_mutect} == false ]] && [[ -n ${gnomad_reference} ]]; then
    #echo "Error: a gnomAD reference can only be provided as an argument (--gnomad_reference [path to gnomAD reference]) or defined in a config file (gnomad_reference=[path to gnomAD reference]) when Mutect is enabled (--mutect is provided as an argument or run_mutect is set to true in a config file (run_mutect=true))."
    #exit 1
#fi

#if [[ ${run_mutect} == true ]] && [[ -z ${gnomad_reference} ]]; then
    #echo "Error: a gnomAD reference must be provided as an argument (--gnomad_reference [path to gnomAD reference]) or defined in a config file (gnomad_reference=[path to gnomAD reference]) when Mutect is enabled (--mutect is provided as an argument or run_mutect is set to true in a config file (run_mutect=true))."
    #exit 1
#fi

#if [[ ${run_mutect} == false ]] && [[ -n ${normal_bam} ]]; then
    #echo "Error: a normal BAM file can only be provided as an argument (--normal_bam [path to normal BAM]) or defined in a config file (normal_bam=[path to normal BAM]) when Mutect is enabled (--mutect is provided as an argument or run_mutect is set to true in a config file (run_mutect=true)). Tumor normal calling is currently not implemented for Varscan"
    #exit 1
#fi

#if [[ -z ${normal_pileups_table} ]] && [[ -n ${normal_bam} ]]; then
    #echo "Warning: a normal sample BAM was provided as an argument (--normal_bam [path to normal BAM]) or defined in a config file (normal_bam=[path to normal BAM]) with no corresponding pileups file (by providing --normal_pileups_table [path to normal pileups table] as an argument or defining normal_pileups_path in a config file (normal_pileups_table=[path to normal pileups table])).  Pileups will be calculated before running Mutect."
#fi

#if [[ -n ${normal_pileups_table} ]] && [[ -z ${normal_bam} ]]; then
    #echo "Error: a normal sample pileups table was provided as an argument (--normal_pileups_table [path to normal pileups table] or defined in a config file (normal_pileups_table=[path to normal pileups table]) with no corresponding normal BAM (by providing --normal_bam [path to normal BAM] as an argument or defining normal_bam in a config file (normal_bam=[path to normal BAM]))."
    #exit 1
#fi

#if [[ ${run_funcotator} == true ]] && [[ -z ${funcotator_sources} ]]; then
    #echo "Error: a Funcotator source must be provided as an argument (--funcotator_sources [path to Funcotator source]) or defined in a config file (funcotator_sources=[path to Funcotator source]) when Funcotator is enabled (--skip_funcotator is NOT provided as argument and/or run_funcotator=true is defined in a config file)."
    #exit 1
#fi

#if [[ ${run_funcotator} == false ]] && [[ -n ${funcotator_sources} ]]; then
    #echo "Error: Funcotator must be enabled (--skip_funcotator is NOT provided as argument and/or run_funcotator=true is defined in a config file) when a Funcotator source is provided as an argument (--funcotator_sources [path to Funcotator source]) or defined in a config file (funcotator_sources=[path to Funcotator source])."
    #exit 1
#fi

#if [[ ${run_funcotator} == false ]] && [[ -n ${transcript_list} ]]; then
    #echo "Error: Funcotator must be enabled (--skip_funcotator is NOT provided as argument and/or run_funcotator=true is defined in a config file) when a transcript list is provided as an argument (--transcript_list [path to transcript list]) or defined in a config file (transcript_list=[path to transcript list])."
    #exit 1
#fi

#if [[ ${run_varscan} == false ]] && [[ -n ${mpileup_interval_bed} ]]; then
    #echo "Error: Varscan must be enabled (--varscan is provided as an argument or run_varscan is set to true in a config file (run_varscan=true)) when a BED file with intervals samtools mpileup is provided as an argument (--mpileup_interval_bed [path to interval BED for samtools mpileup]) or defined in a config file (mpileup_interval_bed=[path to interval BED for samtools mpileup])."
    #exit 1
#fi

#if [[ ${run_varscan} == false ]] && [[ -n ${varscan_min_coverage} ]]; then
    #echo "Error: Varscan must be enabled (--varscan is provided as an argument or run_varscan is set to true in a config file (run_varscan=true)) when minimum coverage for Varscan is provided as an argument (--varscan_min_coverage [minimum coverage of variant for Varscan]) or defined in a config file (varscan_min_coverage=[minimum coverage of variant for Varscan])."
    #exit 1
#fi

#if [[ ${run_varscan} == false ]] && [[ -n ${varscan_min_var_freq} ]]; then
    #echo "Error: Varscan must be enabled (--varscan is provided as an argument or run_varscan is set to true in a config file (run_varscan=true)) when minimum variant allele frequency for Varscan is provided as an argument (--varscan_min_var_freq [minimum variant allele fraction for a variant for Varscan]) or defined in a config file (varscan_min_var_freq=[minimum variant allele fraction for a variant for Varscan])."
    #exit 1
#fi

#if [[ ${run_varscan} == false ]] && [[ -n ${varscan_max_pvalue} ]]; then
    #echo "Error: Varscan must be enabled (--varscan is provided as an argument or run_varscan is set to true in a config file (run_varscan=true)) when maximum p-value for Varscan is provided as an argument (--varscan_max_pvalue [maximum p-value for a variant for Varscan]) or defined in a config file (varscan_max_pvalue=[maximum p-value for a variant for Varscan])."
    #exit 1
#fi

#if [[ ${run_annovar} == true ]] && [[ -n ${annovarroot} ]]; then
    #echo "Error: the path to the Annovar root must be provided as an argument (--annovarroot [path to Annovar root]) or defined in a config file (annovarroot [path to Annovar root]) when Annovar is enabled (--skip_annovar is NOT provided as argument and/or run_annovar=true is defined in  a config file)"
    #exit 1
#fi

#if [[ ${run_haplotypecaller} == true ]] && [[ -n ${germline_snps} ]]; then
    #echo "Error: an interval list of germline SNPs must be provided as an argument (--germline_snps [path to interval list of germline SNPs]) or defined in a config file (germline_snps=[path to interval list of germline SNPs]) when HaplotypeCaller is enabled (--haplotypecaller is provided as argument or run_haplotypecaller is set to true in a config file [run_haplotypecaller=true])"
    #exit 1
#fi

#if [[ -n ${fastq_extension} ]] && [[ ${realign} == true ]]; then
    #echo "Error: realignment of a BAM file cannot be requested (--realign is provided an an argument or realign is set to true in a config file [realign=true]) when a FASTQ extension is provided (provided as an argument with --fastq_extension [FASTQ file extension] or defined in a config file (fastq_extension=[FASTQ file extension])), as this implies the starting files are FASTQs, not BAM file."
    #exit 1
#fi
 
##################################################################################################################################
#############################################--STEP 4: GET NORMAL PILEUPS---#######################################################
##################################################################################################################################
    #if [[ $normal_sample != false ]]; then

        #NORMAL_SAMPLE_BASENAME=$(basename "${normal_sample}")
        #NORMAL_SAMPLE_DIRNAME=$(dirname "${normal_sample}")
        #NORMAL_SAMPLE_NAME=${NORMAL_SAMPLE_BASENAME%.*}

        #run_mutect=false

        #if [[ $mode == "slurm" ]]; then
            #sbatch -o "${output_directory}/Logs/%A_%a.log" `#put into log` \
                #-W `#indicates to the script not to move on until the sbatch operation is complete` \
                #"${code_directory}/mutect_and_pileups.sh" \
                #--bam_file ${bam_file} \
                #--intervals_file
                #${parent_directory} \
                #${output_directory} \
                #${min_coverage} \
                #${min_var_freq} \
                #${p_value} \
                #${get_mutect} \
                #${get_varscan} \
                #${get_haplotype} \
                #${file_extension} \
                #${intervals} \
                #${normal_sample} \
                #${code_directory} \
                #${parameter_file} \
                #${reference_genome} \
                #${panel} \
                #${assembly} \
                #${funcotator_sources} \
                #${transcript_list} \
                #${remove_silent} \
                #${run_funcotator} \
                #${bam_out} \
                #${normal_pileups} \
                #${mode} \
                #${docker_image} \
                #${container_engine} \
                #${split_by_chr} \
                #${sequence_dictionary} \
                #${chr_intervals} \
                #${gnomad_genomes} \
                #${run_mutect}
                #wait
        #else
            #echo "CODE_DIRECTORY: $code_directory" >> $parameter_file
            #TASK_ID=1  bash ${code_directory}/BWA_CHIP.sh \
                #${parent_directory} \
                #${output_directory} \
                #${min_coverage} \
                #${min_var_freq} \
                #${p_value} \
                #${get_mutect} \
                #${get_varscan} \
                #${get_haplotype} \
                #${file_extension} \
                #${intervals} \
                #${normal_sample} \
                #${code_directory} \
                #${parameter_file} \
                #${reference_genome} \
                #${panel} \
                #${assembly} \
                #${funcotator_sources} \
                #${transcript_list} \
                #${remove_silent} \
                #${run_funcotator} \
                #${bam_out} \
                #${normal_pileups} \
                #${mode} \
                #${docker_image} \
                #${container_engine} \
                #${split_by_chr} \
                #${sequence_dictionary} \
                #${chr_intervals} \
                #${gnomad_genomes} \
                #${run_mutect}
        #fi

            #normal_pileups="${output_directory}/NORMAL_PILEUPS/${NORMAL_SAMPLE_NAME}_${assembly}_pileups.table"
            #echo "$normal_pileups" >> $parameter_file

            #if [ ! -f ${normal_pileups} ]; then
                #echo "${normal_pileups} does not exist" >> $parameter_file
            #fi
    #fi

##################################################################################################################################
#############################################--STEP 4: GET ARRAY LENGTHS---#######################################################
##################################################################################################################################

code_directory=$(realpath .)
parent_directory=$(dirname "${input_directory}") #get parent directory of $input_directory
fastq_list="${parent_directory}/fastq_files" #give a path to a file to store the paths to the fastq files in $fastq_directory
bam_list="${parent_directory}/bam_files"
mkdir -p "${output_directory}"

if [[ $bam_extension == "bam" ]] || [[ $bam_extension == "cram" ]]; then
    if [[ $bam_extension == "cram" ]] ; then
        type_file="cram"
        type_index_file="crai"
    else
        type_file="bam"
        type_index_file="bai"
    fi

    find -L "${input_directory}" -type f `#list all files in ${input_directory}` | \
                grep ".*\.${type_file}$" `#only keep files with fastq_extension in name (case insensitive)` | \
                grep -v ".*\.${type_index_file}$" `#remove index files` | \
                sed -e 's/\.bam$//g' | \
                sed -e 's/\.cram$//g' | \
                sort -u  `#sort and remove duplicate names` > "${bam_list}"
        bam_array_length=$(wc -l < "${bam_list}") #get the number of FASTQs

    #if [ $realign = true ]; then

        #if [[ $slurm_mode == true ]]; then
            #sbatch -o "${output_directory}/Logs/%A_%a.log" `#put into log` \
                #-a "1-${bam_array_length}" `#initiate job array equal to the number of fastq files` \
                #-W `#indicates to the script not to move on until the sbatch operation is complete` \
                #"${code_directory}/bam_to_fastq.sh" \
                #"$data_directory" \
                #"$code_directory"
        #fi

        #cp $bam_list $fastq_list
        #array_length=$(wc -l < ${fastq_list}) #get the number of FASTQs
        #file_extension="fastq"

    #else
        array_length=$bam_array_length
    #fi
    array_file=${bam_list}
else
    find -L "${input_directory}" -type f `#list all files in ${fastq_directory}` | \
        grep ".*\.${fastq_extension}$" `#only keep files with fastq_extension in name (case insensitive)` | \
        grep -v "Undetermined" `#remove Undetermined FASTQs` | \
        sed -e 's/_R1.*$//g' | sed -e 's/_R2.*$//g' `#remove _R1/2_fastq.gz file extension`| \
        sort -u  `#sort and remove duplicate names` > "${fastq_list}"
    array_length=$(wc -l < "${fastq_list}") #get the number of FASTQs
    array_file=${fastq_list}
fi
    
##################################################################################################################################
##################################################--STEP 4: BWA_CHIP.sh---########################################################
##################################################################################################################################

passed_args_array=(
    array_file 
    output_directory 
    bam_extension 
    fastq_extension 
    assembly 
    reference_genome 
    slurm_mode 
    run_mutect 
    interval_list 
    split_intervals 
    sequence_dictionary 
    mutect_bam_output 
    gnomad_reference 
    normal_bam 
    normal_pileups_table 
    run_funcotator 
    funcotator_sources 
    transcript_list 
    run_varscan 
    mpileup_interval_bed 
    varscan_min_coverage 
    varscan_min_var_freq 
    varscan_max_pvalue 
    run_annovar 
    annovarroot 
    run_haplotypecaller 
    germline_snps
)

passed_args=$(eval echo "$(printf "%s\n" "${passed_args_array[@]}" | xargs --replace=% echo '--% "${%}"' | tr '\n' ' ')")
read -r -a passed_args_array <<< "${passed_args}"

if [[ $slurm_mode == true ]]; then
    mkdir -p "${output_directory}/logs"
    sbatch --output "${output_directory}/logs/%A_%a.log" `#put into log` \
        --error "${output_directory}/logs/%A_%a.log" `#put into log` \
        --array "1-${array_length}" `#initiate job array equal to the number of fastq files` \
        --wait `#indicates to the script not to move on until the sbatch operation is complete` \
        --time "${slurm_runtime}" \
        --account "${slurm_account}" \
        --cpus-per-task "${slurm_ncpus}" \
        --mem "${slurm_memory}" \
        --job-name "${slurm_jobname}" \
        "${code_directory}/BWA_CHIP.sh" "${passed_args_array[@]}"
else
    # change so log output is grouped
    seq 1 "${array_length}" | parallel --progress -j "${n_jobs}" TASK_ID={} "${code_directory}/BWA_CHIP.sh" "${passed_args_array[@]}"
fi
wait
##################################################################################################################################
################################################--STEP 4: RUN R SCRIPTS---########################################################
##################################################################################################################################
#if [[ $slurm_mode == true ]]; then

    #if [[ $run_mutect = true ]]; then
        #mamba run --no-capture-output -n r Rscript aggregate_variants_mutect.R /labs/sjaiswal/chip_submitted_targets_Twist.xls \
            #"$output_directory" "$panel" "$remove_silent" > "$output_directory/Logs/mutectOutFile.Rout" 2>&1

        #mamba run --no-capture-output -n r Rscript WhiteList/whitelist_mutect_join.R "/labs/sjaiswal/variant_whitelist.xlsx" \
            #"$output_directory/mutect_aggregated_simple.tsv" "$output_directory" > "$output_directory/Logs/annotationOutFile.Rout" 2>&1

    #fi

    #if [[ $run_varscan = true ]]; then
        #mamba run --no-capture-output -n r Rscript aggregate_variants_varscan.R /labs/sjaiswal/chip_submitted_targets_Twist.xls \
            #"$output_directory" > "$output_directory/Logs/mutectOutFile.Rout" 2>&1
        #mamba run --no-capture-output -n r Rscript WhiteList/whitelist_varscan_join.R "/labs/sjaiswal/variant_whitelist.xlsx" \
            #"$output_directory/varscan_aggregated.tsv" "$output_directory" > "$output_directory/Logs/annotationOutFile.Rout" 2>&1
    #fi

    #if [[ $run_haplotype = true ]]; then
    #mamba run --no-capture-output -n r Rscript aggregate_variants_haplotypecaller.R /labs/sjaiswal/chip_submitted_targets_Twist.xls \
            #"$output_directory" > "$output_directory/Logs/haplotypeOutFile.Rout" 2>&1
    #fi

    #mutectRout="$output_directory/Logs/mutectOutFile.Rout"
    #if grep -q "Can't combine" "$mutectRout"; then
        #echo "There is an issue with the list containing maf dataframes. The column type varies from data frame to dataframe within list_of_mafs. Check mutectOutFile.Rout to determine which column(s) are causing this error. All columns sharing the same name must share the same type in order for the bind_rows function line 108 to work."
    #fi
#fi
