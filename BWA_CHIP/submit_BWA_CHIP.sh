#!/usr/bin/env bash

# Run scripts to enable activating conda environments
. "${HOME}/micromamba/etc/profile.d/conda.sh"
. "${HOME}/micromamba/etc/profile.d/mamba.sh"
mamba activate base

# Set bash options for verbose output and to fail immediately on errors or if variables are undefined.
#set -o xtrace -o nounset -o pipefail -o errexit
set -o xtrace -o nounset -o errexit

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
    slurm_partition
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
    annovarroot
    pileup_region_intervals
    germline_snps
    n_jobs
    slurm_mode
    run_mutect
    run_varscan
    run_pileup_region
    run_haplotypecaller
    mutect_bam_output
    varscan_min_coverage
    varscan_min_var_freq
    varscan_max_pvalue
    run_annovar
    run_funcotator
    split_intervals
    realign
    whitelist
    filter_silent
    remove_artifacts
)

eval "$(printf "%s\n" "${options_array[@]}" | xargs --replace=% echo "declare %=none;")"
longoptions="$(echo "${options_array[@]}" | sed -e 's/ /:,/g'):"

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
        --config_file )
            shift 2 ;;
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
            slurm_mode="$2"; shift 2 ;;
        --slurm_runtime )
            slurm_runtime="$2"; shift 2 ;;
        --slurm_account )
            slurm_account="$2"; shift 2 ;;
        --slurm_partition )
            slurm_partition="$2"; shift 2 ;;
        --slurm_ncpus )
            slurm_ncpus="$2"; shift 2 ;;
        --slurm_memory )
            slurm_memory="$2"; shift 2 ;;
        --slurm_jobname )
            slurm_jobname="$2"; shift 2 ;;
        --run_mutect )
            run_mutect="$2"; shift 2 ;;
        --interval_list )
            interval_list="$2"; check_for_file "${1}" "${2}"; shift 2 ;;
        --split_intervals )
            split_intervals="$2"; shift 2 ;;
        --sequence_dictionary )
            sequence_dictionary="$2"; check_for_file "${1}" "${2}"; shift 2 ;;
        --mutect_bam_output )
            mutect_bam_output="$2"; shift 2 ;;
        --gnomad_reference )
            gnomad_reference="$2"; check_for_file "${1}" "${2}"; shift 2 ;;
        --normal_bam )
            normal_bam="$2"; check_for_file "${1}" "${2}"; shift 2 ;;
        --normal_pileups_table )
            normal_pileups_table="$2"; check_for_file "${1}" "${2}"; shift 2 ;;
        --run_funcotator )
            run_funcotator="$2"; shift 2 ;;
        --funcotator_sources )
            funcotator_sources="$2"; check_for_directory "${1}" "${2}"; shift 2 ;;
        --transcript_list )
            transcript_list="$2"; check_for_file "${1}" "${2}"; shift 2 ;;
        --run_varscan )
            run_varscan="$2"; shift 2 ;;
        --mpileup_interval_bed )
            mpileup_interval_bed="$2"; check_for_file "${1}" "${2}"; shift 2 ;;
        --varscan_min_coverage )
            varscan_min_coverage="$2"; shift 2 ;;
        --varscan_min_var_freq )
            varscan_min_var_freq="$2"; shift 2 ;;
        --varscan_max_pvalue )
            varscan_max_pvalue="$2"; shift 2 ;;
        --run_annovar )
            run_annovar="$2"; shift 2 ;;
        --annovarroot )
            annovarroot="$2"; check_for_directory "${1}" "${2}"; shift 2 ;;
        --run_pileup_region )
            run_pileup_region="${2}"; shift 2 ;;
        --pileup_region_intervals )
            pileup_region_intervals="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --run_haplotypecaller )
            run_haplotypecaller="$2"; shift 2 ;;
        --germline_snps )
            germline_snps="$2"; check_for_file "${1}" "${2}"; shift 2 ;;
        --n_jobs )
            n_jobs="$2"; shift 2 ;;
        --realign )
            realign="$2"; shift 2 ;;
        --whitelist )
            whitelist="$2"; check_for_file "${1}" "${2}"; shift 2 ;;
        --filter_silent )
            filter_silent="$2"; shift 2 ;;
        --remove_artifacts )
            remove_artifacts="$2"; shift 2 ;;
        -- )
            shift; break ;;
        * )
            shift; echo "Invalid argument ${1}" >&2; exit 1
    esac
done

##################################################################################################################################
#############################################--STEP 4: GET ARRAY LENGTHS---#######################################################
##################################################################################################################################

code_directory=$(realpath $(dirname ${BASH_SOURCE[0]}))
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
    code_directory
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
    run_pileup_region
    pileup_region_intervals
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
        --partition "${slurm_partition}" \
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
if [[ $slurm_mode == true ]]; then
    if [[ $run_mutect = true ]]; then
        ${code_directory}/filter_variants_targeted.sh \
            --vcf_directory=${output_directory} \
            --output_directory=${output_directory} \
            --filter_silent=${filter_silent} 

        mamba run --no-capture-output -n r Rscript ${code_directory}/annotate_chip_calls.R \
            --output_directory "${output_directory}" \
            --chip_calls "${output_directory}/aggregated.tsv" \
            --pileup_regions "${output_directory}/aggregated.pileup_region" \
            --whitelist "${whitelist}" \
            --remove_artifacts "${remove_artifacts}"
    fi

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
fi
