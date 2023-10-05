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
#############################################---STEP 1: SET UP PARAMETERS---###################################################### 
##################################################################################################################################
options_array=(
    bam_file
    bam_extension
    output_directory
    code_directory
    reference_genome
    interval_list_dir
    exac_reference_dir
)

longoptions=$(echo "${options_array[@]}" | sed -e 's/ /:,/g'):

# Parse command line arguments with getopt
arguments=$(getopt --options a --longoptions "${longoptions}" --name 'scatter_mutect_and_pileups' -- "$@")
eval set -- "${arguments}"

while true; do
    case "${1}" in
        --bam_file )
            bam_file="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --bam_extension )
            bam_extension="${2}"; shift 2 ;;
        --output_directory )
            final_output_directory="${2}"; check_for_directory "${1}" "${2}"; shift 2 ;;
        --code_directory )
            code_directory="${2}"; check_for_directory "${1}" "${2}"; shift 2 ;;
        --reference_genome )
            reference_genome="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --interval_list_dir )
            interval_list_dir="${2}"; check_for_directory "${1}" "${2}"; shift 2 ;;
        --exac_reference_dir )
            exac_reference_dir="${2}"; check_for_dir "${1}" "${2}"; shift 2 ;;
        --split_jobs )
            split_jobs="${2}"; shift 2 ;;
        -- )
            shift; break ;;
        * )
            echo "Invalid argument ${1} ${2}" >&2
            exit 1
    esac
done

sample_name=$(basename "${bam_file//.${bam_extension}/}")
gatk_command="mamba run -n gatk4 gatk"
output_directory=${output_directory}/${sample_name}

num_intervals=$(ls "${interval_list_dir}" | wc -l)
echo "Number of intervals: $num_intervals"
seq 1 "${num_intervals}" | parallel -j${split_jobs} --progress --ungroup \
    "${code_directory}/mutect_and_pileups.sh" \
        --bam_file "${bam_file}" \
        --interval_list_dir "${interval_list_dir}" \
        --interval_number {} \
        --reference_genome "${reference_genome}" \
        --exac_reference_dir "${exac_reference_dir}" \
        --output_directory "${output_directory}" \
        --gatk_command "${gatk_command}"
