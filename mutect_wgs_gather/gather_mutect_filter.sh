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

options_array=(
  input_directory
  output_directory
  min_sequencing_depth
  max_sequencing_depth
  bravo_variants
  sequence_context_window_size
  reference_genome
  sequence_dictionary
  n_jobs
)

longoptions=$(echo "${options_array[@]}" | sed -e 's/ /:,/g' | sed -e 's/$/:/')

# Parse command line arguments with getopt
arguments=$(getopt --options a --longoptions "${longoptions}" --name 'gather_mutect_filter' -- "$@")
eval set -- "${arguments}"

while true; do
    case "${1}" in
        --input_directory )
            input_directory="${2}"; shift 2 ;;
        --output_directory )
            output_directory="${2}"; shift 2 ;;
        --n_jobs )
            n_jobs="${2}"; shift 2;;
        --min_sequencing_depth )
            min_sequencing_depth="${2}"; shift 2 ;;
        --max_sequencing_depth )
            max_sequencing_depth="${2}"; shift 2 ;;
        --bravo_variants )
            bravo_variants="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --sequence_context_window_size )
            sequence_context_window_size="${2}"; shift 2 ;;
        --reference_genome )
            reference_genome="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --sequence_dictionary )
            sequence_dictionary="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        -- )
            shift; break;;
        * )
            echo "Invalid argument ${1} ${2}" >&2
            exit 1
    esac
done

vcf_list=${input_directory}/vcf_list
find -L "${input_directory}" -maxdepth 1 -mindepth 1 -type d  `#list all files in ${vcf_directory}` | \
    sort -u  `#sort and remove duplicate names` > ${vcf_list}
vcf_array_length=$(wc -l < ${vcf_list}) #get the number of FASTQs

mkdir -p "${output_directory}/logs"
code_directory=$(realpath $(dirname ${BASH_SOURCE[0]}))
seq 1 ${vcf_array_length} | parallel --eta -j${n_jobs} TASK_ID={} "${code_directory}/mutect_filter.sh" \
    --array_file "${vcf_list}" \
    --output_directory ${output_directory} \
    --min_sequencing_depth ${min_sequencing_depth} \
    --max_sequencing_depth ${max_sequencing_depth} \
    --bravo_variants ${bravo_variants} \
    --sequence_context_window_size ${sequence_context_window_size} \
    --reference_genome ${reference_genome} \
    --sequence_dictionary ${sequence_dictionary}
