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

options_array=(
    bam_file
    output_directory
    reference_genome
    pileup_region_intervals
)

longoptions=$(echo "${options_array[@]}" | sed -e 's/ /:,/g' | sed -e 's/$/:/')

# Parse command line arguments with getopt
arguments=$(getopt --options a --longoptions "${longoptions}" --name 'funcotator' -- "$@")
eval set -- "${arguments}"

while true; do
    case "${1}" in
        --bam_file )
            bam_file="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --output_directory )
            output_directory="${2}"; check_for_directory "${1}" "${2}"; shift 2 ;;
        --reference_genome )
            reference_genome="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --pileup_region_intervals )
            pileup_region_intervals="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        -- )
            shift; break;;
        * )
            echo "Invalid argument ${1} ${2}" >&2
            exit 1
    esac
done

sample_name=$(basename "${bam_file}" | sed -e 's/.bam$//g'  | sed -e 's/.cram$//g')
pileup_name="${output_directory}/${sample_name}"

if [[ ! -f "${pileup_name}" ]]; then
    mamba run -n pileup_region pileup_region \
        "${pileup_region_intervals}" \
        "${bam_file}" \
        "${reference_genome}" > "${pileup_name}.pileup_region"

    pileup_region_nrows=$(head -n -1 "${pileup_name}.pileup_region" | wc -l)
    paste <(yes "${sample_name}" | head -n "${pileup_region_nrows}") <(head -n -1 "${pileup_name}.pileup_region") | sponge "${pileup_name}.pileup_region"
else
    echo "Region pileups already computed in: ${pileup_name}.pileup_region"
fi
