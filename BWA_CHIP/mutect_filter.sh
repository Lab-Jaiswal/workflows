#!/usr/bin/env bash

# Run scripts to enable activating conda environments
. "${HOME}/micromamba/etc/profile.d/conda.sh"
. "${HOME}/micromamba/etc/profile.d/mamba.sh"
mamba activate base

set -o xtrace -o nounset -o pipefail -o errexit

# Parse command line arguments with getopt
arguments=$(getopt --options a --longoptions vcf_file:,pileups_table:,normal_pileups_table:,orientation_bias_priors:,reference_genome:,gatk_command: --name 'mutect_filter' -- "$@")
eval set -- "${arguments}"

check_for_file() {
    argument_name="${1}"
    file_path="${2}"
    if [[ ${file_path} != false ]] && [[ ! -f ${file_path} ]]; then
        echo "Error: file ${file_path} passed with ${argument_name} does not exist."
        exit 1
    fi
}

while true; do
    case "${1}" in 
        --vcf_file )
            vcf_file="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --pileups_table )
            pileups_table="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --normal_pileups_table )
            normal_pileups_table="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --orientation_bias_priors )
            orientation_bias_priors="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --reference_genome )
            reference_genome="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --gatk_command )
            gatk_command="${2}"; shift 2 ;;
        --)
            shift; break;;
    esac
done

optional_args=""
sample_name=$(echo "${vcf_file}" | sed -e 's/_mutect2.vcf//g' )

# Calculate germline contamination using pileups at known common germline variants.
if [[ ! -f "${sample_name}_contamination.table" ]] ; then
    if [[ "${normal_pileups_table}" != false ]]; then
        optional_args="--matched-normal ${normal_pileups_table}"
    fi

    echo "Getting contamination rate with CalculateContamination..."
    ${gatk_command} CalculateContamination \
        --input "${pileups_table}" \
        ${optional_args} \
        --output "${sample_name}_contamination.table"
    echo "...contamination rate calculated."
else
    echo "Contamination rate already calculated in: ${sample_name}_contamination.table "
fi

# Add FILTER column to Mutect2 VCF to identify variants which pass or fail filters.
if [[ ! -f "${sample_name}_mutect2_filtered.vcf" ]]; then
    echo "Filtering somatic variants with FilterMutectCalls..."
    ${gatk_command} FilterMutectCalls \
        --variant "${vcf_file}" \
        --output "${sample_name}_mutect2_filtered.vcf" \
        --contamination-table "${sample_name}_contamination.table" \
        --ob-priors  "${orientation_bias_priors}" \
        --reference "${reference_genome}"
    echo "...somatic variants filtered."
else
    echo "Mutect2 somatic variants already filtered"
fi
