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
    if [[ ! -f ${file_path} ]]; then
        echo "Error: file ${file_path} passed with ${argument_name} does not exist."
        exit 1
    fi
}

# Parse command line arguments with getopt
arguments=$(getopt --options a --longoptions bam_file:,reference_genome:,germline_snps:,gatk_command: --name 'haplotypecaller' -- "$@")
eval set -- "${arguments}"

while true; do
    case "${1}" in
        --bam_file )
            bam_file="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --reference_genome )
            reference_genome="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --germline_snps )
            germline_snps="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --gatk_command )
            gatk_command="${2}"; shift 2 ;;
        -- )
            shift; break ;;
        * )
            echo "Invalid argument ${1} ${2}" >&2
            exit 1
    esac
done

sample_name=$(echo "${bam_file}" | sed -e 's/\.bam$//g' | sed -e 's/\.cram$//g')

if [[ ! -f "${sample_name}_haplotypecaller.gvcf" ]]; then
    echo "Calling germline variants with HaplotypeCaller..."
    ${gatk_command} HaplotypeCaller \
        --input "${bam_file}" \
        --output "${sample_name}_haplotypecaller.gvcf" \
        --reference "${reference_genome}" \
        --intervals "${germline_snps}" \
        --bamout "${sample_name}_haplotypecaller.bam" \
        --emit-ref-confidence GVCF

    mamba run -n samtools samtools index "${sample_name}_haplotypecaller.bam" "${sample_name}_haplotypecaller.bam.bai"

    echo "...germline variants called."
else
    echo "HaplotypeCaller gVCF already exists"
fi

if [ ! -f "${sample_name}_haplotypecaller.vcf" ]; then
    echo "Genotyping germline variants in gVCF with GenotypeGVCF..."
    ${gatk_command} GenotypeGVCFs \
        --variant "${sample_name}_haplotypecaller.gvcf" \
        --reference "${reference_genome}" \
        --force-output-intervals "${germline_snps}" \
        --output "${sample_name}_haplotypecaller.vcf"
    echo "...germline variants genotyped."
else
    echo "Germline variants already genotyped with GenotypeGVCFs"
fi
