#!/usr/bin/env bash

# Run scripts to enable activating conda environments
. "${HOME}/micromamba/etc/profile.d/conda.sh"
. "${HOME}/micromamba/etc/profile.d/mamba.sh"
mamba activate base
mamba activate bcftools

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
    filtered_vcf
    reference_genome
    funcotator_sources
    transcript_list
    gatk_command
)

longoptions=$(echo "${options_array[@]}" | sed -e 's/ /:,/g' | sed -e 's/$/:/')

# Parse command line arguments with getopt
arguments=$(getopt --options a --longoptions "${longoptions}" --name 'funcotator' -- "$@")
eval set -- "${arguments}"

while true; do
    case "${1}" in
        --filtered_vcf )
            filtered_vcf="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --reference_genome )
            reference_genome="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --funcotator_sources )
            funcotator_sources="${2}"; check_for_directory "${1}" "${2}"; shift 2 ;;
        --transcript_list )
            transcript_list="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --gatk_command )
            gatk_command="${2}"; shift 2 ;;
        -- )
            shift; break;;
        * )
            echo "Invalid argument ${1} ${2}" >&2
            exit 1
    esac
done

declare optional_args=""
sample_name=$(echo "${filtered_vcf}" | sed -e 's/_mutect2_filtered.vcf//g')

if [[ ${transcript_list} != "none" ]]; then
    optional_args="--transcript-list ${transcript_list}"
fi
read -r -a optional_args_array <<< "${optional_args}"

# Remove non-standard chromosomes from VCF header
reheadered_vcf="${sample_name}_mutect2_filtered_reheadered.vcf"
if [[ ! -f "${reheadered_vcf}" ]]; then
    new_header="${sample_name}_new_header"

    grep "^#" < "${filtered_vcf}" | grep -v -E "chrM|chr.*_|HLA|chrEBV" > "${new_header}"

    echo "Removing non-standard chromosomes from VCF header..."
    ${gatk_command} FixVcfHeader \
        --INPUT "${filtered_vcf}" \
        --OUTPUT "${reheadered_vcf}" \
        --HEADER "${new_header}"
    echo "...header fixed."
fi

funcotator_vcf="${sample_name}_mutect2_filtered_funcotator.vcf"
# Annotate filtered Mutect2 VCF with Funcotator annotations
if [[ ! -f "${funcotator_vcf}" ]]; then
    echo "Annotating Mutect2 VCF with Funcotator..."
    ${gatk_command} Funcotator \
        --variant "${reheadered_vcf}" \
        --reference "${reference_genome}" \
        --ref-version hg38 \
        --data-sources-path "${funcotator_sources}" \
        "${optional_args_array[@]}" \
        --output "${funcotator_vcf}" \
        --output-file-format VCF 
    echo "...VCF annotated."
else
    echo "Mutect2 VCF already annotated"
fi

funcotator_header=$(bcftools view --header-only "${funcotator_vcf}" | grep FUNCOTATION | sed 's/.*Description="//g' | sed 's/">//g')
funcotator_replacement_header=$(echo "${funcotator_header}" | \
    sed "s/.*:/Format:/g" | \
    sed "s/Gencode_[0-9]*_variantClassification/Consequence/g" | \
    tr '/' '_' | \
    tr '(' '_' | \
    tr ')' '_')
funcotator_columns=$(echo "${funcotator_replacement_header}" | sed 's/.*: //g' | tr '|' ',')

bcftools view "${funcotator_vcf}" | \
    sed "s/FUNCOTATION,Number=A,Type=String,Description=\".*\">/FUNCOTATION,Number=A,Type=String,Description=\"${funcotator_replacement_header}\">/g" | \
    tr --delete '[' | \
    tr --delete ']' | \
    bgzip >! "${funcotator_vcf}.gz"
tabix --preset vcf --force "${funcotator_vcf}.gz"

bcftools +split-vep --annotation "FUNCOTATION" --columns "${funcotator_columns}" --annot-prefix "funcotator_" "${funcotator_vcf}.gz" | bgzip | sponge "${funcotator_vcf}.gz"
tabix --preset vcf --force "${funcotator_vcf}.gz"
