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
    if [[ -n ${file_path} ]] && [[ ! -f ${file_path} ]]; then
        echo "Error: file ${file_path} passed with ${argument_name} does not exist."
        exit 1
    fi
}

check_for_directory() {
    argument_name="${1}"
    directory_path="${2}"
    if [[ -n ${directory_path} ]] && [[ ! -d ${directory_path} ]]; then
        echo "Error: directory ${directory_path} passed with ${argument_name} does not exist."
        exit 1
    fi
}

options_array=(
    bam_file
    normal_bam_file
    interval_list
    interval_number
    reference_genome
    gnomad_reference
    output_directory
    mutect_bam_output
    run_mutect
    gatk_command
)

longoptions=$(echo "${options_array[@]}" | sed -e 's/ /:,/g' | sed -e 's/$/:/')

# Parse command line arguments with getopt
arguments=$(getopt --options a --longoptions "${longoptions}" --name 'mutect_and_pileups' -- "$@")
eval set -- "${arguments}"

# Set defaults for some variables
declare interval_number

while true; do
    case "${1}" in
        --bam_file )
            bam_file="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --normal_bam_file )
            normal_bam_file="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --interval_list )
            interval_list="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --interval_number )
            interval_number="${2}"; shift 2 ;;
        --reference_genome )
            reference_genome="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --gnomad_reference )
            gnomad_reference="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --output_directory )
            output_directory="${2}"; check_for_directory "${1}" "${2}"; shift 2 ;;
        --mutect_bam_output )
            mutect_bam_output="${2}"; shift 2 ;;
        --gatk_command )
            gatk_command="${2}"; shift 2 ;;
        --run_mutect )
            run_mutect="${2}"; shift 2 ;;
        --)
            shift; break;;
        * )
            echo "Invalid argument ${1} ${2}" >&2
            exit 1
    esac
done

sample_name=$(basename "${bam_file}" | sed -e 's/.bam$//g'  | sed -e 's/.cram$//g')
declare optional_args # Initially blank so that optional arguments can be filled in later

# Extract header from normal BAM/CRAM file and get sample name
if [[ -n ${normal_bam_file} ]]; then
    normal_name=$(mamba run --no-capture-output -n samtools samtools samples -h "${normal_bam_file}" | tail -n 1 | cut -f 1)
    optional_args="${optional_args} --input ${normal_bam_file} --normal ${normal_name}"
fi
         
# If we are splitting by chromosome to run Mutect2 on the whole genome, subset the interval list file to only one interval.
if [[ -n ${interval_number} ]]; then
    interval_line=$(grep -v "@" < "${interval_list}" | sed "${interval_number}q; d" ) # Remove header from interval list and choose appropriate line
    interval_name=$(echo "${interval_line}" | cut -f 1) # Extract chromosome name from interval line
    new_interval_list="${interval_name}.interval_list"
    grep "@" < "${interval_list}" > "${new_interval_list}" # Copy header from original interval list into new interval list
    echo "${interval_line}" >> "${new_interval_list}" # Append line after header in new interval list
    sample_name="${interval_name}_${sample_name}"
    interval_list="${new_interval_list}"

    mkdir -p "${output_directory}/vcfs"
    mkdir -p "${output_directory}/pileups"
    mkdir -p "${output_directory}/f1r2"

    vcf_name="${output_directory}/vcfs/${sample_name}"
    pileup_name="${output_directory}/pileups/${sample_name}"
    f1r2_name="${output_directory}/f1r2/${sample_name}"
else
    vcf_name="${output_directory}/${sample_name}"
    pileup_name="${output_directory}/${sample_name}"
    f1r2_name="${output_directory}/${sample_name}"
fi

if [[ ${mutect_bam_output} == true ]]; then
    optional_args="${optional_args} --bamout ${vcf_name}_mutect2.bam"
else
    echo "Mutect2 BAM output not requested"
fi

# Call somatic variants with Mutect2 and output orientation bias read counts and optionally a BAM file with additional annotations.
if [[ ! -f "${vcf_name}_mutect2.vcf" ]] && [[ ${run_mutect} == true ]]; then
    read -r -a optional_args_array <<< "${optional_args}"
    echo "Calling somatic variants with Mutect2..."
    ${gatk_command} Mutect2 \
        --input "${bam_file}" \
        "${optional_args_array[@]}" \
        --output "${vcf_name}_mutect2.vcf" \
        --reference "${reference_genome}" \
        --intervals "${interval_list}" \
        --dont-use-soft-clipped-bases \
        --f1r2-tar-gz "${f1r2_name}_f1r2.tar.gz" \
        --annotation OrientationBiasReadCounts

    if [[ ${mutect_bam_output} = true ]]; then
        ${gatk_command} BuildBamIndex \
        --INPUT "${vcf_name}_mutect2.bam"
    fi

    echo "...somatic variants called."
else
    echo "Somatic variants already called with Mutect2 in: ${vcf_name}_mutect2.vcf"
fi

# Get pileup summaries at known common germline variants to estimate germline contamination.
if  [[ ! -f "${pileup_name}_pileups.table" ]]; then
    echo "Getting pileup summaries with GetPileupSummaries..."
    ${gatk_command} GetPileupSummaries \
        --input "${bam_file}" \
        --variant "${gnomad_reference}" \
        --intervals "${interval_list}" \
        --reference "${reference_genome}" \
        --output "${pileup_name}_pileups.table"
    echo "...pileup summaries calculated."
else
    echo "Pileup summaries already calculated in: ${pileup_name}_pileups.table"
fi
