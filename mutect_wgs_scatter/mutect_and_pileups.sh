#!/usr/bin/env bash

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
    bam_file
    bam_extension
    interval_list_directory
    interval_number
    reference_genome
    exac_reference_directory
    output_directory
)

longoptions=$(echo "${options_array[@]}" | sed -e 's/ /:,/g' | sed -e 's/$/:/')

# Parse command line arguments with getopt
arguments=$(getopt --options a --longoptions "${longoptions}" --name 'mutect_and_pileups' -- "$@")
eval set -- "${arguments}"

while true; do
    case "${1}" in
        --bam_file )
            bam_file="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --bam_extension )
            bam_extension="${2}"; shift 2 ;;
        --interval_list_directory )
            interval_list_directory="${2}"; check_for_directory "${1}" "${2}"; shift 2 ;;
        --interval_number )
            interval_number="${2}"; shift 2 ;;
        --reference_genome )
            reference_genome="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --exac_reference_directory )
            exac_reference_directory="${2}"; check_for_directory "${1}" "${2}"; shift 2 ;;
        --output_directory )
            output_directory="${2}"; shift 2 ;;
        --)
            shift; break;;
        * )
            echo "Invalid argument ${1} ${2}" >&2
            exit 1
    esac
done

pushd ${interval_list_directory}
    interval_list_basename=$(ls *.interval_list | sort -h | sed "${interval_number}q; d")
    interval_list="${interval_list_directory}/${interval_list_basename}"
popd
pushd ${exac_reference_directory}
    exac_reference_basename=$(ls *.vcf.gz | sort -h | sed "${interval_number}q; d" ) 
    exac_reference="${exac_reference_directory}/${exac_reference_basename}"
popd
chromosome_name="chr$(echo ${interval_list_basename//.interval_list/} | sed "s/_.*$//g")"
sample_name="$(basename ${bam_file//.${bam_extension}/})"
sample_name="${sample_name}_${chromosome_name}"

mkdir -p "${output_directory}/vcfs"
mkdir -p "${output_directory}/pileups"
mkdir -p "${output_directory}/f1r2"

vcf_name="${output_directory}/vcfs/${sample_name}"
pileup_name="${output_directory}/pileups/${sample_name}"
f1r2_name="${output_directory}/f1r2/${sample_name}"

# Get pileup summaries at known common germline variants to estimate germline contamination.
echo "Getting pileup summaries with GetPileupSummaries..."
gatk GetPileupSummaries \
    --input "${bam_file}" \
    --variant "${exac_reference}" \
    --intervals "${exac_reference}" \
    --reference "${reference_genome}" \
    --output "${pileup_name}_pileups.table"
echo "...pileup summaries calculated."

# Call somatic variants with Mutect2 and output orientation bias read counts
echo "Calling somatic variants with Mutect2..."
gatk Mutect2 \
    --input "${bam_file}" \
    --output "${vcf_name}_mutect2.vcf" \
    --reference "${reference_genome}" \
    --intervals "${interval_list}" \
    --dont-use-soft-clipped-bases \
    --f1r2-tar-gz "${f1r2_name}_f1r2.tar.gz" \
    --annotation OrientationBiasReadCounts
echo "...somatic variants called."
