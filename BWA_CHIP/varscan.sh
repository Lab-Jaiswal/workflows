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
    bam_file
    reference_genome
    varscan_min_coverage
    varscan_min_var_freq
    varscan_max_pvalue
    annovarroot
    mpileup_interval_bed
    run_annovar
)

longoptions=$(echo "${options_array[@]}" | sed -e 's/ /:,/g' | sed -e 's/$/:/')

# Parse command line arguments with getopt
arguments=$(getopt --options a --longoptions "${longoptions}" --name 'varscan' -- "$@")
eval set -- "${arguments}"

while true; do
    case "${1}" in
        --bam_file )
            bam_file="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --reference_genome )
            reference_genome="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --varscan_min_coverage )
            varscan_min_coverage="${2}"; shift 2 ;;
        --varscan_min_var_freq )
            varscan_min_var_freq="${2}"; shift 2 ;;
        --varscan_max_pvalue )
            varscan_max_pvalue="${2}"; shift 2 ;;
        --run_annovar )
            run_annovar="${2}"; shift 2 ;;
        --annovarroot )
            annovarroot="${2}"; check_for_directory "${1}" "${2}"; shift 2 ;;
        --mpileup_interval_bed )
            mpileup_interval_bed="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        -- )
            shift; break ;;
        * )
            echo "Invalid argument ${1} ${2}" >&2
            exit 1
    esac
done

sample_name=$(echo "${bam_file}" | sed -e 's/\.bam//g')
declare optional_args=""

if [[ ${mpileup_interval_bed} != "none" ]]; then
    optional_args="--positions ${mpileup_interval_bed}"
fi
read -r -a optional_args_array <<< "${optional_args}"

if [[ ! -f "${sample_name}.pileup" ]]; then
    echo "Generating pileup from BAM..."
    mamba run -n samtools samtools mpileup \
        --count-orphans \
        --max-depth 0 \
        --adjust-MQ 50 \
        --fasta-ref "${reference_genome}" \
        "${optional_args_array[@]}" \
        "${bam_file}" > "${sample_name}.pileup"
    echo "...pileup generated"
else
    echo "Pileup already generated"
fi

if [[ ! -f "${sample_name}_varscan2.vcf" ]]; then
    echo "Calling variants from pileup..."
    mamba run -n varscan varscan mpileup2cns \
    "${sample_name}.pileup" \
        -Xmx64g \
        --min-coverage "${varscan_min_coverage}" \
        --min-var-freq "${varscan_min_var_freq}" \
        --p-value "${varscan_max_pvalue}" \
        --output-vcf 1 > "${sample_name}_varscan2.vcf"
    echo "...variants called"
else
    echo "Variants already called"
fi

if [[ ! -f "${sample_name}_varscan2_filter.vcf" ]]; then
    echo "Filtering variants in VCF..."
    mamba run -n varscan varscan filter \
        "${sample_name}_varscan2.vcf" \
        --output-file "${sample_name}_varscan2_filter.vcf" \
        --min-coverage "${varscan_min_coverage}" \
        --min-var-freq "${varscan_min_var_freq}" \
        --p-value "${varscan_max_pvalue}"
    echo "...variants filtered"
else
    echo "Variants already filtered"
fi

annovar_vcf="${sample_name}_varscan2_filter_annovar.hg38_multianno.vcf"
if [[ ! -f "${annovar_vcf}" ]] && [[ $run_annovar == true ]]; then
    echo "Annotating VCF with Annovar..."
    perl "${annovarroot}/table_annovar.pl" \
        "${sample_name}_varscan2_filter.vcf" \
        "${annovarroot}/humandb" \
        --buildver hg38 \
        --remove \
        --outfile "${sample_name}_varscan2_filter_annovar" \
        --protocol ensGene \
        --operation g \
        --nastring '.' \
        --vcfinput \
        --thread 1
    echo "...VCF annotated"

    bgzip --stdout "${annovar_vcf}" > "${annovar_vcf}.gz"
    tabix --preset vcf --force "${annovar_vcf}.gz"

    fixed_aachange_columns="${annovar_vcf//.vcf/_fixed_aachange_columns}"
    paste <(bcftools query --format "%CHROM\t%POS\t%REF\t%ALT\n" "${annovar_vcf}.gz") \
        <(bcftools query --format "%AAChange.ensGene\n" "${annovar_vcf}.gz" | tr ':' '|') | \
        bgzip --stdout > "${fixed_aachange_columns}"
    tabix --force --sequence 1 --begin 2 --end 2 "${fixed_aachange_columns}"

    aachange_columns="Transcript.ensGene|Exon.ensGene|cDNAChange.ensGene|proteinChange.ensGene"
    bcftools view "${annovar_vcf}.gz" | \
        sed "s/AAChange.ensGene annotation provided by ANNOVAR/Format: Consequence|${aachange_columns}/" | \
        bcftools annotate --annotations "${fixed_aachange_columns}" --columns "CHROM,POS,REF,ALT,AAChange.ensGene" | \
        bgzip | sponge "${annovar_vcf}.gz"

    bcftools +split-vep --annotation "AAChange.ensGene" --columns "${aachange_columns//|/,}" --annot-prefix "" ${annovar_vcf}.gz | \
        bcftools norm --multiallelics - | \
        bgzip | sponge "${annovar_vcf}.gz"
else
    echo "VCF already annotated"
fi
