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

##################################################################################################################################
#############################################---STEP 1: SET UP PARAMETERS---###################################################### 
##################################################################################################################################
options_array=(
    array_file
    output_directory
    code_directory
    reference_genome
    sequence_dictionary
    bravo_variants
    min_sequencing_depth
    max_sequencing_depth
    sequence_context_window_size
)

longoptions=$(echo "${options_array[@]}" | sed -e 's/ /:,/g'):

# Parse command line arguments with getopt
arguments=$(getopt --options a --longoptions "${longoptions}" --name 'somatic_variant_pipeline' -- "$@")
eval set -- "${arguments}"

while true; do
    case "${1}" in
        --array_file )
            array_file="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --output_directory )
            output_directory="${2}"; shift 2 ;;
        --code_directory )
            code_directory="${2}"; check_for_directory "${1}" "${2}"; shift 2 ;;
        --reference_genome )
            reference_genome="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --sequence_dictionary )
            sequence_dictionary="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --bravo_variants )
            bravo_variants="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --min_sequencing_depth )
            min_sequencing_depth="${2}"; shift 2 ;;
        --max_sequencing_depth )
            max_sequencing_depth="${2}"; shift 2 ;;
        --sequence_context_window_size )
            sequence_context_window_size="${2}"; shift 2 ;;
        -- )
            shift; break ;;
        * )
            echo "Invalid argument ${1} ${2}" >&2
            exit 1
    esac
done

line_number=${TASK_ID}
input_directory="$(sed "${line_number}q; d" "${array_file}")" #extract only the line number corresponding to $SLURM_ARRAY_TASK_ID
sample_name=$(basename "${input_directory}")
gatk_command="mamba run -n gatk4 gatk"

output_directory="${output_directory}/${sample_name}"
mkdir -p ${output_directory}

# Gather pileup summaries
if [[ ! -f "${output_directory}/${sample_name}_pileups.table" ]]; then
    pileup_tables=$(find "${input_directory}/pileups" -type f | grep -E ".*_pileups.table$" | sort -V | sed -e 's/^/--I /g' | tr '\n' ' ')
    read -r -a pileup_tables_array <<< "${pileup_tables}"
    mamba run -n gatk4 gatk GatherPileupSummaries \
        --sequence-dictionary "${sequence_dictionary}" \
        "${pileup_tables_array[@]}" \
        -O "${output_directory}/${sample_name}_pileups.table"
fi

# Gather VCFs
if [[ ! -f "${output_directory}/${sample_name}_mutect2.vcf" ]]; then
    vcf_files=$(find "${input_directory}/vcfs" -type f | grep -E ".*.vcf$" | sort -V | sed -e 's/^/--INPUT /g' | tr '\n' ' ')
    read -r -a vcf_files_array <<< "${vcf_files}"
    mamba run -n gatk4 gatk MergeVcfs \
        "${vcf_files_array[@]}" \
        --OUTPUT "${output_directory}/${sample_name}_mutect2.vcf"
fi

# Gather VCF stats
if [[ ! -f "${output_directory}/${sample_name}_mutect2.vcf.stats" ]]; then
    vcf_stats=$(find "${input_directory}/vcfs" -type f | grep -E ".*.vcf.stats$" | sort -V | sed -e 's/^/--stats /g' | tr '\n' ' ')
    read -r -a vcf_stats_array <<< "${vcf_stats}"
    mamba run -n gatk4 gatk MergeMutectStats \
        "${vcf_stats_array[@]}" \
        --output "${output_directory}/${sample_name}_mutect2.vcf.stats"
fi

# Learn read orientatation bias model
if [[ ! -f "${output_directory}/${sample_name}_mutect2_artifact_prior.tar.gz" ]]; then
    echo "Learning read orientation bias model..."
    f1r2_files=$(find "${input_directory}/f1r2" -type f | sed -e 's/^/-I /g' | tr '\n' ' ')
    read -r -a f1r2_files_array <<< "${f1r2_files}"
    mamba run -n gatk4 gatk LearnReadOrientationModel \
        "${f1r2_files_array[@]}" \
        --output "${output_directory}/${sample_name}_mutect2_artifact_prior.tar.gz"
    echo "...read orientation bias model trained."
fi
    
# Calculate cross sample contamination using pileups at known common germline variants.
if [[ ! -f "${output_directory}/${sample_name}_contamination.table" ]]; then
    echo "Getting contamination rate with CalculateContamination..."
    mamba run -n gatk4 gatk CalculateContamination \
        --input "${output_directory}/${sample_name}_pileups.table" \
        --output "${output_directory}/${sample_name}_contamination.table"
    echo "...contamination rate calculated."
fi

# Add FILTER column to Mutect2 VCF to identify variants which pass or fail filters.
filtered_vcf="${output_directory}/${sample_name}_mutect2_filtered.vcf"
if [[ ! -f "${filtered_vcf}" ]]; then
    echo "Filtering somatic variants with FilterMutectCalls..."
    mamba run -n gatk4 gatk FilterMutectCalls \
        --variant "${output_directory}/${sample_name}_mutect2.vcf" \
        --output "${filtered_vcf}" \
        --contamination-table "${output_directory}/${sample_name}_contamination.table" \
        --ob-priors  "${output_directory}/${sample_name}_mutect2_artifact_prior.tar.gz" \
        --reference "${reference_genome}"
    echo "...somatic variants filtered."
fi

# Filter variants that do not pass filter
filtered_vcf_filter_pass="${filtered_vcf//.vcf/_filter_pass.vcf.gz}"
if [[ ! -f "${filtered_vcf_filter_pass}" ]]; then
    filter_array=(
            'FILTER="PASS"'
            "INFO/DP>=${min_sequencing_depth}"
            "INFO/DP<=${max_sequencing_depth}"
            'TYPE="snp"'
            '((REF="C" & ALT="T") || (REF="T" & ALT="C"))'
        )
    filter_string=$(printf '%s & ' "${filter_array[@]}" | sed 's/ & $//g')
    bcftools view "${filtered_vcf}" --include "${filter_string}" --min-alleles 2 --max-alleles 2 --output "${filtered_vcf_filter_pass}"
    tabix --force --preset vcf "${filtered_vcf_filter_pass}"
fi

# Remove variants which are BRAVO TOPMed germline variants
filtered_vcf_bravo_pass="${filtered_vcf_filter_pass//.vcf.gz/_bravo_pass.vcf.gz}"
if [[ ! -f "${filtered_vcf_bravo_pass}" ]]; then
    bcftools view "${filtered_vcf_filter_pass}" --targets-file "^${bravo_variants}" --output "${filtered_vcf_bravo_pass}"
    tabix --force --preset vcf "${filtered_vcf_bravo_pass}"
fi

# Create BED file from VCF with window around reference
filtered_bed="${filtered_vcf//.vcf/.bed}"
if [[ ! -f "${filtered_bed}" ]]; then
    sequence_context_window_padding=$(( sequence_context_window_size / 2 ))
    bcftools view -H "${filtered_vcf_bravo_pass}" | \
        cut --fields 1-2 | \
        awk "{print \$1,\$2-${sequence_context_window_size}-1,\$2+${sequence_context_window_padding}}" | \
        tr ' ' '\t' > "${filtered_bed}"
fi

# Retrieve sequence context from reference FASTA using BED file
sequence_context="${filtered_bed//.bed/.txt.gz}"
if [[ ! -f "${sequence_context}" ]]; then
    paste <(bcftools view --no-header "${filtered_vcf_bravo_pass}" | cut --fields 1-2) \
        <(bedtools getfasta -fi "${reference_genome}" -bed "${filtered_bed}" | grep -v ">") | \
        bgzip > "${sequence_context}"
    tabix --sequence 1 --begin 2 --end 2 --force "${sequence_context}"
fi

# Create header line for sequence_context annotation
sequence_context_header="${filtered_bed//.bed/_sequence_context_header}"
echo "##INFO=<ID=sequence_context,Number=1,Type=String,Description=\"Sequence context\">" > "${sequence_context_header}"

# Annotate VCF with sequence context information
bcftools annotate "${filtered_vcf_bravo_pass}" -a "${sequence_context}" -c "CHROM,POS,sequence_context" -h "${sequence_context_header}" | bgzip | sponge "${filtered_vcf_bravo_pass}"
tabix --force --preset vcf "${filtered_vcf_bravo_pass}"
rm "${filtered_bed}" "${sequence_context_header}" "${sequence_context}" "${sequence_context}".tbi

# Remove variants which contain N in the sequence context
filtered_vcf_sequence_context_pass="${filtered_vcf_bravo_pass//.vcf.gz/_sequence_context_pass.vcf.gz}"
bcftools view "${filtered_vcf_bravo_pass}" --include 'INFO/sequence_context!~"N"' --output ${filtered_vcf_sequence_context_pass}
tabix --force --preset vcf "${filtered_vcf_sequence_context_pass}"
