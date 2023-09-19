#!/usr/bin/env bash

# Run scripts to enable activating conda environments
. "${HOME}/micromamba/etc/profile.d/conda.sh"
. "${HOME}/micromamba/etc/profile.d/mamba.sh"
mamba activate base
mamba activate bcftools

# Set bash options for verbose output and to fail immediately on errors or if variables are undefined.
set -o xtrace -o nounset -o errexit

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
    vcf_directory
    vcf_directory_list
    minimum_dp
    minimum_ad
    minimum_af
    filter_silent
    filter_pass
    filter_exclusions
    filter_sb
    reannotate
    funcotator_sources
    transcript_list
    reference_genome
    output_directory
)

eval "$(printf "%s\n" "${options_array[@]}" | xargs --replace=% echo "declare %=none;")"
longoptions=$(echo "${options_array[@]}" | sed -e 's/ /:,/g' | sed -e 's/$/:/')

# Parse command line arguments with getopt
arguments=$(getopt --options a --longoptions "${longoptions}" --name 'filter_variants_targeted' -- "$@")
eval set -- "${arguments}"

while true; do
    case "${1}" in
        --vcf_directory )
            vcf_directory="${2}"; check_for_directory "${1}" "${2}"; shift 2 ;;
        --vcf_directory_list )
            vcf_directory_list="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --minimum_dp )
            minimum_dp="${2}"; shift 2 ;;
        --minimum_ad )
            minimum_ad="${2}"; shift 2 ;;
        --minimum_af )
            minimum_af="${2}"; shift 2 ;;
        --filter_silent )
            filter_silent="${2}"; shift 2 ;;
        --filter_pass )
            filter_pass="${2}"; shift 2 ;;
        --filter_exclusions )
            filter_exclusions="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --filter_sb )
            filter_sb="${2}"; shift 2 ;;
        --reannotate )
            reannotate="${2}"; shift 2 ;;
        --funcotator_sources )
            funcotator_sources="${2}"; shift 2 ;;
        --transcript_list )
            transcript_list="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --reference_genome )
            reference_genome="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --output_directory )
            output_directory="${2}"; shift 2 ;;
        -- )
            shift; break;;
        * )
            echo "Invalid argument ${1} ${2}" >&2
            exit 1
    esac
done

reannotate() {
    vcf_file="${1}"
    reference_genome="${2}"
    funcotator_sources="${3}"
    transcript_list="${4}"
    output_directory="${5}"

    vcf_basename="$(basename ${vcf_file})"
    filter_vcf="${output_directory}/${vcf_basename//_funcotator/}"
    reheader_vcf="${output_directory}/${vcf_basename//_funcotator/_reheader}"
    new_header="${output_directory}/${vcf_basename//_funcotator.vcf.gz/_new_header}"
    funcotator_vcf="${output_directory}/${vcf_basename}"
    echo ${filter_vcf}

    funcotator_cols=$(bcftools view --header-only ${vcf_file} | grep "INFO=" | grep funcotator | sed 's/.*<ID=/INFO\//g' | sed 's/,.*$//g' | tr '\n' ',' | sed 's/,$//g')
    bcftools annotate --remove "INFO/FUNCOTATION,${funcotator_cols}" ${vcf_file} --output ${filter_vcf}
    tabix --preset vcf --force ${filter_vcf}

    bcftools view --header-only ${filter_vcf} | grep -v Funcotator > ${new_header}
    bcftools reheader --header ${new_header} ${filter_vcf} --output ${reheader_vcf}
    tabix --preset vcf --force ${reheader_vcf}

    mamba run -n gatk4 gatk Funcotator \
        --variant "${reheader_vcf}" \
        --reference "${reference_genome}" \
        --ref-version hg38 \
        --data-sources-path "${funcotator_sources}" \
        --transcript-list "${transcript_list}" \
        --output "${funcotator_vcf}" \
        --output-file-format VCF 
    tabix --preset vcf --force ${funcotator_vcf}

    funcotator_header=$(bcftools view --header-only "${funcotator_vcf}" | grep "ID=FUNCOTATION" | sed 's/.*Description="//g' | sed 's/">//g')
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
        bgzip --stdout | sponge "${funcotator_vcf}"
    tabix --preset vcf --force "${funcotator_vcf}"

    bcftools +split-vep --annotation "FUNCOTATION" --columns "${funcotator_columns}" --annot-prefix "funcotator_" "${funcotator_vcf}" | \
        sed "s/Number=\./Number=A/g" | \
        bgzip --stdout | sponge "${funcotator_vcf}"
    tabix --preset vcf --force "${funcotator_vcf}"
}

export -f reannotate

query_vcf() {
    vcf_file="${1}"
    minimum_dp="${2}"
    minimum_ad="${3}"
    minimum_af="${4}"
    filter_silent="${5}"
    filter_pass=${6}
    filter_exclusions="${7}"
    filter_sb="${8}"
    sample_tsv="${vcf_file//.vcf.gz/.tsv}"

    readarray -t info_columns < <(bcftools view --header "${vcf_file}" | \
        grep "INFO=" | \
        grep --invert-match --extended-regexp "funcotator|FUNCOTATION" | \
        sed 's/##INFO=<ID=//g' | \
        sed 's/,.*$//g' | \
        sort --unique)
    readarray -t funcotator_columns < <(bcftools view --header "${vcf_file}" | \
        grep "INFO=" | \
        grep "funcotator" | \
        grep --extended-regexp "Gencode|Consequence" | \
        grep --invert-match "otherTranscripts" | \
        sed 's/##INFO=<ID=//g' | \
        sed 's/,.*$//g')
    readarray -t format_columns < <(bcftools view --header "${vcf_file}" | \
        grep "FORMAT=" | \
        sed 's/##FORMAT=<ID=//g' | \
        sed 's/,.*$//g')

    info_string=$(printf '%%INFO/%s\\t' "${info_columns[@]}" | sed 's/\\t$//g')
    funcotator_string=$(printf '%%INFO/%s\\t' "${funcotator_columns[@]}" | sed 's/\\t$//g')
    format_string=$(printf '%%%s\\t' "${format_columns[@]}" | sed 's/\\t$//g')
    query_string=$(printf '[%s]\\t%s\\t%s\\t[%s]\\t%s' "%SAMPLE" "%FILTER" "${funcotator_string}" "${format_string}" "${info_string}")

    echo "${query_string}" | sed 's/\\t/\t/g' | \
        tr --delete '%[]' | \
        sed 's?INFO/funcotator_??g' | \
        sed 's?INFO/DP?info_dp?g' | \
        sed 's?INFO/??g' | \
        tr '[:upper:]' '[:lower:]' > "${sample_tsv}"

    declare filter_condition=()

    if [[ ${minimum_dp} != "none" ]]; then
        filter_condition+=("INFO/DP < ${minimum_dp} || FORMAT/DP < ${minimum_dp}")
    fi

    if [[ ${minimum_ad} != "none" ]]; then
        filter_condition+=("FORMAT/AD[0:*] < ${minimum_ad}")
    fi

    if [[ ${minimum_af} != "none" ]]; then
        filter_condition+=("FORMAT/AF < ${minimum_af}")
    fi

    if [[ ${filter_silent} == true ]]; then
        filter_condition+=('(INFO/funcotator_Consequence !~ "NONSENSE" && INFO/funcotator_Consequence !~ "MISSENSE" && INFO/funcotator_Consequence !~ "SPLICE_SITE" && INFO/funcotator_Consequence !~ "FRAME_SHIFT_DEL" && INFO/funcotator_Consequence !~ "FRAME_SHIFT_INS" && INFO/funcotator_Consequence !~ "IN_FRAME_DEL" && INFO/funcotator_Consequence !~ "IN_FRAME_INS")')
    fi

    if [[ ${filter_sb} == true ]]; then
        filter_condition+=("FORMAT/SB[0:*] == 0 || FORMAT/F1R2[0:*] == 0 || FORMAT/F2R1[0:*] == 0")
    fi

    if [[ ${filter_pass} == true ]]; then
        filter_condition+=('FILTER != "PASS"')
    fi

    if [[ ${filter_exclusions} != "none" ]]; then
        filter_condition+=("$(printf 'FILTER ~ "%s" || ' $(cat "${filter_exclusions}") | sed 's/ || $//g')")
    fi

    declare filter_string=""
    if [[ ((${#filter_condition[@]} -gt 0 )) ]]; then
        filter_string_concat="$(printf '%s || ' "${filter_condition[@]}" | sed 's/ || $//g')"
        filter_string=("--exclude" "${filter_string_concat}")
    fi

    bcftools query "${filter_string[@]}" "${vcf_file}" --format "${query_string}\\n" >> "${sample_tsv}"
    sample_tsv_nrows=$(tail -n +2 ${sample_tsv} | wc --lines)
    sample_file_name="$(basename ${vcf_file//.vcf/})"
    paste <(echo -e "file_name\n$(yes "${sample_file_name}" | head -n ${sample_tsv_nrows})") "${sample_tsv}" | sponge "${sample_tsv}"
}

export -f query_vcf

query_vcf_varscan() {
    vcf_file="${1}"
    minimum_dp="${2}"
    minimum_ad="${3}"
    minimum_af="${4}"
    filter_silent="${5}"
    filter_pass=${6}
    filter_exclusions="${7}"
    filter_sb="${8}"
    sample_tsv="${vcf_file//.vcf.gz/.tsv}"

    readarray -t info_columns < <(bcftools view --header "${vcf_file}" | \
        grep "INFO=" | \
        grep --invert-match --extended-regexp "AAChange.ensGene|ALLELE_END" | \
        sed 's/##INFO=<ID=//g' | \
        sed 's/,.*$//g' | \
        sort --unique)
    readarray -t format_columns < <(bcftools view --header "${vcf_file}" | \
        grep "FORMAT=" | \
        sed 's/##FORMAT=<ID=//g' | \
        sed 's/,.*$//g')

    info_string=$(printf '%%INFO/%s\\t' "${info_columns[@]}" | sed 's/\\t$//g')
    format_string=$(printf '%%%s\\t' "${format_columns[@]}" | sed 's/\\t$//g')
    query_string=$(printf '%s\\t%s\\t%s\\t%s\\t[%s]\\t%s' "%FILTER" "%CHROM" "%REF" "%ALT" "${format_string}" "${info_string}")

    echo "${query_string}" | sed 's/\\t/\t/g' | \
        tr --delete '%[]' | \
        sed 's?INFO/??g' | \
        tr '[:upper:]' '[:lower:]' > "${sample_tsv}"

    declare filter_condition=()

    if [[ ${minimum_dp} != "none" ]]; then
        filter_condition+=("INFO/ADP < ${minimum_dp} || FORMAT/DP < ${minimum_dp}")
    fi

    if [[ ${minimum_ad} != "none" ]]; then
        filter_condition+=("FORMAT/AD[0:*] < ${minimum_ad}")
    fi

    if [[ ${minimum_af} != "none" ]]; then
        filter_condition+=("FORMAT/AF < ${minimum_af}")
    fi

    if [[ ${filter_silent} == true ]]; then
        filter_condition+=('(INFO/Func.ensGene !~ "exonic" && INFO/Func.ensGene !~ "splice")')
    fi

    if [[ ${filter_sb} == true ]]; then
        filter_condition+=("FORMAT/RDF[0:*] == 0 || FORMAT/RDR[0:*] == 0 || FORMAT/ADF[0:*] == 0 || FORMAT/ADR[0:*] == 0")
    fi

    if [[ ${filter_pass} == true ]]; then
        filter_condition+=('FILTER != "PASS"')
    fi

    if [[ ${filter_exclusions} != "none" ]]; then
        filter_condition+=("$(printf 'FILTER ~ "%s" || ' $(cat "${filter_exclusions}") | sed 's/ || $//g')")
    fi

    declare filter_string=""
    if [[ ((${#filter_condition[@]} -gt 0 )) ]]; then
        filter_string_concat="$(printf '%s || ' "${filter_condition[@]}" | sed 's/ || $//g')"
        filter_string=("--exclude" "${filter_string_concat}")
    fi

    #bcftools query "${filter_string[@]}" "${vcf_file}" --format "${query_string}\\n" >> "${sample_tsv}"
    bcftools query "${vcf_file}" --format "${query_string}\\n" >> "${sample_tsv}"
    sample_tsv_nrows=$(tail -n +2 ${sample_tsv} | wc --lines)
    sample_file_name="$(basename ${vcf_file//.vcf.gz/})"
    paste <(echo -e "file_name\n$(yes "${sample_file_name}" | head -n ${sample_tsv_nrows})") "${sample_tsv}" | sponge "${sample_tsv}"
}

export -f query_vcf_varscan

add_file_name() {
    vcf_file="${1}"
    pileup_region_file="${1}"
    fixed_pileup_region="${pileup_region_file//.pileup_region/_file_name.pileup_region}"
    pileup_region_file_name="$(basename ${pileup_region_file//.vcf/})"
    pileup_region_nrows="$(wc -l "${pileup_region_file}" | cut -d ' ' -f 1)"
    
    paste <(yes "${pileup_region_file_name}" | head -n ${pileup_region_nrows}) "${pileup_region_file}" > ${fixed_pileup_region}
}

export -f add_file_name

if [[ ${vcf_directory} == "none" ]]; then
    line_number=${SLURM_ARRAY_TASK_ID}
    vcf_directory="$(sed "${line_number}q; d" "${vcf_directory_list}")"
    aggregated_file="${output_directory}/aggregated_${line_number}.tsv"
    aggregated_file_varscan="${output_directory}/aggregated_varscan_${line_number}.tsv"
    aggregated_pileup_file="${output_directory}/aggregated_${line_number}.pileup_region" 
else
    aggregated_file="${output_directory}/aggregated.tsv"
    aggregated_file_varscan="${output_directory}/aggregated_varscan.tsv"
    aggregated_pileup_file="${output_directory}/aggregated.pileup_region" 
fi

if [[ ${reannotate} == true ]]; then
    vcf_directory_dirname=$(dirname ${vcf_directory})
    vcf_directory_basename=$(basename ${vcf_directory})
    reannotated_directory="${vcf_directory_dirname}_reannotated/${vcf_directory_basename}_reannotated"
    mkdir -p ${reannotated_directory}
    find "${vcf_directory}" -name "*.vcf.gz" | sort --unique | \
        xargs --replace=% bash -c "reannotate % ${reference_genome} ${funcotator_sources} ${transcript_list} ${reannotated_directory}"
    vcf_directory=${reannotated_directory}
fi

find "${vcf_directory}" -name "*_funcotator.vcf.gz" | sort --unique | \
    xargs --replace=% bash -c "query_vcf % ${minimum_dp} ${minimum_ad} ${minimum_af} ${filter_silent} ${filter_pass} ${filter_exclusions} ${filter_sb}"

find "${vcf_directory}" -name "*_multianno.vcf.gz" | sort --unique | \
    xargs --replace=% bash -c "query_vcf_varscan % ${minimum_dp} ${minimum_ad} ${minimum_af} ${filter_silent} ${filter_pass} ${filter_exclusions} ${filter_sb}"

find "${vcf_directory}" -name "*.pileup_region" | grep --invert-match "file_name" | sort -u | \
    xargs --replace=% bash -c "add_file_name %"

mkdir -p ${output_directory}
head -n 1 "$(find "${vcf_directory}" -name "*filtered_funcotator.tsv" | sort --unique | head -n 1)" > "${aggregated_file}"
find "${vcf_directory}" -name "*filtered_funcotator.tsv" -print0 | xargs --null --replace=% bash -c "tail -n +2 % >> ${aggregated_file}"

head -n 1 "$(find "${vcf_directory}" -name "*_multianno.tsv" | sort --unique | head -n 1)" > "${aggregated_file_varscan}"
find "${vcf_directory}" -name "*_multianno.tsv" -print0 | xargs --null --replace=% bash -c "tail -n +2 % >> ${aggregated_file_varscan}"

echo -e "file_name\tsample\tchr\tpos\tref\talt\tprotein_change\tref_count\talt_count" > "${aggregated_pileup_file}"
find "${vcf_directory}" -name "*_file_name.pileup_region" -print0 | xargs --null --replace=% bash -c "cat % >> ${aggregated_pileup_file}"
