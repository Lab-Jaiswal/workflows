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

##################################################################################################################################
#############################################---STEP 1: SET UP PARAMETERS---###################################################### 
##################################################################################################################################
options_array=(
    array_file
    output_directory
    annotated_output_directory
    varscan_min_coverage
    varscan_min_var_freq
    varscan_max_pvalue
    annovarroot
    mpileup_interval_bed
    run_varscan
    run_mutect
    run_funcotator
    run_annovar
    run_haplotypecaller
    bam_extension
    fastq_extension
    interval_list
    split_intervals
    normal_bam
    normal_pileups_table
    reference_genome
    gnomad_reference
    germline_snps
    assembly
    funcotator_sources
    transcript_list
    sequence_dictionary
    mutect_bam_output
    slurm_mode
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
            final_output_directory="${2}"; check_for_directory "${1}" "${2}"; shift 2 ;;
        --annotated_output_directory )
            annotated_output_directory="${2}"; check_for_directory "${1}" "${2}"; shift 2 ;;
        --bam_extension )
            bam_extension="${2}"; shift 2 ;;
        --fastq_extension )
            fastq_extension="${2}"; shift 2 ;;
        --assembly )
            assembly="${2}"; shift 2 ;;
        --reference_genome )
            reference_genome="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --slurm_mode )
            slurm_mode="${2}"; shift 2 ;;
        --run_mutect )
            run_mutect="${2}"; shift 2 ;;
        --interval_list )
            interval_list="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --split_intervals )
            split_intervals="${2}"; shift 2 ;;
        --sequence_dictionary )
            sequence_dictionary="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --mutect_bam_output )
            mutect_bam_output="${2}"; shift 2 ;;
        --gnomad_reference )
            gnomad_reference="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --normal_bam )
            normal_bam_file="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --normal_pileups_table )
            normal_pileups_table="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --run_funcotator )
            run_funcotator="${2}"; shift 2 ;;
        --funcotator_sources )
            funcotator_sources="${2}"; check_for_directory "${1}" "${2}"; shift 2 ;;
        --transcript_list )
            transcript_list="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --run_varscan )
            run_varscan="${2}"; shift 2 ;;
        --mpileup_interval_bed )
            mpileup_interval_bed="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
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
        --run_haplotypecaller )
            run_haplotypecaller="${2}"; shift 2 ;;
        --germline_snps )
            germline_snps="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        -- )
            shift; break ;;
        * )
            echo "Invalid argument ${1} ${2}" >&2
            exit 1
    esac
done

if [[ ${slurm_mode} == true ]]; then
    line_number=${SLURM_ARRAY_TASK_ID} #get index of which file to process from $SLURM_ARRAY_TASK_ID provided by SLURM
else
    line_number=${TASK_ID}
fi

array_prefix="$(sed "${line_number}q; d" "${array_file}")" #extract only the line number corresponding to $SLURM_ARRAY_TASK_ID
sample_name=$(basename "${array_prefix}")
bam_file="${array_prefix}.${bam_extension}"
gatk_command="mamba run -n gatk4 gatk"
code_directory=$(realpath $(dirname ${BASH_SOURCE[0]}))

if [[ $final_output_directory != "none" ]]; then
    final_output_directory=${final_output_directory}/${sample_name}
fi

if [[ ${slurm_mode} == true ]]; then
    input_directory="${TMPDIR}/inputs"
    output_directory="${TMPDIR}/outputs/${sample_name}"
    references_directory="${TMPDIR}/references"
    if [[ $normal_bam_file != "none" ]]; then
        normal_input_directory="${TMPDIR}/inputs/normal"
    fi

    rsync_command="rsync --human-readable --links --mkpath --progress --recursive --times --update --verbose"
    mkdir -p "${final_output_directory}"
    ${rsync_command} --exclude "logs" "${final_output_directory}/" "${output_directory}"

    #if [[ ${fastq_extension} != "none" ]]; then
        #original_bam_file="${bam_file}"
        #bam_file="${input_directory}/$(basename ${original_bam_file})"
        #${rsync_command} "${original_bam_file}" "${bam_file}"
        #check_for_file "${rsync_command}" "${bam_file}"

        #if [[ ${bam_extension} == "bam" ]]; then
            #original_index_file="${original_bam_file}.bai"
        #else
            #original_index_file="${original_bam_file%.*}.crai"
        #fi
        #index_file="${input_directory}/$(basename ${original_index_file})"
        #${rsync_command} "${original_index_file}" "${index_file}"
        #check_for_file "${rsync_command}" "${index_file}"
    #fi

    #original_reference_genome="${reference_genome}"
    #reference_genome="${references_directory}/$(basename ${original_reference_genome})"
    #${rsync_command} ${original_reference_genome%.*}* "${references_directory}"
    #check_for_file "${rsync_command}" "${reference_genome}"

    #if [[ ${interval_list} != "none" ]]; then
        #original_interval_list="${interval_list}"
        #interval_list="${references_directory}/$(basename ${original_interval_list})"
        #${rsync_command} "${original_interval_list}" "${interval_list}"
        #check_for_file "${rsync_command}" "${interval_list}"
    #fi

    #if [[ ${gnomad_reference} != "none" ]]; then
        #original_gnomad_reference="${gnomad_reference}"
        #gnomad_reference="${references_directory}/$(basename ${original_gnomad_reference})"
        #${rsync_command} "${original_gnomad_reference}" "${gnomad_reference}"
        #${rsync_command} "${original_gnomad_reference}.tbi" "${gnomad_reference}.tbi"
        #check_for_file "${rsync_command}" "${gnomad_reference}"
        #check_for_file "${rsync_command}" "${gnomad_reference}.tbi"
    #fi

    #if [[ ${sequence_dictionary} != "none" ]]; then
        #original_sequence_dictionary="${sequence_dictionary}"
        #sequence_dictionary="${references_directory}/$(basename ${original_sequence_dictionary})"
        #${rsync_command} "${original_sequence_dictionary}" "${sequence_dictionary}"
        #check_for_file "${rsync_command}" "${sequence_dictionary}"
    #fi

    #if [[ ${transcript_list} != "none" ]]; then
        #original_transcript_list="${transcript_list}"
        #transcript_list="${references_directory}/$(basename ${original_transcript_list})"
        #${rsync_command} "${original_transcript_list}" "${transcript_list}"
        #check_for_file "${rsync_command}" "${transcript_list}"
    #fi
    
    #if [[ ${funcotator_sources} != "none" ]]; then
        #original_funcotator_sources="${funcotator_sources}"
        #funcotator_sources="${references_directory}/funcotator_sources/$(basename ${original_funcotator_sources})"
        #${rsync_command} "${original_funcotator_sources}/" "${funcotator_sources}"
        #check_for_directory "${rsync_command}" "${funcotator_sources}"
    #fi

    #if [[ ${normal_bam_file} != "none" ]]; then
        #original_normal_bam_file="${normal_bam_file}"
        #normal_bam_file="${normal_input_directory}/$(basename ${original_normal_bam_file})"
        #${rsync_command} "${original_normal_bam_file}" "${normal_bam_file}"
        #check_for_file "${rsync_command}" "${normal_bam_file}"

        #if [[ ${bam_extension} == "bam" ]]; then
            #original_normal_index_file="${original_normal_bam_file}.bai"
        #else
            #original_normal_index_file="${original_normal_bam_file%.*}.crai"
        #fi
        #normal_index_file="${normal_input_directory}/$(basename ${original_normal_index_file})"
        #${rsync_command} "${original_normal_index_file}" "${normal_index_file}"
        #check_for_file "${rsync_command}" "${normal_index_file}"

        #original_normal_pileups_table="${normal_pileups_table}"
        #normal_pileups_table="${normal_input_directory}/$(basename ${original_normal_pileups_table})"
        #${rsync_command} "${original_normal_pileups_table}" "${normal_pileups_table}"
        #check_for_file "${rsync_command}" "${normal_pileups_table}"
    #fi
else
    output_directory=${final_output_directory}
fi

##################################################################################################################################
#################################################---STEP 2: MUTECT.sh---########################################################## 
##################################################################################################################################       

if [[ ${fastq_extension} != "none" ]]; then
    fastq_read1="${array_prefix}_R1_001.fastq.gz"
    fastq_read2="${array_prefix}_R2_001.fastq.gz"
    read_group="@RG\tID:${sample_name}\tLB:${sample_name}\tPL:illumina\tSM:${sample_name}"
    sample_name="${sample_name}_${assembly}"

    "${code_directory}/fastq_to_bam.sh" \
        --sample_name "${sample_name}" \
        --read_group "${read_group}" \
        --reference_genome "${reference_genome}" \
        --fastq_read1 "${fastq_read1}" \
        --fastq_read2 "${fastq_read2}" \
        --output_directory "${output_directory}"

    bam_file="${output_directory}/${sample_name}.bam"

    if [[ ${slurm_mode} == true ]]; then
        ${rsync_command} "${output_directory}/" "${final_output_directory}"
        check_for_file "${rsync_command}" "${final_output_directory}/${sample_name}.bam"
    fi
fi

if [[ ${run_mutect} == true ]]; then
    if [[ ${split_intervals} == true ]]; then
        num_intervals=$(grep -c -v "@" < "${interval_list}")
        echo "Number of intervals: $num_intervals"
        seq 1 "${num_intervals}" | parallel -j8 --progress --ungroup \
            "${code_directory}/mutect_and_pileups.sh" \
                --bam_file "${bam_file}" \
                --normal_bam_file "${normal_bam_file}" \
                --interval_list "${interval_list}" \
                --interval_number {} \
                --reference_genome "${reference_genome}" \
                --gnomad_reference "${gnomad_reference}" \
                --output_directory "${output_directory}" \
                --mutect_bam_output "${mutect_bam_output}" \
                --gatk_command "${gatk_command}"
           
        if [[ ${slurm_mode} == true ]]; then
            ${rsync_command} "${output_directory}/" "${final_output_directory}"
            check_for_file "${rsync_command}" "${final_output_directory}/${sample_name}_mutect2.vcf"
        fi
        
        if [[ ! -f "${output_directory}/${sample_name}_pileups.table" ]]; then
            pileup_tables=$(find "${output_directory}/pileups" -type f | grep -E ".*_pileups.table$" | sort -V | sed -e 's/^/--input /g' | tr '\n' ' ')
            read -r -a pileup_tables_array <<< "${pileup_tables}"
            ${gatk_command} GatherPileupSummaries \
                --sequence-dictionary "${sequence_dictionary}" \
                "${pileup_tables_array[@]}" \
                --output "${output_directory}/${sample_name}_pileups.table"

            if [[ ${slurm_mode} == true ]]; then
                ${rsync_command} "${output_directory}/" "${final_output_directory}"
                check_for_file "${rsync_command}" "${final_output_directory}/${sample_name}_pileups.table"
            fi
        fi

        if [[ ! -f "${output_directory}/${sample_name}_mutect2.vcf" ]]; then
            vcf_files=$(find "${output_directory}/vcfs" -type f | grep -E ".*.vcf$" | sort -V | sed -e 's/^/--INPUT /g' | tr '\n' ' ')
            read -r -a vcf_files_array <<< "${vcf_files}"
            ${gatk_command} MergeVcfs \
                "${vcf_files_array[@]}" \
                --OUTPUT "${output_directory}/${sample_name}_mutect2.vcf"

            if [[ ${slurm_mode} == true ]]; then
                ${rsync_command} "${output_directory}/" "${final_output_directory}"
                check_for_file "${rsync_command}" "${final_output_directory}/${sample_name}_mutect2.vcf"
            fi
        fi

        if [[ ! -f "${output_directory}/${sample_name}_mutect2.vcf.stats" ]]; then
            vcf_stats=$(find "${output_directory}/vcfs" -type f | grep -E ".*.vcf.stats$" | sort -V | sed -e 's/^/--stats /g' | tr '\n' ' ')
            read -r -a vcf_stats_array <<< "${vcf_stats}"
            ${gatk_command} MergeMutectStats \
                "${vcf_stats_array[@]}" \
                --output "${output_directory}/${sample_name}_mutect2.vcf.stats"

            if [[ ${slurm_mode} == true ]]; then
                ${rsync_command} "${output_directory}/" "${final_output_directory}"
                check_for_file "${rsync_command}" "${final_output_directory}/${sample_name}_mutect2.vcf.stats"
            fi
        fi

        if [[ ! -f "${output_directory}/${sample_name}_mutect2_artifact_prior.tar.gz" ]]; then
            f1r2_files=$(find "${output_directory}/f1r2" -type f | sed -e 's/^/-I /g' | tr '\n' ' ')
            read -r -a f1r2_files_array <<< "${f1r2_files}"
            ${gatk_command} LearnReadOrientationModel \
                "${f1r2_files_array[@]}" \
                --output "${output_directory}/${sample_name}_mutect2_artifact_prior.tar.gz"

            if [[ ${slurm_mode} == true ]]; then
                ${rsync_command} "${output_directory}/" "${final_output_directory}"
                check_for_file "${rsync_command}" "${final_output_directory}/${sample_name}_mutect2_artifact_prior.tar.gz"
            fi
        fi
    else
        "${code_directory}/mutect_and_pileups.sh" \
            --bam_file "${bam_file}" \
            --normal_bam_file "${normal_bam_file}" \
            --interval_list "${interval_list}" \
            --reference_genome "${reference_genome}" \
            --gnomad_reference "${gnomad_reference}" \
            --output_directory "${output_directory}" \
            --mutect_bam_output "${mutect_bam_output}" \
            --gatk_command "${gatk_command}" \
            --run_mutect "${run_mutect}"
        
        if [[ ! -f "${output_directory}/${sample_name}_mutect2_artifact_prior.tar.gz" ]]; then
            ${gatk_command} LearnReadOrientationModel \
                --input "${output_directory}/${sample_name}_f1r2.tar.gz" \
                --output "${output_directory}/${sample_name}_mutect2_artifact_prior.tar.gz"
        fi

        if [[ $slurm_mode == true ]]; then
            ${rsync_command} "${output_directory}/" "${final_output_directory}"
            check_for_file "${rsync_command}" "${final_output_directory}/${sample_name}_pileups.table"
            check_for_file "${rsync_command}" "${final_output_directory}/${sample_name}_mutect2.vcf"
            check_for_file "${rsync_command}" "${final_output_directory}/${sample_name}_mutect2.vcf.stats"
            check_for_file "${rsync_command}" "${final_output_directory}/${sample_name}_mutect2_artifact_prior.tar.gz"
        fi
    fi
    
    "${code_directory}/mutect_filter.sh" \
        --vcf_file "${output_directory}/${sample_name}_mutect2.vcf" \
        --pileups_table "${output_directory}/${sample_name}_pileups.table" \
        --normal_pileups_table "${normal_pileups_table}" \
        --orientation_bias_priors "${output_directory}/${sample_name}_mutect2_artifact_prior.tar.gz" \
        --reference_genome "${reference_genome}" \
        --gatk_command "${gatk_command}" \

    if [[ ${slurm_mode} == true ]]; then
        ${rsync_command} "${output_directory}/" "${final_output_directory}"
        check_for_file "${rsync_command}" "${final_output_directory}/${sample_name}_mutect2_filtered.vcf"
    fi

    if [[ ${run_funcotator} == true ]]; then
        "${code_directory}/funcotator.sh" \
            --filtered_vcf "${output_directory}/${sample_name}_mutect2_filtered.vcf" \
            --annotated_output_directory= "${annotated_output_directory}" \
            --reference_genome "${reference_genome}" \
            --funcotator_sources "${funcotator_sources}" \
            --transcript_list "${transcript_list}" \
            --gatk_command "${gatk_command}"
        if [[ ${slurm_mode} == true ]]; then
            ${rsync_command} "${output_directory}/" "${final_output_directory}"
            check_for_file "${rsync_command}" "${final_output_directory}/${sample_name}_mutect2_filtered_funcotator.vcf"
        fi
    fi
        
else
    echo "Mutect2 calling of somatic variants not requested"
fi
            
##################################################################################################################################
############################################---STEP 3: HAPLOTYPECALLER.sh---######################################################
##################################################################################################################################    
if [[ $run_haplotypecaller = true ]] ; then
    "${code_directory}/haplotypecaller.sh" \
        --bam_file "${bam_file}" \
        --reference_genome "${reference_genome}" \
        --germline_snps "${germline_snps}" \
        --gatk_command "${gatk_command}"

    if [[ ${slurm_mode} == true ]]; then
        ${rsync_command} "${output_directory}/" "${final_output_directory}"
        check_for_file "${rsync_command}" "${final_output_directory}/${sample_name}_haplotypecaller.bam"
        check_for_file "${rsync_command}" "${final_output_directory}/${sample_name}_haplotypecaller.vcf"
    fi
else
    echo "HaplotypeCaller analysis not requested"
fi

##################################################################################################################################
################################################---STEP 4: VARSCAN.sh---########################################################## 
##################################################################################################################################    
if [[ ${run_varscan} = true ]]; then
    echo "Varscan analysis requested..."
    "${code_directory}/varscan.sh" \
        --bam_file "${bam_file}" \
        --reference_genome "${reference_genome}" \
        --varscan_min_coverage "${varscan_min_coverage}" \
        --varscan_min_var_freq "${varscan_min_var_freq}" \
        --varscan_max_pvalue "${varscan_max_pvalue}" \
        --run_annovar "${run_annovar}" \
        --annovarroot "${annovarroot}" \
        --mpileup_interval_bed "${mpileup_interval_bed}"
    echo "...Varscan analysis complete"

    if [[ ${slurm_mode} == true ]]; then
        ${rsync_command} "${output_directory}/" "${final_output_directory}"
        check_for_file "${rsync_command}" "${final_output_directory}/${sample_name}.pileup"
        check_for_file "${rsync_command}" "${final_output_directory}/${sample_name}_varscan2.vcf"
        check_for_file "${rsync_command}" "${final_output_directory}/${sample_name}_varscan2_filter.vcf"
        check_for_file "${rsync_command}" "${final_output_directory}/${sample_name}_varscan2_filter_annovar.hg38_multianno.vcf"
    fi
else 
    echo "Varscan analysis not requested"
fi
