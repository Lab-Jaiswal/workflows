#!/usr/bin/env bash

function run_job() {
    set -o pipefail -o errexit
    export PATH="${HOME}/.local/bin:${PATH}"
    curl micro.mamba.pm/install.sh | bash
    micromamba shell init --shell=bash "${HOME}/micromamba"
    export MAMBA_EXE="${HOME}/.local/bin/micromamba"
    export MAMBA_ROOT_PREFIX="${HOME}/micromamba"
    . "${HOME}/micromamba/etc/profile.d/micromamba.sh"
    micromamba activate
    micromamba config prepend channels defaults
    micromamba config prepend channels bioconda
    micromamba config prepend channels conda-forge
    micromamba config prepend channels dnachun
    micromamba config set auto_activate_base true
    micromamba install --yes \
        bash \
        coreutils \
        findutils \
        git \
        grep \
        parallel \
        sed \
        tar \
        util-linux
    micromamba create -n gatk gatk4
    conda init bash
    . "${HOME}/micromamba/etc/profile.d/conda.sh"
    mamba init bash
    conda config --set auto_stack 1
    set -o nounset -o xtrace

    cd "${HOME}"
    git clone --branch "${git_branch}" https://github.com/Lab-Jaiswal/workflows
    batch_size=$(( 2*batch_size ))

    local_output_directory="${HOME}/${output_directory}"
    logs_directory=${local_output_directory}/logs
    mkdir -p "${logs_directory}"
    log_file="${logs_directory}/job_${array_number}.log"
    
    # Download references
    local_references_directory="${HOME}/${references_directory}"
    dx download --recursive --overwrite "${project_id}:/${references_directory}/" --output "${local_references_directory}" &> "${log_file}"
    
    # Download list of files if provided
    # If a meta list is provided instead, download the file list from the meta list by extracting the row corresponding the array number.
    local_sample_list_directory="${HOME}/${sample_list_directory}"
    mkdir -p "${local_sample_list_directory}"
    local_sample_list="${local_sample_list_directory}/initial_sample_list"

    if [[ -z ${sample_list} ]] && [[ -n ${meta_sample_list} ]]; then
        sample_list=$(sed "${array_number}"q; d "${meta_sample_list}")
    fi
    dx download --recursive --overwrite "${project_id}:/${sample_list}" --output "${local_sample_list}" &> "${log_file}"

    local_input_directory="${HOME}/${input_directory}"
    mkdir -p "${local_input_directory}"
    filtered_list_of_samples="${local_input_directory}/filtered_list_of_samples"

    if [[ ${batch_size} != 0 ]]; then
        head -n ${batch_size} "${local_sample_list}" > "${filtered_list_of_samples}"
    else
        cp "${local_sample_list}" "${filtered_list_of_samples}"
    fi

    xargs --replace=% bash -c "dx download --overwrite ${project_id}:/${sample_directory}/% --output ${local_input_directory}/%" < "${filtered_list_of_samples}" &> "${log_file}"
        
    "${HOME}/workflows/BWA_CHIP/submit_BWA_CHIP.sh" \
        --input_directory "${local_input_directory}" \
        --output_directory "${local_output_directory}" \
        --run_mutect "${run_mutect}" \
        --mutect_bam_out "${mutect_bam_out}" \
        --run_haplotypecaller "${run_haplotypecaller}" \
        --run_varscan "${run_varscan}" \
        --reference_genome "${local_references_directory}/${reference_genome}" \
        --gnomad_reference "${local_references_directory}/${gnomad_reference}" \
        --funcotator_sources "${local_references_directory}/${funcotator_sources}" \
        --interval_list "${interval_list}" \
        --n_jobs "${n_jobs}" \
        --bam_extension "cram" \
        --assembly "GRCh38" &> "${log_file}"
    output_tar="${HOME}/outputs_${array_number}.tar"
    tar --create --file "${output_tar}" "${local_output_directory}"
    outputs_folder_id=$(dx upload "${output_tar}" --brief)
    dx-jobutil-add-output outputs_tar "${outputs_folder_id}" --class=file
}

function submit_job() {
    array_number=${1}
    job_args=${2}
    instance_type=${3}

    job_id=$(dx-jobutil-new-job run_job -i"${array_number}" "${job_args}" --instance-type="${instance_type}")
    dx-jobutil-add-output outputs_folder "${job_id}:outputs_tar" --class=jobref --array
}

function main() {
    set -o pipefail -o errexit
    export PATH="${HOME}/.local/bin:${PATH}"
    curl micro.mamba.pm/install.sh | bash
    micromamba shell init --shell=bash "${HOME}/micromamba"
    export MAMBA_EXE="${HOME}/.local/bin/micromamba"
    export MAMBA_ROOT_PREFIX="${HOME}/micromamba"
    . "${HOME}/micromamba/etc/profile.d/micromamba.sh"
    micromamba activate
    micromamba config prepend channels defaults
    micromamba config prepend channels bioconda
    micromamba config prepend channels conda-forge
    micromamba config prepend channels dnachun
    micromamba config set auto_activate_base true
    micromamba install --yes \
        bash \
        coreutils \
        findutils
    conda init bash
    . "${HOME}/micromamba/etc/profile.d/conda.sh"
    mamba init bash
    conda config --set auto_stack 1
    set -o nounset -o xtrace

    dx-download-all-inputs 
    
    if [[ ${number_of_batches} == 0 && -n ${meta_sample_list} ]]; then
        local_sample_list_directory="${HOME}/${sample_list_directory}"
        mkdir -p "${local_sample_list_directory}"

        meta_sample_list="${sample_list_directory}/${meta_sample_list}"
        local_meta_sample_list="${local_sample_list_directory}/${meta_sample_list}"

        dx download --recursive --overwrite "${project_id}:/${meta_sample_list}" --output "${local_meta_sample_list}"
        
        number_of_batches=$( wc -l "${meta_sample_list}" )
    elif [[ ${number_of_batches} != 0 && -n ${meta_sample_list} ]]; then
        echo "Error: a meta sample list must be provided if the number of batches is greater than 1"
        exit 1
    fi

    arg_array=(
        batch_size
        sample_list
        meta_sample_list
        sample_directory
        sample_list_directory
        references_directory
        input_directory
        output_directory
        project_id
        n_jobs
        git_branch
        run_mutect
        mutect_bam_out
        run_haplotypecaller
        run_varscan
        reference_genome
        gnomad_reference
        funcotator_sources
        interval_list
    )

    job_args=$(eval echo "$(printf "%s\n" "${arg_array[@]}" | xargs --replace=% echo '-i%="${%}' | tr '\n' ' ')")
                
    seq 1 "${number_of_batches}" | xargs --replace=% bash -c "submit_job % ${job_args} ${instance_type}"
}
