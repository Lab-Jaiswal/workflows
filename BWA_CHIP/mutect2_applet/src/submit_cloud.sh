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
        mamba \
        parallel \
        sed \
        tar \
        util-linux
    micromamba create -n gatk4 gatk4
    conda init bash
    . "${HOME}/micromamba/etc/profile.d/conda.sh"
    mamba init bash
    conda config --set auto_stack 1
    set -o nounset -o xtrace

    git clone --branch "${git_branch}" https://github.com/Lab-Jaiswal/workflows
    
    # Download references
    local_references_directory="${HOME}/${references_directory}"
    dx download --recursive --overwrite "${project_id}:/${references_directory}/" --output "${local_references_directory}"
    
    # Download list of samples
    local_sample_list_directory="${HOME}/${sample_list_directory}"
    mkdir -p "${local_sample_list_directory}"
    local_sample_list="${local_sample_list_directory}/sample_list"
    sample_list=${sample_list_directory}/sample_batches/sample_list_${array_number}
    dx download --recursive --overwrite "${project_id}:/${sample_list}" --output "${local_sample_list}"

    local_input_directory="${HOME}/${input_directory}"
    mkdir -p "${local_input_directory}"
    local_output_directory="${HOME}/${output_directory}"
    mkdir -p "${local_output_directory}"

    sample_directory=$(echo "${sample_directory}" | tr '?' ' ')
    xargs --replace=% bash -c "dx download --overwrite ${project_id}:/\"${sample_directory}/%\" --output ${local_input_directory}/\$(basename %)" < "${local_sample_list}"

    passed_args=(
        bam_extension
        fastq_extension
        assembly
        run_mutect
        run_haplotypecaller
        run_varscan
        mutect_bam_out
        run_funcotator
        run_annovar
        split_intervals
        realign
        varscan_min_coverage
        varscan_min_var_freq
        varscan_max_pvalue
        n_jobs
    )

    reference_args=(
        reference_genome
        interval_list
        sequence_dictionary
        gnomad_reference
        normal_bam
        normal_pileups_table
        funcotator_sources
        transcript_list
        mpileup_interval_bed
        annovarroot
    )

    read -r -a passed_args_array <<< "$(printf "%s\n" "${passed_args[@]}" | \
        xargs --replace=% bash -c 'echo --% "${%}"' | \
        tr '\n' ' ')"
    #read -r -a reference_args_array <<< "$(printf "%s\n" "${reference_args[@]}" | xargs --replace=% bash -c 'echo --% "${local_references_directory}/${%}"' | sed -e 's/.*none$/none/' | tr '\n' ' ')"
    read -r -a reference_args_array <<< "$(printf "%s\n" "${reference_args[@]}" | \
        xargs --replace=% bash -c 'echo --% "${%}"' | \
        sed -e "s? ? ${local_references_directory}/?" | \
        sed -e 's/.*none$/none/' | \
        tr '\n' ' ')"

    "${HOME}/workflows/BWA_CHIP/submit_BWA_CHIP.sh" \
        --input_directory "${local_input_directory}" \
        --output_directory "${local_output_directory}" \
        "${passed_args_array[@]}" \
        "${reference_args_array[@]}"

    output_tar="${HOME}/outputs_${array_number}.tar"
    tar --create --file "${output_tar}" "${local_output_directory}"
    outputs_folder_id=$(dx upload "${output_tar}" --brief)
    dx-jobutil-add-output outputs_tar "${outputs_folder_id}" --class=file
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
        findutils \
        grep \
        mamba
    conda init bash
    . "${HOME}/micromamba/etc/profile.d/conda.sh"
    mamba init bash
    conda config --set auto_stack 1
    set -o nounset -o xtrace

    # Fix me later
    #function submit_job() {
        #set -o pipefail -o errexit -o nounset -o xtrace
        #array_number=${1}
        #job_args=${2}
        #instance_type=${3}
        #read -r -a job_args_array <<< "${job_args}"

        #job_id=$(dx-jobutil-new-job run_job -iarray_number="${array_number}" "${job_args_array[@]}" --instance-type="${instance_type}")
        #dx-jobutil-add-output outputs_folder "${job_id}:outputs_tar" --class=jobref --array
    #}
    #export -f submit_job
    
    batch_size=$(( 2*batch_size ))

    local_sample_list_directory="${HOME}/${sample_list_directory}"
    mkdir -p "${local_sample_list_directory}"
    all_samples_list="${local_sample_list_directory}/all_samples"

    # Get list of folders in main CRAM folder, then get the list of files in each folder.
    # dx tree is too slow to iterate over all the files, so we prepend the folder name to each file name.
    # The sed command cannot use '/' as the substitution separator because the directory paths contain '/', so '?' is used instead.
    # This output must be saved to a file, because the subsequent xargs command cannot parse lines that are too long.
    # We only run this step if the file "all_samples" doesn't already in the sample_lists directory
    #if ! dx ls "${project_id}":"${sample_list_directory}" | cut -d ' ' -f 1 | grep -q "^all_samples$"; then
        #dx upload --path "${project_id}":"${sample_list_directory}/" "${all_samples_list}"
    #else
        #dx download "${project_id}":"${sample_list_directory}/all_samples" --output "${all_samples_list}"
    #fi

    dx ls "${project_id}":"${sample_directory}" | xargs --replace=% bash -c "dx ls ${project_id}:\"${sample_directory}\"/% | sed -e 's?^?%?g'" > "${all_samples_list}"

    filtered_sample_list="${local_sample_list_directory}/filtered_sample_list"
    if [[ ${sample_list} != "none" ]]; then
        local_sample_list="${local_sample_list_directory}/${sample_list}"
        dx download "${project_id}":"${sample_list_directory}/${sample_list}" --output "${local_sample_list}"
        grep -f "${local_sample_list}" < "${all_samples_list}" > "${filtered_sample_list}"
    else
        cp "${all_samples_list}" "${filtered_sample_list}"
    fi

    filtered_samples_batches="${local_sample_list_directory}/filtered_samples_batches"
    xargs --max-args=${batch_size} < "${filtered_sample_list}" > "${filtered_samples_batches}"
    filtered_samples_batches_length=$(wc -l "${filtered_sample_list}" | cut -d ' ' -f 1)

    # Check if the number of bat
    if [[ ${number_of_batches} == 0 ]]; then
        number_of_batches="${filtered_samples_batches_length}"
    elif (( number_of_batches > filtered_samples_batches_length )); then
        filtered_sample_list_length=$(wc -l "${filtered_sample_list}")
        echo "Error: you have requested ${number_of_batches} but the largest number of batches possible with ${filtered_sample_list_length} samples and a batch size of ${batch_size} is ${filtered_samples_batches_length}."
        echo "Please reduce the number of batches or batch size until they are compatible with each other and the number of samples you wish to include."
        exit 1
    fi

    sample_batches_directory="${local_sample_list_directory}/sample_batches"
    mkdir -p "${sample_batches_directory}"
    seq 1 "${number_of_batches}" | xargs --replace=% bash -c "sed '%q; d' ${filtered_samples_batches} | tr ' ' '\n' > ${sample_batches_directory}/sample_list_%"
    dx rm --recursive --force --all "${project_id}:${sample_list_directory}/sample_batches"
    dx upload --recursive --path "${project_id}:${sample_list_directory}/sample_batches" "${sample_batches_directory}" 

    # HACK: The CRAM folder contains spaces in the file names.  This causes numerous problems, including that it breaks bash array logic.
    # To circumvent this, we replace the spaces with question marks, then fix these later in the download step
    sample_directory=$(echo "${sample_directory}" | tr ' ' '?')
    
    arg_array=(
        sample_directory
        sample_list_directory
        references_directory
        input_directory
        output_directory
        project_id
        bam_extension
        fastq_extension
        assembly
        git_branch
        run_mutect
        run_haplotypecaller
        run_varscan
        mutect_bam_out
        run_funcotator
        run_annovar
        split_intervals
        realign
        n_jobs
        varscan_min_coverage
        varscan_min_var_freq
        varscan_max_pvalue
        reference_genome
        interval_list
        sequence_dictionary
        gnomad_reference
        normal_bam
        normal_pileups_table
        funcotator_sources
        transcript_list
        mpileup_interval_bed
        annovarroot
    )

    read -r -a job_args <<< "$(printf "%s\n" "${arg_array[@]}" | xargs --replace=% bash -c 'echo -i%="${%}"' | tr '\n' ' ')"
                
    #seq 1 "${number_of_batches}" | xargs --replace=% bash -c "submit_job % ${instance_type} ${job_args[@]}"

    for array_number in $(seq 1 "${number_of_batches}"); do
        job_id=$(dx-jobutil-new-job run_job -iarray_number="${array_number}" "${job_args[@]}" --instance-type="${instance_type}")
        dx-jobutil-add-output outputs_folder "${job_id}:outputs_tar" --class=jobref --array
    done
}
