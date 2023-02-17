#!/usr/bin/env bash

# Set bash options for verbose output and to fail immediately on errors or if variables are undefined.
set -o xtrace -o nounset -o pipefail -o errexit

function run_job() {
    export PATH="${HOME}/.local/bin:${PATH}"
    curl micro.mamba.pm/install.sh | bash
    micromamba shell init --shell=bash ${HOME}/micromamba
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
        git \
        grep \
        parallel \
        sed \
        tar \
        util-linux
    micromamba create -n gatk gatk4
    conda init bash
    . $HOME/micromamba/etc/profile.d/conda.sh
    mamba init bash
    conda config --set auto_stack 1

    cd "${HOME}"
    git clone --branch refactor1 https://github.com/Lab-Jaiswal/workflows
    batch_size=$(( 2*batch_size ))
    batch_size=$(( 1+batch_size ))

    # Download references
    local_references_directory="${HOME}/${references_directory}"
    dx download --recursive --overwrite "${project_id}:/${references_directory}/" --output "${local_references_directory}"
    
    # Download list of files if provided
    # If a meta list is provided instead, download the file list from the meta list by extracting the row corresponding the array number.
    local_sample_list_directory="${HOME}/${sample_list_directory}"
    mkdir -p "${local_sample_list_directory}"
    local_sample_list="${local_sample_list_directory}/initial_sample_list"

    if [[ -z ${sample_list} ]] && [[ -n ${meta_sample_list} ]]; then
        sample_list=$(sed "${array_number}"q; d "${meta_sample_list}")
    fi
    dx download --recursive --overwrite "${project_id}:/${sample_list}" --output "${local_sample_list}"

    local_input_directory="${HOME}/${input_directory}"
    mkdir -p "${local_input_directory}"
    filtered_list_of_samples="${local_input_directory}/filtered_list_of_samples"

    if [[ ${batch_size} != 0 ]]; then
        head -n ${batch_size} "${local_sample_list}" > "${filtered_list_of_samples}"
    else
        cp "${local_sample_list}" "${filtered_list_of_samples}"
    fi

    xargs --replace=% bash -c "dx download --overwrite ${project_id}:/${sample_directory}/% --output ${local_input_directory}/%" < "${filtered_list_of_samples}"
       
    #cd ~
    #LOGS=~/Outputs/Logs
    #if [ ! -d $LOGS ]; then
             #mkdir -p ${LOGS}
    #fi
    #log_file="${LOGS}/job_${array_number}.log"
        
    "${HOME}"/workflows/BWA_CHIP/submit_BWA_CHIP.sh --n_jobs "${n_jobs}"
    local_output_directory="${HOME}/${output_directory}"
    output_tar="${HOME}/outputs_${array_number}.tar"
    tar --create --file "${output_tar}" "${local_output_directory}"
    outputs_folder_id=$(dx upload "${output_tar}" --brief)
    dx-jobutil-add-output outputs_tar "${outputs_folder_id}" --class=file
}

function main() {
    dx-download-all-inputs 
    
    if [[ ${number_of_batches} == 0 && -n ${meta_sample_list} ]]; then
        local_sample_list_directory="${HOME}/${sample_list_directory}"
        mkdir -p "${local_sample_list_directory}"

        meta_sample_list="${sample_list_directory}/${meta_sample_list}"
        local_meta_sample_list="${local_sample_list_directory}/${meta_sample_list}"

        dx download --recursive --overwrite "${project_id}:/${meta_sample_list}" --output "${local_meta_sample_list}"
        
        number_of_batches=$( wc -l "${meta_sample_list}" )
    elif [[ ${number_of_batches} != 0 && -z ${meta_sample_list} ]]; then
        echo "Error: a meta sample list must be provided if the number of batches is greater than 1"
        exit 1
    fi
                
    for array_number in $(seq 1 "${number_of_batches}"); do
        job_id=$(dx-jobutil-new-job run_job -iarray_number="${array_number}" -ibatch_size="${batch_size}" -isample_list="${sample_list}" -imeta_sample_list="${meta_sample_list}" -isample_directory="${sample_directory}" -isample_list_directory="${sample_list_directory}" -ireference_directory="${references_directory}" -iinput_directory="${input_directory}" -iproject_id="${project_id}" -in_jobs="${n_jobs}" --instance-type mem3_ssd1_v2_x96)
        dx-jobutil-add-output outputs_folder "${job_id}:outputs_tar" --class=jobref --array
    done
}
