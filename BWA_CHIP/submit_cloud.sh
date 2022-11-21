#!/bin/bash

function main() {

        array_number=$1
        cd ~
        File_Lists=~/file_lists
        if [ ! -p ${File_Lists} ]; then
                mkdir -p ${File_Lists}
                cd ${File_Lists}
        fi

        meta_file_list=~/file_lists/meta_filelist.txt

        if [ ! -f ${meta_file_list} ]; then
                dx download -r project-G5B07V8JPkg740v9GjfF9PzV:/File_Lists/meta_filelist.txt
        fi

        ARRAY_PREFIX="$(sed "${array_number}q; d" "${meta_file_list}")"
        ARRAY_PREFIX="$(basename $ARRAY_PREFIX)"

        list_of_samples=~/file_lists/${ARRAY_PREFIX}
        if [ ! -f ${list_of_samples} ]; then
                cd ${File_Lists}
                dx download -r "project-G5B07V8JPkg740v9GjfF9PzV:/File_Lists/${ARRAY_PREFIX}"
        fi
        cd ~
        LOGS=~/Outputs/Logs
        if [ ! -p $LOGS ]; then
               mkdir -p ${LOGS}
        fi

        log_file="${LOGS}/log_file_${array_number}.sh"
        cd ~/workflows/BWA_CHIP
        echo "submitting BWA_CHIP"
        ./submit_BWA_CHIP.sh --working_dir ~ --file_extension cram --mutect --container_engine docker --mode cloud --array_prefix ${list_of_samples} &>${log_file}
}

main 1
