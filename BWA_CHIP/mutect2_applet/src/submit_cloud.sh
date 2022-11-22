#!/bin/bash
function run_job() {
        array_number=${1}
        batch_size=${2}
        batch_size=$(( 2*batch_size ))
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
        echo "number: $array_number"
        ARRAY_PREFIX="$(basename ${ARRAY_PREFIX})"

        list_of_samples=~/file_lists/${ARRAY_PREFIX}
        if [ ! -f ${list_of_samples} ]; then
                 cd ${File_Lists}
                 dx download -r "project-G5B07V8JPkg740v9GjfF9PzV:/File_Lists/${ARRAY_PREFIX}"
                 filtered_list_of_samples=~/file_lists/filtered_list_of_samples
                 head -n $batch_size $list_of_samples > $filtered_list_of_samples
        fi
        cd ~
        LOGS=~/Outputs/Logs
        if [ ! -p $LOGS ]; then
                 mkdir -p ${LOGS}
        fi
      
       ./submit_BWA_CHIP.sh --working_dir ~ --file_extension cram --mutect --container_engine docker --mode cloud --array_prefix ${filtered_list_of_samples} &>${log_file}
       
       
       Outputs_folder=$(dx upload ~/Outputs --brief)
       dx-jobutil-add-output Outputs_folder "${Outputs_folder}" --class=file --array
}

function main() {

        number_of_batches=${1}
        batch_size=${2}
        
        for i in $(seq 1 ${array_number}); do
        
                echo "submitting BWA_CHIP"
                
               dx-jobutil-new-job run_job -i ${list_of_samples} -i ${batch_size}

        done

}

#main $1
