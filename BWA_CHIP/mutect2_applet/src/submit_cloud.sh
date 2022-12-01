#!/bin/bash
set -u

function run_job() {
        cd ~
        git clone https://github.com/Lab-Jaiswal/workflows
        git checkout build_applet
        
        array_number=${1}
        batch_size=${2}
        batch_size=$(( 2*batch_size ))
        
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
                 if [[ $batch_size != 0 ]]; then
                        filtered_list_of_samples=~/file_lists/filtered_list_of_samples
                        head -n $batch_size $list_of_samples > $filtered_list_of_samples
                 else
                        filtered_list_of_samples=$list_of_samples
                 fi
                 
        fi
       
        cd ~
        LOGS=~/Outputs/Logs
        if [ ! -p $LOGS ]; then
                 mkdir -p ${LOGS}
        fi
      
       bash ~/workflows/submit_BWA_CHIP.sh --working_dir ~ --file_extension cram --mutect --container_engine docker --mode cloud --array_prefix ${filtered_list_of_samples} #&>${log_file}
       Output_Dir=~/Outputs
       echo "Output Dir: $Output_Dir"
       Output_tar=Outputs_${array_number}.tar
       tar cf $Output_tar $Output_Dir
       Outputs_folder_id=$(dx upload $Output_tar --brief)
       dx-jobutil-add-output Outputs_folder "${Outputs_folder_id}" --class=file
}

function main() {
        dx-download-all-inputs 
        
        
        if [[ ${number_of_batches} = 0 ]]; then
                File_Lists=~/file_lists
                if [ ! -p ${File_Lists} ]; then
                        mkdir -p ${File_Lists}
                        cd ${File_Lists}
                fi
        
                meta_file_list=~/file_lists/meta_filelist.txt

                if [ ! -f ${meta_file_list} ]; then
                        dx download -r project-G5B07V8JPkg740v9GjfF9PzV:/File_Lists/meta_filelist.txt
                fi
                
                number_of_batches=$( wc -l $meta_file_list )
         fi
                
        for i in $(seq 1 ${number_of_batches}); do
        
                echo "submitting BWA_CHIP"
                
               dx-jobutil-new-job run_job -i${i} -i${batch_size}

        done

}

#main $1
