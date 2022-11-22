#!/bin/bash

function main() {

        array_number=${1}
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

        for i in $(seq 1 ${array_number}); do

                ARRAY_PREFIX="$(sed "${i}q; d" "${meta_file_list}")"
                echo "number: $i"
                ARRAY_PREFIX="$(basename ${ARRAY_PREFIX})"

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
                #cd ~/workflows/BWA_CHIP
                echo "submitting BWA_CHIP"


               docker run --rm -v /home/dnanexus:/home/dnanexus \
               -v /usr/local/bin/dx:/usr/local/bin/dx \
               -v /usr/local/bin/dx-jobutil-add-output:/usr/local/bin/dx-jobutil-add-output \
               -v /usr/local/bin/:/usr/local/bin \
               -v /usr/local/reference/:/usr/local/reference/ \
               -w /home/dnanexus us.gcr.io/broad-gatk/gatk:4.1.6.0 \ /bin/bash /usr/local/bin/submit_BWA_CHIP.sh --working_dir ~ --file_extension cram --mutect --container_engine docker --mode cloud --array_prefix ${list_of_samples} &>${log_file}

        done
         #filtered_vcf_output=$(dx upload "${filtered_vcf[i]}" --brief)
         #dx-jobutil-add-output filtered_vcf_output "${filtered_vcf_output}" --class=file --array

./submit_BWA_CHIP.sh --working_dir ~ --file_extension cram --mutect --container_engine docker --mode cloud --array_prefix ${list_of_samples} &>${log_file}
}

#main $1
