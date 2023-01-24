#!/bin/bash
set -u

function run_job() {
        cd $HOME
        git clone https://github.com/Lab-Jaiswal/workflows
        #git checkout build_applet
        #next two were commented out
        #curl -L -O https://github.com/Lab-Jaiswal/workflows/archive/build_applet.tar.gz
        #tar xzf build_applet.tar.gz 
        #array_number=${1}
        #batch_size=${2}
        batch_size=$(( 2*batch_size ))
        batch_size=$(( 1+batch_size ))
        
        #dx download -r project-G5B07V8JPkg740v9GjfF9PzV:/workflows
        #if [ ! -d ~/workflows/BWA_CHIP ]; then
            #echo "it was not downloaded"
            #exit 1
        #fi
        
        #if [ ! -f ~/workflows/BWA_CHIP/submit_BWA_CHIP.sh ]; then
            #echo "submit_BWA_CHIP.sh was not downloaded"
            #exit 1
        #fi


        File_Lists=~/file_lists
        if [ ! -d ${File_Lists} ]; then
                mkdir -p ${File_Lists}
                cd ${File_Lists}
        fi
        
        #meta_file_list=~/file_lists/meta_filelist.txt
        meta_file_list=$file_list

        #if [ ! -f ${meta_file_list} ]; then
                #dx download -r project-G5B07V8JPkg740v9GjfF9PzV:/File_Lists/meta_filelist.txt
                dx download -r "$file_list" -f
		echo "just downloaded file_list"
                file_list_name=$(ls -Art | tail -n 1)
                if grep -q sh "$file_list_name"; then
                        meta_list=true
                else
                        meta_list=false
                fi
        #fi
        
        if [[ meta_list == true ]]; then
                ARRAY_PREFIX="$(sed "${array_number}q; d" "${meta_file_list}")"
                echo "number: $array_number"
                ARRAY_PREFIX="$(basename ${ARRAY_PREFIX})"

                list_of_samples=~/file_lists/${ARRAY_PREFIX}
                if [ ! -f ${list_of_samples} ]; then
                         cd ${File_Lists}
			 echo "about to download"
                         dx download -r "project-G5B07V8JPkg740v9GjfF9PzV:/File_Lists/${ARRAY_PREFIX}"
                         echo "just downloaded"
			 ls
			 mkdir -p ~/Inputs
                        filtered_list_of_samples=~/Inputs/${ARRAY_PREFIX}
			
                         if [[ $batch_size != 0 ]]; then
                                head -n $batch_size $list_of_samples > $filtered_list_of_samples
                         else
                                cp $list_of_samples $filtered_list_of_samples
                         fi

                fi
        else
                mkdir -p ~/Inputs
                filtered_list_of_samples=~/Inputs/$file_list_name
                if [[ $batch_size != 0 ]]; then
			ls
			real_path=$(realpath .)
			echo "real_path: $real_path"
			head -n $batch_size $file_list_name > $filtered_list_of_samples
                else
                        cp $file_list_name $filtered_list_of_samples
                fi
                                
        fi
       
        cd ~
        LOGS=~/Outputs/Logs
        if [ ! -d $LOGS ]; then
                 mkdir -p ${LOGS}
        fi
        log_file="${LOGS}/job_${array_number}.log"
        
        echo "command to run:" >> $log_file
        echo "bash $HOME/workflows/BWA_CHIP/submit_BWA_CHIP.sh --working_dir $HOME --file_extension cram --mutect --container_engine docker --mode cloud --array_prefix '${filtered_list_of_samples}'" >> $log_file
      #commented out below because not needed for test
       #cd ~/workflows/BWA_CHIP
       bash $HOME/workflows/BWA_CHIP/submit_BWA_CHIP.sh --working_dir $HOME --file_extension cram --mutect --container_engine docker --mode cloud --array_prefix "${filtered_list_of_samples}" --n_jobs 96 #&>${log_file}
       cd $HOME
       echo "$HOME" >> $log_file
       echo "filtered_list_of_samples: $filtered_list_of_samples" >> $log_file
       Output_Dir=~/Outputs
       #below line is a test to see if the file lists are even downloading
       #Output_Dir=~/workflows
       echo "Output Dir: $Output_Dir" >> $log_file
       Output_tar=Outputs_${array_number}.tar
       tar cf $Output_tar $Output_Dir
       Outputs_folder_id=$(dx upload $Output_tar --brief)
       dx-jobutil-add-output Outputs_tar "${Outputs_folder_id}" --class=file
}

function main() {
        dx-download-all-inputs 
        
        
        if [[ ${number_of_batches} = 0 ]]; then
                File_Lists=~/file_lists
                if [ ! -d ${File_Lists} ]; then
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
                
               job_id=$(dx-jobutil-new-job run_job -iarray_number=${i} -ibatch_size=${batch_size} -ifile_list=${file_list} --instance-type mem3_ssd1_v2_x96)
               echo "job_id: $job_id"
               dx-jobutil-add-output Outputs_folder "${job_id}:Outputs_tar" --class=jobref --array

        done

}

#main $1
