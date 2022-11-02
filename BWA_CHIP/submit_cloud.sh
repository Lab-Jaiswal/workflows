#!/bin/bash

main() {
#like on batch node
array_number=$1

##################################################################################################################################
##########################################---STEP 1: GET NECESSARY CODE/ TOOLS---################################################# 
##################################################################################################################################

cd ~
git clone https://github.com/Lab-Jaiswal/workflows
sudo apt-get update
sudo apt-get install -y parallel

##################################################################################################################################
##############################################---STEP 2: DOWNLOAD FILES---######################################################## 
##################################################################################################################################

#list of files
File_Lists=~/file_lists
mkdir -p ${File_Lists}
cd ${File_Lists}
dx download -r project-G5B07V8JPkg740v9GjfF9PzV:/File_Lists/meta_filelist.txt
ARRAY_PREFIX="$(sed "${array_number}q; d" "${file_list}")" 
ARRAY_PREFIX="$(basename $ARRAY_PREFIX)"
dx download -r "project-G5B07V8JPkg740v9GjfF9PzV:/File_Lists/${ARRAY_PREFIX}"
cd ..

#References
dx download -r project-G5B07V8JPkg740v9GjfF9PzV:/References/
cd ~/References/cloud_references
mv GRCh* ../
mv 201* ../
cd ../..

#Inputs
cd ${File_Lists}
Folder_Number=$(echo $ARRAY_PREFIX | grep -oP '(?<=_).*(?=_)')
sed -e "1 ! s@^@dx\ download\ project-G5B07V8JPkg740v9GjfF9PzV:/Bulk/Exome\\\ sequences/Exome\\\ OQFE\\\ CRAM\\\ files/${Folder_Number}@" ${ARRAY_PREFIX} > download_file.sh



















cd ~/workflows/BWA_CHIP

#dx_download meta_filelist
file_list="${FILELIST}/meta_filelist.txt"

#add dx download to download the folder with the lists
#pregenerated and on dnanexus
FILELIST=~/filelists
mkdir -p ${FILELIST}
#pass Array prefix to submit bwa chip so it can run with a list of samples
#it will need to download within bwa chip

#this function is what is passed to dx-jobutil-new-job
#install vim (if run in one go) if it is not in there and run this script line by line interactively to see if it works. 
}
#like on login node
#dx download list of files in applet command
#we call dx-jobutil-new-job in this script 500 times (once for every list in list of list using xargs)
#the applet calls this script (NOT THE FUNCTION)

if [[ -z "$(ls -A ${FILELIST})" ]]; then
        cd ${FILELIST}
        seq 10 59 | parallel -j8 --progress --ungroup "~/workflows/BWA_CHIP/split_files.sh {}"
        find . -type f -not -name "*.*" | xargs -I % mv % %.sh
        sed -i '1 i\#!/bin/bash' *.sh
fi

find -L ${FILELIST} -type f `#list all files in ${fastq_directory}` | \
 grep ".*\.sh$" `#only keep files with FASTQ in name (case insensitive)` | \
                             #sed -e 's/\.sh$//g' | \
                             sort -u  `#sort and remove duplicate names` > ${file_list}
