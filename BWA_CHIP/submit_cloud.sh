#!/bin/bash
cd ~
mkdir filelists
cd filelists
seq 10 59 | parallel -j8 --progress --ungroup "~/workflows/BWA_CHIP/split_files.sh {}"
find . -type f -not -name "*.*" | xargs -I % mv % %.sh
sed -i '1 i\#!/bin/bash' *.sh

cd ~/workflows/BWA_CHIP

file_list=~/filelists/meta_filelist.txt

find -L "~/filelists" -type f `#list all files in ${fastq_directory}` | \
                             grep ".*\.${sh}$" `#only keep files with FASTQ in name (case insensitive)` | \
                             sed -e 's/\.sh$//g' | \
                             sort -u  `#sort and remove duplicate names` > ${file_list}
