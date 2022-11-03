#!/bin/bash

output_directory=$1
mkdir -p ${output_directory}
code_directory=$(realpath .)

cd ${output_directory}
        seq 10 59 | parallel -j8 --progress --ungroup "${code_directory}/split_files.sh {}"
        find . -type f -not -name "*.*" | xargs -I % mv % %.sh
        sed -i '1 i\#!/bin/bash' *.sh

file_list="${output_directory}/meta_filelist.txt"
find -L ${output_directory} -type f `#list all files in ${fastq_directory}` | \
        grep ".*\.sh$" `#only keep files with FASTQ in name (case insensitive)` | \
        sort -u  `#sort and remove duplicate names` > ${file_list}

output_dir=c=${output_directory%/} # Remove possible trailing
dx upload -r $output_dir

#dx mv "*folder*" File_Lists/
