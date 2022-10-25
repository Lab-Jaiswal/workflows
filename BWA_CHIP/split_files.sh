#!/bin/bash
#run using seq 10 59 | parallel -j8 --progress --ungroup "./split_files.sh {}"
array_length=$1
dx ls project-G5B07V8JPkg740v9GjfF9PzV:/Bulk/Exome\ sequences/Exome\ OQFE\ CRAM\ files/${array_length}/ > folder_${array_length}.txt

split -l 1000 -d folder_${array_length}.txt folder_${array_length}_split
