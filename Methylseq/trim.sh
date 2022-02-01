#!/bin/bash

#script to trim paired reads with error checking code wrapped around it. 
#Is set to autodetect cores available and use them.

while getopts r:R:t:T: flag
do
    case "${flag}" in
        r) read1=${OPTARG};;
        R) read2=${OPTARG};;
        t) read1_trimmed=${OPTARG};;
        T) read2_trimmed=${OPTARG};;
    esac
done


#TO DO: add code to check all inputs are present and if not to give error message with help information about all the script inputs

echo ""
echo "trim command used:"
echo "$0 -r $read1 -R $read2 -t $read1_trimmed -T $read2_trimmed"
echo ""

module load cutadapt/3.4

#if statement to stopcode from running this section if already done

#detect parent dir path of trimmed file 
#check if directory of trimmed file exists yet, if it does go there, else throw error message to fix upstream code.
parent_directory="$(dirname "$read1_trimmed")"

echo "Parent directory for trimming is '$parent_directory'"

if [ -d "$parent_directory" ] 
	then
		#trim code
		if [ ! -f "$read1_trimmed" ] && [ ! -f "$read2_trimmed" ] 
			then
			#I don't think this part of the code needs to cd to the trim directory, so it doesn't. 
			
			echo "$(basename "$read1") and $(basename "$read2") read pair not yet trimmed. Trimming now."
			cutadapt -q 20 -m 15 -a AGATCGGAAGAGC -A AAATCAAAAAAAC -o $read1_trimmed -p $read2_trimmed $read1 $read2 -j 0
			echo "finished trimming $(basename $read1) and $(basename $read2) read pair."

			#optional rsync to copy files back to seq_path directory
			#rsync -vur --exclude "main_genome" --exclude "unmethyl_genome" --exclude "hydroxymethyl_genome" --exclude "methyl_genome" $temp_path/ $seq_path

			else
		        echo "$(basename "$read1") and $(basename "$read2") have already been trimmed"
		fi

	else
		echo "Error: Directory '$parent_directory' does not exist. Directory '$parent_directory' must be created before running '$0'"
fi

#TO DO: This script could use testing code that tests each 'then' and 'else' sections of each 'if' statement! So a total of 4 test cases.

#TO DO: halt code with any error to prevent errors from snowballing and not being detected by end user or by development team. http://web.archive.org/web/20110314180918/http://www.davidpashley.com/articles/writing-robust-shell-scripts.html

#TO DO: #test code: bash trim.sh -r hi -R yo -t dude -T pumpkin
#should output that it is trimming the files. it's looking for the relative path . which exists always. It should trim and output the trimmed files to the current working directory.









 




