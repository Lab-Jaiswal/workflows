#!/bin/bash
echo ""
echo "entering extract_methyl script"
echo ""

bismark_extraction_input=$1
bismark_extraction_report=$2
output_temp_directory=$3
output_directory=$4
genome_fasta_path=$5
cores=$6
parameter_file=$7

echo "extract_methyl command used the following parameters:
$0 $1 $2 $3 $4 $5 $6 $7"

if [ $SLURM_ARRAY_TASK_ID -eq 1 ]; then
        echo "bismark_extraction_input=$1
        bismark_extraction_report=$2
        output_temp_directory=$3
        output_directory=$4
        genome_fasta_path=$5
        cores=$6
        parameter_file=$7
        " >> $parameter_file
fi

module load bismark/0.22.3

#copy over output files if they are already made
bismark_extraction_report_name=$(basename "${index_output}")
bismark_extraction_input_name=$(basename "${index_input}")


N_cores=$((cores/3))

#if there is not the bismark extraction output in the output directory, then transfer files needed to make it and make it. else skip
if [ ! -f "$output_directory/$bismark_extraction_report_name" ]; then #TO DO: should add requirement to also have the extracted methylation files too! Not just the report.
    #transfer files needed if not present yet in temp directory.
    if [ ! -f "$bismark_extraction_input" ]; then
        rsync -vur --include="${bismark_extraction_input_name}" --exclude="*" "$output_directory/" $output_temp_directory
    fi

        echo "$bismark_extraction_report does not exist yet"
    bismark_methylation_extractor --gzip --cytosine_report --bedGraph --genome_folder "$genome_fasta_path" $bismark_extraction_input -o $output_temp_directory --multicore $N_cores --paired-end
        echo ""
        echo "extract_methylation now complete for this genome"
        echo ""

        rsync -vur $output_temp_directory/ $output_directory

else
    echo ""
    echo "methylation extraction for this genome already created."
    echo ""
fi


