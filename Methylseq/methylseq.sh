#!/bin/bash

#SBATCH --time=72:00:00
#SBATCH --account=sjaiswal
#SBATCH --cpus-per-task=20
#SBATCH --mem=256GB
#SBATCH --job-name=methylseq

##################################################################################################################################
#############################################---STEP 1: SET UP ARGUMENTS---####################################################### 
##################################################################################################################################
data_path=$1
output_path=$2
genetic_locations=$3
cores=$4
log_name=$5
parameter_file=$6
code_directory=$7
Logs=$8
parameter_file=$9
initial_path=${10}

line_number=$SLURM_ARRAY_TASK_ID #get index of which file to process from $SLURM_ARRAY_TASK_ID provided by SLURM
if [ $SLURM_ARRAY_TASK_ID -eq 1 ]; then
    most_recent=$(ls $Logs -c | head -n 1 | sed 's/.*\///' | cut -d'_' -f1 | sed 's/[^0-9]*//g')
    if [ $most_recent -gt $SLURM_JOBID ]; then
        echo "log file name: ${log_name}_${most_recent}_#.log" 
        >> $parameter_file
    else
        echo "log file name ${log_name}_${SLURM_JOB_ID}_#.log"
        >> $parameter_file
    fi
fi 

fastq_file="${data_path}/fastq/FASTQs"                                                         #provide path to file containing list of fastq files
fastq_path="$(sed "${line_number}q; d" $fastq_file)"                                           #extract only the line number corresponding to $SLURM_ARRAY_TASK_ID
sample_name="${fastq_path##*/}"
fastq_temp=$(basename "${fastq_path}")

line_count=$( wc -l < "${genetic_locations}" )
total_genomes=$(bc -l <<< "scale=0; (($line_count / 3) - 1)")                                   #the total number of genomes requested to be mapped against

temp_path=$(mktemp -d)
echo "temp_path is: " $temp_path
echo "copying FASTQs from the data path..."                                                     #copy data from data_path to the temp_path
rsync -vur "$data_path/fastq/" $temp_path

cd $temp_path

##################################################################################################################################
###########################################---STEP 2: LOAD MODULES---######################################################
##################################################################################################################################
module load bismark/0.22.3
module load samtools
echo "Bismark and samtools modules have been loaded"

##################################################################################################################################
###############################################---STEP 3: RUN trim.sh---##########################################################
##########################################---Copy output data (if any exist)---###################################################
##################################################################################################################################
read1="${temp_path}/${sample_name}_R1_001.fastq.gz"
read2="${temp_path}/${sample_name}_R2_001.fastq.gz"
read1_trimmed=$(echo $read1 | sed 's/fastq.gz/trimmed.fastq.gz/')
read2_trimmed=$(echo $read2 | sed 's/fastq.gz/trimmed.fastq.gz/')

echo "running trim.sh"
${code_directory}/trim.sh $read1 $read2 $read1_trimmed $read2_trimmed $data_path $parameter_file

echo "copying data from the output file (if there is any)"                                      #copy data from output_path to temp_path
rsync -vur --exclude "Logs" --exclude "Parameters" $initial_path/ $temp_path                    #done after trim.sh to simplify rsync step

##################################################################################################################################
#########################################---STEP 5: RUN map_and_deduplicate.sh---#################################################
################################################---sort_and_index.sh---###########################################################
################################################---extract_methyl.sh---###########################################################
##################################################################################################################################
for i in $(seq 0 $total_genomes); do
    number1=$(bc -l <<< "scale=0; (($i * 3) +1)")
    number2=$(bc -l <<< "scale=0; (($i * 3) +2)")
    number3=$(bc -l <<< "scale=0; (($i * 3) +3)") 
            
    genome_name=$(sed -n ${number1}'p' $genetic_locations)
    genome_fasta_path=$(sed -n ${number2}'p' $genetic_locations)
    deduplicate=$(sed -n ${number3}'p' $genetic_locations)

    k=$(bc -l <<< "scale=0; ($i + 1)")
    
    if [ $SLURM_ARRAY_TASK_ID -eq 1 ]; then
        echo "LOOP $k INFORMATION:
                genome name: $genome_name
                genome fasta path: $genome_fasta_path
                deduplicate status: $deduplicate
                " >> $parameter_file
    fi

    echo "LOOP $k INFORMATION:
        genome name: $genome_name
        genome fasta path: $genome_fasta_path
        deduplicate status: $deduplicate "            
       
    output_temp_directory=$(find $temp_path -type d -name "$genome_name")
    output_directory=$(find $initial_path -type d -name "$genome_name")
            
    j=$(bc -l <<< "scale=0; ($i - 1)")
    if [ $j -gt -1 ]; then            
        number4=$(bc -l <<< "scale=0; (($j * 3) +1)")
        genome_name_input=$(sed -n ${number4}'p' $genetic_locations)
        input_directory=$(find $initial_path -type d -name "$genome_name_input")
        input_temp_directory=$(find $temp_path -type d -name "$genome_name_input")
    else
        input_temp_directory=$temp_path
    fi
    
    echo "copying transcriptome to the temporary file directory"
    rsync -vur "$genome_fasta_path/" "$temp_path/${genome_name}_fasta"
    temp_genome="$temp_path/${genome_name}_fasta"
        
    read1_addition="_unmapped_reads_1.fq.gz"
    read2_addition="_unmapped_reads_2.fq.gz"
    bismark_addition="_unmapped_reads_1"

    if [ $i -eq 0 ]; then
        echo "level is 0"
        read1_filename="${sample_name}_R1_001.trimmed.fastq.gz"
        read2_filename="${sample_name}_R2_001.trimmed.fastq.gz"
        read1_trimmed="${sample_name}_R1_001.trimmed"
        read2_trimmed="${sample_name}_R2_001.trimmed"
        bismark_output="${output_temp_directory}/${read1_trimmed}_bismark_bt2_PE_report.txt"
    else
        read1_ending=${read1_ending}${read1_addition}
        read2_ending=${read2_ending}${read2_addition}
        read1_trimmed="${sample_name}_R1_001.trimmed.fastq.gz"
        read2_trimmed="${sample_name}_R2_001.trimmed.fastq.gz"
        read1_filename="${read1_trimmed}${read1_ending}"
        read2_filename="${read2_trimmed}${read2_ending}"
        bismark_file=$(echo $read1_filename | sed 's/\(.*\).fq.gz/\1/')
        bismark_filename="${bismark_file}_bismark_bt2_PE_report.txt"
        bismark_output="${output_temp_directory}/${bismark_filename}"
    fi
       
    read1_input="${input_temp_directory}/${read1_filename}"
    read2_input="${input_temp_directory}/${read2_filename}"
        
    dedup_input=$(echo $bismark_output | sed 's/PE_report.txt/pe.bam/')
    dedup_output=$(echo $bismark_output | sed 's/PE_report.txt/pe.deduplicated.bam/')
    
    ${code_directory}/map_and_deduplicate.sh $read1_input $read2_input $bismark_output $dedup_input $dedup_output $temp_genome $output_temp_directory $output_directory $cores $deduplicate $parameter_file

    cd ${code_directory}
    
    if [ $deduplicate == TRUE ] || [ "$deduplicate" == "true" ] || [ "$deduplicate" == "TRUE" ]; then
            sort_input=$(echo $bismark_output | sed 's/PE_report.txt/pe.deduplicated.bam/')
            index_input="${dedup_output}.sorted.bam"
            index_output="${index_input}.bai"
    else
            sort_input=$(echo $bismark_output | sed 's/PE_report.txt/pe.bam/')
            index_input="${sort_input}.sorted.bam"
            index_output="${index_input}.bai"
    fi

    ${code_directory}/sort_and_index.sh $sort_input $index_input $index_output $output_temp_directory $output_directory $parameter_file

    bismark_input=$index_input
    bismark_methyl_output=$(echo $bismark_input | sed 's/\(.*\).bam/\1_splitting_report.txt/')
        
    ${code_directory}/extract_methyl.sh $bismark_input $bismark_methyl_output $output_temp_directory $output_directory $temp_genome $cores $parameter_file

done 

##################################################################################################################################
##################################---STEP 10: PREVIOUSLY insert_size_analysis.sh---###############################################
##################################################################################################################################
picard_output="${dedup_input}_picard_insert_size_plot.pdf"
if [ ! -f $picard_output ]; then 
    module load R
    module load picard/2.9.5
    picard CollectInsertSizeMetrics INPUT=$dedup_input OUTPUT=$dedup_input\_picard_insert_size_metrics.txt HISTOGRAM_FILE=$dedup_input\_picard_insert_size_plot.pdf METRIC_ACCUMULATION_LEVEL=ALL_READS
        echo "picard insert_size_analysis complete"
    rsync -vur $output_temp_directory/ $output_directory
else
    echo "picard insert size analysis already complete"
fi
