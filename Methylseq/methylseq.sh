#!/bin/bash

#SBATCH --time=32:00:00
#SBATCH --account=sjaiswal
#SBATCH --cpus-per-task=20
#SBATCH --mem=256GB
#SBATCH --job-name=methylseq

##################################################################################################################################
#############################################---STEP 1: SET UP ARGUMENTS---####################################################### 
##################################################################################################################################
data_path=$1
output_path=$2
seq_path="$1/fastq"
unmethyl_control_fasta=$3
chmod 775 $unmethyl_control_fasta
unmethyl_control=$4
hydroxymethyl_control_fasta=$5
chmod 775 $hydroxymethyl_control_fasta
hydroxymethyl_control=$6
methyl_control_fasta=$7
chmod 755 $methyl_control_fasta
methyl_control=$8
main_genome=$9
chmod 775 $main_genome
phix_path=${10}
cores="${11}"
log_name="${12}"
parameter_file="${13}"

line_number=$SLURM_ARRAY_TASK_ID #get index of which file to process from $SLURM_ARRAY_TASK_ID provided by SLURM

if [ $SLURM_TASK_ARRAY_ID -eq 1 ]; then
    most_recent=$(ls $output_path/Logs -c | head -n 1 | sed 's/.*\///' | cut -d'_' -f1 | sed 's/[^0-9]*//g')
    if [ $most_recent -gt $SLURM_JOBID ]; then
        echo "log file name: ${log_name}_${most_recent}_#.log" 
        >> $parameter_file
    else
        echo "log file name ${log_name}_${SLURM_JOB_ID}_#.log"
        >> $parameter_file
    fi
fi


##################################################################################################################################
############################################---STEP 2: COPY TO TEMP PATH---####################################################### 
##################################################################################################################################
fastq_file="${data_path}/fastq/FASTQs"                                                         #provide path to file containing list of fastq files
fastq_path="$(sed "${line_number}q; d" $fastq_file)"                                           #extract only the line number corresponding to $SLURM_ARRAY_TASK_ID

sample_prefix=$(ls $fastq_path | grep fastq | head -n 1 | sed -e 's/_S[0-9]*_L[0-9]*_[IR][0-9]_[0-9]*.fastq.gz//g')
sample_name="${fastq_path##*/}"
filename=$(basename "${fastq_path}")
fastq_temp="${filename}"
PREFIX=$filename

seq_path="$data_path/fastq"

echo "copying FASTQs..."
temp_path=$(mktemp -d /tmp/tmp.XXXXXXXXXX)
echo "temp_path is: " $temp_path
#copy fastq files to temp_path
rsync -vur "$data_path/fastq/" $temp_path
rsync -a --relative "$output_path/./$unmethyl_control/$hydroxymethyl_control/$methyl_control/$genome_alignment/" "$temp_path/"
echo "FASTQs have been copied to the temporary file"

echo "copying reference transcriptome..."
rsync -vur "$main_genome/" "$temp_path/main_genome"
rsync -vur "$unmethyl_control_fasta/" "$temp_path/unmethyl_genome"
rsync -vur "$hydroxymethyl_control_fasta/" "$temp_path/hydroxymethyl_genome"
rsync -vur "$methyl_control_fasta/" "$temp_path/methyl_genome"

echo "Transcriptomes have been copied to the temporary file directory"


##################################################################################################################################
###########################################---STEP 3: CREATE SAMPLE NAMES---######################################################
##################################################################################################################################
module load bismark/0.22.3
module load samtools
"Bismark and samtools modules have been loaded"

R1="${fastq_temp}_R1_001.fastq.gz"
R2="${fastq_temp}_R2_001.fastq.gz"

R1_trimmed="${fastq_temp}_R1_001.trimmed.fastq.gz"
R2_trimmed="${fastq_temp}_R2_001.trimmed.fastq.gz"

trimmed_R1_file_path=$(echo $R1_trimmed| sed 's/.fastq.gz//')
trimmed_R2_file_path=$(echo $R2_trimmed| sed 's/.fastq.gz//')

trimmed_R1_file_name="${trimmed_R1_file_path##*/}"
trimmed_R2_file_name="${trimmed_R2_file_path##*/}"
R1_trimmed_bismark_PE_report="${trimmed_R1_file_name}_bismark_bt2_PE_report.txt"

R1_unmethyl="${temp_path}/$unmethyl_control/${trimmed_R1_file_name}.fastq.gz_unmapped_reads_1.fq.gz"
R2_unmethyl="${temp_path}/$unmethyl_control/${trimmed_R2_file_name}.fastq.gz_unmapped_reads_2.fq.gz"

R1_hydroxy="${temp_path}/$unmethyl_control/$hydroxymethyl_control/${trimmed_R1_file_name}.fastq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1.fq.gz"
R2_hydroxy="${temp_path}/$unmethyl_control/$hydroxymethyl_control/${trimmed_R2_file_name}.fastq.gz_unmapped_reads_2.fq.gz_unmapped_reads_2.fq.gz"

trimmed_R1_bismark_bam="${trimmed_R1_file_path}_bismark_bt2_pe.bam"

unmethyl_R1_file_path="${temp_path}/$unmethyl_control/$hydroxymethyl_control/${trimmed_R1_file_name}.fastq.gz_unmapped_reads_1"
unmethyl_R1_file_name="${unmethyl_R1_file_path##*/}"
R1_unmethyl_bismark_PE_report="${unmethyl_R1_file_name}_bismark_bt2_PE_report.txt"

unmethyl_R1_bismark_bam="${unmethyl_R1_file_path}_bismark_bt2_pe.bam"

hydroxy_R1_file_path="${temp_path}/$unmethyl_control/$hydroxymethyl_control/$methyl_control/${trimmed_R1_file_name}.fastq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1"
hydroxy_R1_file_name="${hydroxy_R1_file_path##*/}"
R1_hydroxy_bismark_PE_report2="${hydroxy_R1_file_name}_bismark_bt2_PE_report.txt"

hydroxy_R1_bismark_bam="${hydroxy_R1_file_path}_bismark_bt2_pe.bam"

methyl_R1="${temp_path}/$unmethyl_control/$hydroxymethyl_control/$methyl_control/${trimmed_R1_file_name}.fastq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1.fq.gz"
methyl_R2="${temp_path}/$unmethyl_control/$hydroxymethyl_control/$methyl_control/${trimmed_R2_file_name}.fastq.gz_unmapped_reads_2.fq.gz_unmapped_reads_2.fq.gz_unmapped_reads_2.fq.gz"

trimmed_genome_R1_bam="${temp_path}/$unmethyl_control/$hydroxymethyl_control/$methyl_control/genome_alignment/${trimmed_R1_file_name}.fastq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1_bismark_bt2_pe.deduplicated.bam.sorted.bam"

##################################################################################################################################
###########################################---STEP 4: PREVIOUSLY trim.sh---######################################################
##################################################################################################################################
if [ ! -f "${temp_path}/$R1_trimmed" ]; then
    cd $temp_path
        echo "$R1_trimmed does not exist yet"
        cutadapt -q 20 -m 15 -a AGATCGGAAGAGC -A AAATCAAAAAAAC -o $R1_trimmed -p $R2_trimmed $R1 $R2 --cores=${SLURM_CPUS_PER_TASK}
    #copy files back to seq_path directory
        rsync -vur --exclude "main_genome" --exclude "unmethyl_genome" --exclude "hydroxymethyl_genome" --exclude "methyl_genome" $temp_path/ $output_path

        echo "Fastq files have been trimmed"
else
        echo "Fastq files have already been trimmed"
fi

##################################################################################################################################
####################################---STEP 5: PREVIOUSLY map_to_control_seqs.sh---###############################################
##################################################################################################################################
#map to unmethyl control
if [ ! -f "${temp_path}/${unmethyl_control}/${trimmed_R1_file_name}_bismark_bt2_PE_report.txt" ]; then
    cd $temp_path
    	echo "${temp_path}/${unmethyl_control}/${trimmed_R1_file_name}_bismark_bt2_PE_report.txt does not exist yet" 
    	echo "Run bismark on $R1_trimmed and $R2_trimmed"
    bismark --bam --maxins 800 "$temp_path/unmethyl_genome" -1 $R1_trimmed -2 $R2_trimmed -o $temp_path/$unmethyl_control --unmapped --nucleotide_coverage --multicore $cores
        echo "map_to_control_seqs.sh complete for unmethyl control"
    rsync -vur --exclude "main_genome" --exclude "unmethyl_genome" --exclude "hydroxymethyl_genome" --exclude "methyl_genome" $temp_path/$unmethyl_control/ $output_path/$unmethyl_control
else
    	echo "Bismark mapping to control sequences for unmethyl control already completed"
fi

#map to hydroxymethyl control
unmapped_path=$temp_path/$unmethyl_control

if [ ! -f "$unmapped_path/$hydroxymethyl_control/${unmethyl_R1_file_name}_bismark_bt2_PE_report.txt" ]; then
    cd $unmapped_path
    	echo "$unmapped_path/$hydroxymethyl_control/${unmethyl_R1_file_name}_bismark_bt2_PE_report.txt does not exist yet" 
    	echo Run bismark on $R1_unmethyl and $R2_unmethyl
    bismark --bam --maxins 800 "$temp_path/hydroxymethyl_genome" -1 $R1_unmethyl -2 $R2_unmethyl -o $unmapped_path/$hydroxymethyl_control --unmapped --nucleotide_coverage --multicore $cores
        echo "map_to_control_seqs.sh complete for hydroxymethyl control"
    rsync -vur --exclude "main_genome" --exclude "unmethyl_genome" --exclude "hydroxymethyl_genome" --exclude "methyl_genome" $temp_path/$unmethyl_control/ $output_path/$unmethyl_control
else
    	echo "Bismark mapping to control sequences for hydroxymethyl control already completed"
fi


#map to methyl control
if [ ! -f "$unmapped_path/$hyrdoxymethyl_control/$methyl_control/${hydroxy_R1_file_name}_bismark_bt2_PE_report.txt" ]; then
    cd $unmapped_path/$hydroxymethyl_control
        echo "$unmapped_path/$hydroxymethyl_control/$methyl_control/${hydroxy_R1_file_name}_bismark_bt2_PE_report.txt does not exist yet" 
        echo Run bismark on $R1_hydroxy and $R2_hydroxy
    bismark --bam --maxins 800 "$temp_path/methyl_genome" -1 $R1_hydroxy -2 $R2_hydroxy -o $unmapped_path/$hydroxymethyl_control/$methyl_control --unmapped --nucleotide_coverage --multicore $cores
        echo "map_to_control_seqs.sh complete for methyl control"
    rsync -vur --exclude "main_genome" --exclude "unmethyl_genome" --exclude "hydroxymethyl_genome" --exclude "methyl_genome" $temp_path/$unmethyl_control/ $output_path/$unmethyl_control
else
        echo "Bismark mapping to control sequences for methyl control already completed"
fi

echo "map to control seqs complete"

##################################################################################################################################
###############################---STEP 6: PREVIOUSLY extract_methylation_controls.sh---###########################################
##################################################################################################################################
if [ ! -f "$temp_path/$unmethyl_control/${trimmed_R1_file_name}_bismark_bt2_pe_splitting_report.txt" ]; then  
    cd $temp_path/$unmethyl_control
        echo "$temp_path/$unmethyl_control/${trimmed_R1_file_name}_bismark_bt2_pe_splitting_report.txt does not exist yet"
    bismark_methylation_extractor --gzip --cytosine_report --bedGraph --genome_folder "$temp_path/unmethyl_genome" $trimmed_R1_bismark_bam --multicore $cores
        echo "extract_methylation_controls complete for unmethyl control"
    rsync -vur --exclude "main_genome" --exclude "unmethyl_genome" --exclude "hydroxymethyl_genome" --exclude "methyl_genome" $temp_path/$unmethyl_control/ $output_path/$unmethyl_control
else
    echo "methylation control extraction for the unmethyl control found and already created"
fi

if [ ! -f "$temp_path/$unmethyl_control/$hydroxymethyl_control/${unmethyl_R1_file_name}_bismark_bt2_pe_splitting_report.txt" ]; then 
    cd $temp_path/$unmethyl_control/$hydroxymethyl_control
        echo "$temp_path/$unmethyl_control/${unmethyl_R1_file_name}_bismark_bt2_pe_splitting_report.txt does not exist yet"
    bismark_methylation_extractor --gzip --cytosine_report --bedGraph --genome_folder "$temp_path/hydroxymethyl_genome" $unmethyl_R1_bismark_bam --multicore $cores
        echo "extract_methylation_controls complete for hydroxymethyl"
    rsync -vur --exclude "main_genome" --exclude "unmethyl_genome" --exclude "hydroxymethyl_genome" --exclude "methyl_genome" $temp_path/$unmethyl_control/ $output_path/$unmethyl_control
else
    echo "methylation control extraction for the hydroxymethyl control found and already created"
fi

if [ ! -f "$temp_path/$unmethyl_control/$hydroxymethyl_control/$methyl_control/${hydroxy_R1_file_name}_bismark_bt2_pe_splitting_report.txt" ]; then 
    cd $temp_path/$unmethyl_control/$hydroxymethyl_control/$methyl_control
        echo "$temp_path/$unmethyl_control/$hydroxymethyl_control/$methyl_control/${hydroxy_R1_file_name}_bismark_bt2_pe_splitting_report.txt does not exist yet"
    bismark_methylation_extractor --gzip --cytosine_report --bedGraph --genome_folder "$temp_path/methyl_genome" $hydroxy_R1_bismark_bam --multicore $cores
        echo "extract_methylation_controls complete for methyl"
    rsync -vur --exclude "main_genome" --exclude "unmethyl_genome" --exclude "hydroxymethyl_genome" --exclude "methyl_genome" $temp_path/$unmethyl_control/ $output_path/$unmethyl_control
else
    echo "methylation control extraction for the methyl control found and already created"
fi

echo "extract unmethyl, hydroxymethyl, and methyl controls complete"

##################################################################################################################################
####################################---STEP 7: PREVIOUSLY map_to_genome_seqs.sh---################################################
##################################################################################################################################
if [ ! -f "$temp_path/$unmethyl_control/$hydroxymethyl_control/$methyl_control/genome_alignment/${unmethyl_R1_file_name}.fq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1_bismark_bt2_pe.bam" ]; then   
    cd $temp_path/$unmethyl_control/$hydroxymethyl_control
        echo "$temp_path/$unmethyl_control/$hydroxymethyl_control/$methyl_control/genome_alignment/${unmethyl_R1_file_name}.fq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1_bismark_bt2_pe.bam does not exist yet"
    bismark --bam --maxins 1000 "$temp_path/main_genome" -1 $methyl_R1 -2 $methyl_R2 -o "$temp_path/$unmethyl_control/$hydroxymethyl_control/$methyl_control/genome_alignment" --unmapped --nucleotide_coverage --multicore $cores
        #--parallel ${SLURM_CPUS_PER_TASK}
    
    #copy all results back to #seq_path and updates to the genome
    rsync -vur --exclude "main_genome" --exclude "unmethyl_genome" --exclude "hydroxymethyl_genome" --exclude "methyl_genome" $temp_path/$unmethyl_control/ $output_path/$unmethyl_control  
    rsync -vur "$temp_path/main_genome/" $main_genome
    rsync -vur "$temp_path/hydroxymethyl_genome/" $hydroxymethyl_control_fasta
    rsync -vur "$temp_path/methyl_genome/" $methyl_control_fasta
    rsync -vur "$temp_path/unmethyl_genome/" $unmethyl_control_fasta
    echo "map_to_genome_seqs.sh complete"
else
   echo "genome alignment found and already completed"
fi

echo "map_to_genome_seqs complete"

##################################################################################################################################
#####################################---STEP 8: PREVIOUSLY remove_duplicates.sh---################################################
##################################################################################################################################
if [ ! -f "$temp_path/$unmethyl_control/$hydroxymethyl_control/$methyl_control/genome_alignment/${trimmed_R1_file_name}.fastq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1_bismark_bt2_pe.deduplication_report.txt" ]; then 
    cd "$temp_path/$unmethyl_control/$hydroxymethyl_control/$methyl_control/genome_alignment"
        echo "$temp_path/$unmethyl_control/$hydroxymethyl_control/$methyl_control/genome_alignment/${trimmed_R1_file_name}.fastq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1_bismark_bt2_pe.deduplication_report.txt does not exist"
            deduplicate_bismark -p --bam "${trimmed_R1_file_name}.fastq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1_bismark_bt2_pe.bam"
        echo "remove_duplicates complete"
    rsync -vur --exclude "main_genome" --exclude "unmethyl_genome" --exclude "hydroxymethyl_genome" --exclude "methyl_genome" $temp_path/$unmethyl_control/ $output_path/$unmethyl_control
else
    echo "duplicates already removed"
fi

if [ ! -f "$temp_path/$unmethyl_control/$hydroxymethyl_control/$methyl_control/genome_alignment/${trimmed_R1_file_name}.fastq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1_bismark_bt2_pe.deduplicated.bam.sorted.bam" ]; then
    cd "$temp_path/$unmethyl_control/$hydroxymethyl_control/$methyl_control/genome_alignment"
        echo "$temp_path/$unmethyl_control/$hydroxymethyl_control/$methyl_control/genome_alignment/${trimmed_R1_file_name}.fastq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1_bismark_bt2_pe.deduplicated.bam.sorted.bam does not yet exist"
    samtools sort "$temp_path/$unmethyl_control/$hydroxymethyl_control/$methyl_control/genome_alignment/${trimmed_R1_file_name}.fastq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1_bismark_bt2_pe.deduplicated.bam" -o  "$temp_path/$unmethyl_control/$hydroxymethyl_control/$methyl_control/genome_alignment/${trimmed_R1_file_name}.fastq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1_bismark_bt2_pe.deduplicated.bam.sorted.bam"
    samtools index "$temp_path/$unmethyl_control/$hydroxymethyl_control/$methyl_control/genome_alignment/${trimmed_R1_file_name}.fastq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1_bismark_bt2_pe.deduplicated.bam.sorted.bam"
    rsync -vur --exclude "main_genome" --exclude "unmethyl_genome" --exclude "hydroxymethyl_genome" --exclude "methyl_genome" $temp_path/$unmethyl_control/ $output_path/$unmethyl_control
    echo "bam file has been indexed and sorted"
else
    echo "indexed and sorted bam file already exists"
fi

##################################################################################################################################
####################################---STEP 9: PREVIOUSLY extract_methylation.sh---###############################################
##################################################################################################################################
if [ ! -f "$temp_path/$unmethyl_control/$hydroxymethyl_control/$methyl_control/genome_alignment/${trimmed_R1_file_name}.fastq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1_bismark_bt2_pe.bam.bismark.cov.gz" ]; then 
    cd "$temp_path/$unmethyl_control/$hydroxymethyl_control/$methyl_control/genome_alignment"
        echo "$temp_path/$unmethyl_control/$hydroxymethyl_control/$methyl_control/genome_alignment/${trimmed_R1_file_name}.fastq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1_bismark_bt2_pe.bam.bismark.cov.gz does not exist yet"
    #remove --cytosine_report
    bismark_methylation_extractor --buffer_size 20G --gzip --cytosine_report --bedGraph --genome_folder "$temp_path/main_genome" "$temp_path/$unmethyl_control/$hydroxymethyl_control/$methyl_control/genome_alignment/${trimmed_R1_file_name}.fastq.gz__unmapped_reads_1.fq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1_bismark_bt2_pe.deduplicated.bam.sorted.bam" --multicore $cores
        echo "extract_methylation complete"
    rsync -vur --exclude "main_genome" --exclude "unmethyl_genome" --exclude "methyl_genome" --exclude "hydroxymethyl_genome" $temp_path/$unmethyl_control/ $output_path/$unmethyl_control
else 
    echo "methlyation extraction already completed"
fi

##################################################################################################################################
##################################---STEP 10: PREVIOUSLY insert_size_analysis.sh---###############################################
##################################################################################################################################
if [ ! -f "$temp_path/$unmethyl_control/$hydroxymethyl_control/$methyl_control/genome_alignment/${trimmed_R1_file_name}.fastq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1_bismark_bt2_pe.bam_picard_insert_size_plot.pdf" ]; then 
   cd "$temp_path/$unmethyl_control/$hydroxymethyl_control/$methyl_control/genome_alignment"
        echo "$temp_path/$unmethyl_control/$hydroxymethyl_control/$methyl_control/genome_alignment/${trimmed_R1_file_name}.fastq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1_bismark_bt2_pe.bam_picard_insert_size_plot.pdf"
    module load R
    module load picard/2.9.5
    picard CollectInsertSizeMetrics INPUT=$trimmed_genome_R1_bam OUTPUT=$trimmed_genome_R1_bam\_picard_insert_size_metrics.txt HISTOGRAM_FILE=$trimmed_genome_R1_bam\_picard_insert_size_plot.pdf METRIC_ACCUMULATION_LEVEL=ALL_READS
        echo "picard insert_size_analysis complete"
    #copy files back to seq_path directory
    rsync -vur --exclude "main_genome" --exclude "unmethyl_genome" --exclude "hydroxymethyl_genome" --exclude "methyl_genome" $temp_path/$unmethyl_control/ $output_path/$unmethyl_control
else
    echo "picard insert size analysis already complete"
fi
