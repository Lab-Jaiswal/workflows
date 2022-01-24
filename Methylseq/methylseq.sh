#!/bin/bash

#SBATCH --time=32:00:00
#SBATCH --account=sjaiswal
#SBATCH --cpus-per-task=20
#SBATCH --mem=256GB
#SBATCH --job-name=methylseq

data_path=$1
seq_path="$1/fastq"
unmethyl_control_fasta=$2
chmod 775 $unmethyl_control_fasta
unmethyl_control=$3
hydroxymethyl_control_fasta=$4
chmod 775 $hydroxymethyl_control_fasta
hydroxymethyl_control=$5
methyl_control_fasta=$6
chmod 755 $methyl_control_fasta
methyl_control=$7
main_genome=$8
chmod 775 $main_genome
phix_path=$9
cores="${10}"

line_number=$SLURM_ARRAY_TASK_ID #get index of which file to process from $SLURM_ARRAY_TASK_ID provided by SLURM
fastq_file="${data_path}/fastq/FASTQs" #provide path to file containing list of fastq files
fastq_path="$(sed "${line_number}q; d" $fastq_file)" #extract only the line number corresponding to $SLURM_ARRAY_TASK_ID

sample_prefix=$(ls $fastq_path | grep fastq | head -n 1 | sed -e 's/_S[0-9]*_L[0-9]*_[IR][0-9]_[0-9]*.fastq.gz//g')
sample_name="${fastq_path##*/}"
filename=$(basename "${fastq_path}")
fastq_temp="${filename}"
PREFIX=$filename

seq_path="$data_path/fastq"

echo "copying FASTQs..."
temp_path=$(mktemp -d /tmp/tmp.XXXXXXXXXX)
echo "temp_path is: " $temp_path
mkdir $temp_path
#copy fastq files to temp_path
rsync -vur "$data_path/fastq/" $temp_path
echo "FASTQs have been copied to the temporary file"

echo "copying reference transcriptome..."
mkdir /tmp/genomes
mkdir -p /tmp/genomes/main_genome
mkdir -p /tmp/genomes/unmethyl_genome
mkdir -p /tmp/genomes/hydroxy_genome
mkdir -p /tmp/genomes/methyl_genome
rsync -vur "$main_genome/" "$temp_path/main_genome"
rsync -vur "$unmethyl_control_fasta/" "$temp_path/unmethyl_genome"
rsync -vur "$hydroxymethyl_control_fasta/" "$temp_path/hydroxymethyl_genome"
rsync -vur "$methyl_control_fasta/" "$temp_path/methyl_genome"

echo "Transcriptomes have been copied to the temporary file directory"

module load bismark/0.22.3
module load samtools
"Bismark and samtools modules have been loaded"

R1="${fastq_temp}_R1_001.fastq.gz"
R2="${fastq_temp}_R1_001.fastq.gz"

R1_trimmed="${fastq_temp}_R1_001.trimmed.fastq.gz"
R2_trimmed="${fastq_temp}_R2_001.trimmed.fastq.gz"
directory="$(dirname $R1_trimmed)"
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

#unmapped_R1_file_path=$(echo $R1_unmapped| sed 's/.fq.gz//')
unmethyl_R1_file_path="${temp_path}/$unmethyl_control/$hydroxymethyl_control/${trimmed_R1_file_name}.fastq.gz_unmapped_reads_1"
unmethyl_R1_file_name="${unmethyl_R1_file_path##*/}"
R1_unmethyl_bismark_PE_report="${unmethyl_R1_file_name}_bismark_bt2_PE_report.txt"

unmethyl_R1_bismark_bam="${unmethyl_R1_file_path}_bismark_bt2_pe.bam"

hydroxy_R1_file_path="${temp_path}/$unmethyl_control/$hydroxymethyl_control/$methyl_control/${trimmed_R1_file_name}.fastq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1"
hydroxy_R1_file_name="${hydroxy_R1_file_path##*/}"
R1_hydroxy_bismark_PE_report2="${hydroxy_R1_file_name}_bismark_bt2_PE_report.txt"

hydroxy_R1_bismark_bam="${hydroxy_R1_file_path}_bismark_bt2_pe.bam"

#unmapped_unmapped_R1="${temp_path}/$unmethyl_control/$hydroxymethyl_control/${trimmed_R1_file_name}.fastq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1.fq.gz"
#unmapped_unmapped_R2="${temp_path}/$unmethyl_control/$hydroxymethyl_control/${trimmed_R2_file_name}.fastq.gz_unmapped_reads_2.fq.gz_unmapped_reads_2.fq.gz"
methyl_R1="${temp_path}/$unmethyl_control/$hydroxymethyl_control/$methyl_control/${trimmed_R1_file_name}.fastq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1.fq.gz"
methyl_R2="${temp_path}/$unmethyl_control/$hydroxymethyl_control/$methyl_control/${trimmed_R2_file_name}.fastq.gz_unmapped_reads_2.fq.gz_unmapped_reads_2.fq.gz_unmapped_reads_2.fq.gz"

#trimmed_unmapped_R1_bam="${temp_path}/$unmethyl_control/$hydroxymethyl_control/genome_alignment/${trimmed_R1_file_name}.fastq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1_bismark_bt2_pe.bam"
trimmed_genome_R1_bam="${temp_path}/$unmethyl_control/$hydroxymethyl_control/$methyl_control/genome_alignment/${trimmed_R1_file_name}.fastq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1_bismark_bt2_pe.bam"

echo "temp_path is $temp_path"

#####################previously trim.sh###########################
###############################################################################
if [ ! -f "${temp_path}/$R1_trimmed" ]; then
    cd $temp_path
        echo "$R1_trimmed does not exist yet"
        cutadapt -q 20 -m 15 -a AGATCGGAAGAGC -A AAATCAAAAAAAC -o $R1_trimmed -p $R2_trimmed $R1 $R2 --cores=${SLURM_CPUS_PER_TASK}
    #copy files back to seq_path directory
        rsync -vur --exclude "main_genome" --exclude "unmethyl_genome" --exclude "hydroxymethyl_genome" --exclude "methyl_genome" $temp_path/ $seq_path

        echo "Fastq files have been trimmed"
else
        echo "Fastq files have already been trimmed"
fi

#####################previously map_to_control_seqs.sh###########################
###############################################################################

#map to unmethyl control
if [ ! -f "${temp_path}/${unmethyl_control}/${trimmed_R1_file_name}_bismark_bt2_PE_report.txt" ]; then
    cd $temp_path
    	echo "${temp_path}/${unmethyl_control}/${trimmed_R1_file_name}_bismark_bt2_PE_report.txt does not exist yet" 
    	echo "Run bismark on $R1_trimmed and $R2_trimmed"
    bismark --bam --maxins 800 "$temp_path/unmethyl_genome" -1 $R1_trimmed -2 $R2_trimmed -o $temp_path/$unmethyl_control --unmapped --nucleotide_coverage --multicore $cores
        echo "map_to_control_seqs.sh complete for unmethyl control"
    rsync -vur --exclude "main_genome" --exclude "unmethyl_genome" --exclude "hydroxymethyl_genome" --exclude "methyl_genome" $temp_path/ $seq_path
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
    rsync -vur --exclude "main_genome" --exclude "unmethyl_genome" --exclude "hydroxymethyl_genome" --exclude "methyl_genome" $temp_path/ $seq_path
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
    rsync -vur --exclude "main_genome" --exclude "unmethyl_genome" --exclude "hydroxymethyl_genome" --exclude "methyl_genome" $temp_path/ $seq_path
else
        echo "Bismark mapping to control sequences for methyl control already completed"
fi

echo "map to control seqs complete"


#####################previously extract_methylation_controls.sh###################
################################################################################

if [ ! -f "$temp_path/$unmethyl_control/${trimmed_R1_file_name}_bismark_bt2_pe_splitting_report.txt" ]; then  
    cd $temp_path/$unmethyl_control
        echo "$temp_path/$unmethyl_control/${trimmed_R1_file_name}_bismark_bt2_pe_splitting_report.txt does not exist yet"
    bismark_methylation_extractor --gzip --cytosine_report --bedGraph --genome_folder "$temp_path/unmethyl_genome" $trimmed_R1_bismark_bam --multicore $cores
        echo "extract_methylation_controls complete for unmethyl control"
    rsync -vur --exclude "main_genome" --exclude "unmethyl_genome" --exclude "hydroxymethyl_genome" --exclude "methyl_genome" $temp_path/ $seq_path
else
    echo "methylation control extraction for the unmethyl control found and already created"
fi

if [ ! -f "$temp_path/$unmethyl_control/$hydroxymethyl_control/${unmethyl_R1_file_name}_bismark_bt2_pe_splitting_report.txt" ]; then 
    cd $temp_path/$unmethyl_control/$hydroxymethyl_control
        echo "$temp_path/$unmethyl_control/${unmethyl_R1_file_name}_bismark_bt2_pe_splitting_report.txt does not exist yet"
    bismark_methylation_extractor --gzip --cytosine_report --bedGraph --genome_folder "$temp_path/hydroxymethyl_genome" $unmethyl_R1_bismark_bam --multicore $cores
        echo "extract_methylation_controls complete for hydroxymethyl"
    rsync -vur --exclude "main_genome" --exclude "unmethyl_genome" --exclude "hydroxymethyl_genome" --exclude "methyl_genome" $temp_path/ $seq_path
else
    echo "methylation control extraction for the hydroxymethyl control found and already created"
fi

if [ ! -f "$temp_path/$unmethyl_control/$hydroxymethyl_control/$methyl_control/${hydroxy_R1_file_name}_bismark_bt2_pe_splitting_report.txt" ]; then 
    cd $temp_path/$unmethyl_control/$hydroxymethyl_control/$methyl_control
        echo "$temp_path/$unmethyl_control/$hydroxymethyl_control/$methyl_control/${hydroxy_R1_file_name}_bismark_bt2_pe_splitting_report.txt does not exist yet"
    bismark_methylation_extractor --gzip --cytosine_report --bedGraph --genome_folder "$temp_path/methyl_genome" $hydroxy_R1_bismark_bam --multicore $cores
        echo "extract_methylation_controls complete for methyl"
    rsync -vur --exclude "main_genome" --exclude "unmethyl_genome" --exclude "hydroxymethyl_genome" --exclude "methyl_genome" $temp_path/ $seq_path
else
    echo "methylation control extraction for the methyl control found and already created"
fi

echo "extract unmethyl, hydroxymethyl, and methyl controls complete"

#####################previously report_controls.sh################################
################################################################################

if [ ! -f "$temp_path/$unmethyl_control/bismark_summary_report.txt" ]; then 
    cd $temp_path/$unmethyl_control
        echo "$temp_path/$unmethyl_control/bismark_summary_report.txt does not exist yet"
    #http://felixkrueger.github.io/Bismark/Docs/
    #bismark report options:
    #--alignment_report FILE
    #--dedup_report FILE
    #--splitting_report FILE
    #--mbias_report FILE
    #--nucleotide_report FILE
    
    #https://rawgit.com/FelixKrueger/Bismark/master/Docs/Bismark_User_Guide.html
    bismark2report
    bismark2summary
    #need to specify a nucleotide coverage report file in the above command!
    rsync -vur --exclude "main_genome" --exclude "unmethyl_genome" --exclude "hydroxymethyl_genome" --exclude "methyl_genome" $temp_path/ $seq_path
        echo "created report for unmethyl control"
else 
        echo "bismark summary already completed for unmethyl control"
fi
    
if [ ! -f "$temp_path/$unmethyl_control/$hydroxymethyl_control/bismark_summary_report.txt" ]; then
    cd $temp_path/$unmethyl_control/$hydroxymethyl_control
        echo "$temp_path/$unmethyl_control/$hydroxymethyl_control/bismark_summary_report.txt does not exist yet"
    bismark2report
    bismark2summary
    #need to specify a nucleotide coverage report file in the above command! ^
    rsync -vur --exclude "main_genome" --exclude "unmethyl_genome" --exclude "hydroxymethyl_genome" --exclude "methyl_genome" $temp_path/ $seq_path
        echo "created report for hydroxymethyl control"
else
    echo "bismark summary already completed for hydroxymethyl control"
fi

if [ ! -f "$temp_path/$unmethyl_control/$hydroxymethyl_control/$methyl_control/bismark_summary_report.txt" ]; then
    cd $temp_path/$unmethyl_control/$hydroxymethyl_control/$methyl_control
        echo "$temp_path/$unmethyl_control/$hydroxymethyl_control/$methyl_control/bismark_summary_report.txt does not exist yet"
    bismark2report
    bismark2summary
    #need to specify a nucleotide coverage report file in the above command! ^
    rsync -vur --exclude "main_genome" --exclude "unmethyl_genome" --exclude "hydroxymethyl_genome" --exclude "methyl_genome" $temp_path/ $seq_path
        echo "created report for methyl control"
else
    echo "bismark summary already completed for methyl control"
fi

echo "report_controls complete for unmethyl, hydroxymethyl, and methyl control sequences"

#####################previously map_to_genome_seqs.sh#############################
################################################################################
#if [ ! -f "$temp_path/$unmethyl_control/$hydroxymethyl_control/genome_alignment/${unmapped_R1_file_name}.fq.gz_unmapped_reads_1_bismark_bt2_pe.bam" ]; then   
#    cd $temp_path/$unmethyl_control/$hydroxymethyl_control
#        echo "$temp_path/$unmethyl_control/$hydroxymethyl_control/genome_alignment/${unmapped_R1_file_name}.fq.gz_unmapped_reads_1_bismark_bt2_pe.bam does not exist yet"
if [ ! -f "$temp_path/$unmethyl_control/$hydroxymethyl_control/$methyl_control/genome_alignment/${unmethyl_R1_file_name}.fq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1_bismark_bt2_pe.bam" ]; then   
    cd $temp_path/$unmethyl_control/$hydroxymethyl_control
        echo "$temp_path/$unmethyl_control/$hydroxymethyl_control/$methyl_control/genome_alignment/${unmethyl_R1_file_name}.fq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1_bismark_bt2_pe.bam does not exist yet"
    bismark --bam --maxins 1000 "$temp_path/main_genome" -1 $methyl_R1 -2 $methyl_R2 -o "$temp_path/$unmethyl_control/$hydroxymethyl_control/$methyl_control/genome_alignment" --unmapped --nucleotide_coverage --multicore $cores
        #--parallel ${SLURM_CPUS_PER_TASK}
    
    #copy all results back to #seq_path and updates to the genome
    rsync -vur --exclude "main_genome" --exclude "unmethyl_genome" --exclude "hydroxymethyl_genome" --exclude "methyl_genome" $temp_path/ $seq_path  
    rsync -vur "$temp_path/main_genome/" $main_genome
    rsync -vur "$temp_path/hydroxymethyl_genome/" $hydroxymethyl_control_fasta
    rsync -vur "$temp_path/methyl_genome/" $methyl_control_fasta
    rsync -vur "$temp_path/unmethyl_genome/" $unmethyl_control_fasta
    echo "map_to_genome_seqs.sh complete"
else
   echo "genome alignment found and already completed"
fi

echo "map_to_genome_seqs complete"

#####################previously remove_duplicates.sh##############################
################################################################################
if [ ! -f "$temp_path/$unmethyl_control/$hydroxymethyl_control/$methyl_control/genome_alignment/${trimmed_R1_file_name}.fastq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1_bismark_bt2_pe.deduplication_report.txt" ]; then 
    cd "$temp_path/$unmethyl_control/$hydroxymethyl_control/$methyl_control/genome_alignment"
        echo "$temp_path/$unmethyl_control/$hydroxymethyl_control/$methyl_control/genome_alignment/${trimmed_R1_file_name}.fastq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1_bismark_bt2_pe.deduplication_report.txt"
            deduplicate_bismark -p --bam $trimmed_genome_R1_bam  
        echo "remove_duplicates complete"
    rsync -vur --exclude "main_genome" --exclude "unmethyl_genome" --exclude "hydroxymethyl_genome" --exclude "methyl_genome" $temp_path/ $seq_path
else
    echo "duplicates already removed"
fi

if [ ! -f "$temp_path/$unmethyl_control/$hydroxymethyl_control/$methyl_control/genome_alignment/${trimmed_R1_file_name}.fastq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1_bismark_bt2_pe.deduplicated.bam.sorted.bam" ]; then
    cd "$temp_path/$unmethyl_control/$hydroxymethyl_control/$methyl_control/genome_alignment"
        echo "$temp_path/$unmethyl_control/$hydroxymethyl_control/$methyl_control/genome_alignment/${trimmed_R1_file_name}.fastq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1_bismark_bt2_pe.deduplicated.bam.sorted.bam does not yet exist"
    samtools sort "$temp_path/$unmethyl_control/$hydroxymethyl_control/$methyl_control/genome_alignment/${trimmed_R1_file_name}.fastq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1_bismark_bt2_pe.deduplicated.bam" -o  "$temp_path/$unmethyl_control/$hydroxymethyl_control/$methyl_control/genome_alignment/${trimmed_R1_file_name}.fastq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1_bismark_bt2_pe.deduplicated.bam.sorted.bam"
    samtools index "$temp_path/$unmethyl_control/$hydroxymethyl_control/$methyl_control/genome_alignment/${trimmed_R1_file_name}.fastq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1_bismark_bt2_pe.deduplicated.bam.sorted.bam"
    rsync -vur --exclude "main_genome" --exclude "unmethyl_genome" --exclude "hydroxymethyl_genome" --exclude "methyl_genome" $temp_path/ $seq_path
    echo "bam file has been indexed and sorted"
else
    echo "indexed and sorted bam file already exists"
fi


if [ ! -f "$temp_path/$unmethyl_control/$hydroxymethyl_control/$methyl_control/genome_alignment/split_bams/${trimmed_R1_file_name}.fastq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1_bismark_bt2_pe.deduplicated.bam.1.bam" ]; then
    cd "$temp_path/$unmethyl_control/$hydroxymethyl_control/$methyl_control/genome_alignment"
    echo "$temp_path/$unmethyl_control/$hydroxymethyl_control/$methyl_control/genome_alignment/split_bams/${trimmed_R1_file_name}.fastq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1_bismark_bt2_pe.deduplicated.bam.1.bam does not yet exist"
    #this depends if it is mouse or human genome!
    #chromosomes=$(seq 20)
    chromosomes=$(seq 22)
    chromosomes=(${chromosomes[@]} "X" "Y")
    echo $chromosomes
    mkdir split_bams
    chromosome="chr"
    for i in ${chromosomes[@]}
        do
            chr="${chromosome}${i}"
            echo "$chr"
            samtools view -b "${trimmed_R1_file_name}.fastq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1_bismark_bt2_pe.deduplicated.bam.sorted.bam" $chr > "split_bams/${trimmed_R1_file_name}.fastq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1_bismark_bt2_pe.deduplicated.bam.${chr}.${sample_name}.bam"  
        done
    rsync -vur --exclude "main_genome" --exclude "unmethyl_genome" --exclude "hydroxymethyl_genome" --exclude "methyl_genome" $temp_path/ $seq_path
    echo "bam file has been split up by chromosome"
else
    echo "bam file was already split up by chromosome"
fi

#####################previously insert_size_analysis.sh###########################
################################################################################
if [ ! -f "$temp_path/$unmethyl_control/$hydroxymethyl_control/$methyl_control/genome_alignment/${trimmed_R1_file_name}.fastq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1_bismark_bt2_pe.bam_picard_insert_size_plot.pdf" ]; then 
   cd "$temp_path/$unmethyl_control/$hydroxymethyl_control/$methyl_control/genome_alignment"
        echo "$temp_path/$unmethyl_control/$hydroxymethyl_control/$methyl_control/genome_alignment/${trimmed_R1_file_name}.fastq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1_bismark_bt2_pe.bam_picard_insert_size_plot.pdf"
    module load R
    module load picard/2.9.5
        #java -Xmx4g -XX:+UseG1GC -XX:ParallelGCThreads=8 -jar ~/picard.jar CollectInsertSizeMetrics INPUT=$sample OUTPUT={$sample}_picard_insert_size_metrics.txt HISTOGRAM_FILE={$sample}_picard_insert_size_plot.pdf METRIC_ACCUMULATION_LEVEL=ALL_READS
    picard CollectInsertSizeMetrics INPUT=$trimmed_genome_R1_bam OUTPUT=$trimmed_genome_R1_bam\_picard_insert_size_metrics.txt HISTOGRAM_FILE=$trimmed_genome_R1_bam\_picard_insert_size_plot.pdf METRIC_ACCUMULATION_LEVEL=ALL_READS
        echo "picard insert_size_analysis complete"
    #copy files back to seq_path directory
    rsync -vur --exclude "main_genome" --exclude "unmethyl_genome" --exclude "hydroxymethyl_genome" --exclude "methyl_genome" $temp_path/ $seq_path
else
    echo "picard insert size analysis already complete"
fi

#####################previously report.sh#########################################
################################################################################
if [ ! -f "$temp_path/$unmethyl_control/$hydroxymethyl_control/$methyl_control/genome_alignment/bismark_summary_report.txt" ]; then 
    #https://rawgit.com/FelixKrueger/Bismark/master/Docs/Bismark_User_Guide.html
    
    cd $temp_path/$unmethyl_control/$hydroxymethyl_control/$methyl_control/genome_alignment
    
    bismark2report
    bismark2summary
        echo "report complete"
    rsync -vur --exclude "main_genome" --exclude "unmethyl_genome" --exclude "hydroxymethyl_genome" --exclude "methyl_genome" $temp_path/ $seq_path
else
    echo "bismark summary found and already created"
fi
