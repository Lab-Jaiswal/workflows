#!/bin/bash

#SBATCH --time=1:00:00
#SBATCH --account=sjaiswal
#SBATCH --cpus-per-task=8
#SBATCH --mem=64GB
#SBATCH --job-name=methylseq
        
data_path=$1
seq_path="$1/fastq"
unmethyl_control_fasta=$2
unmethyl_control=$3
hydroxymethyl_control_fasta=$4
hydroxymethyl_control=$5
genome_path=$6
phix_path=$7
temp_path=$8
cores=$9

#####################previous map_to_control_seqs.sh###########################
###############################################################################
unmethyl_control_count=$(find "$data_path/fastq" -type d | grep "$unmethyl_control" | wc -l)

if [ $unmethyl_control_count -lt 1 ]; then
    
    module load bismark/0.22.3
    
    temp_genomes_path=/tmp/genomes
    
    #WHAT IS RIGHT HERE?
    genomes_path=$(dirname $genome_path)
    rsync -vur "$genomes_path/$unmethyl_control_fasta" $temp_genomes_path
    #rsync -vur /oak/stanford/groups/smontgom/kameronr/methylseq/JG97/$unmethyl_control_fasta $temp_genomes_path
    rsync -vur "$genomes_path/$hydroxymethyl_control_fasta" $temp_genomes_path
    
    #rsync -vur /oak/stanford/groups/smontgom/kameronr/methylseq/JG97/$hydroxymethyl_control_fasta $temp_genomes_path
    
    cd $temp_path
    
    for read1 in *R1_001.trimmed.fastq.gz
    do 
        read2=$(echo $read1| sed 's/R1/R2/')   
        file_name=$(echo $read1| sed 's/.fastq.gz//')
        file_name_check=$file_name\_bismark_bt2_PE_report.txt
        if [ ! -f $temp_path/$unmethyl_control/$file_name_check ]; then
    	echo $file_name_check does not exist yet
    	echo Run bismark on $read1 and $read2
    	bismark --bam --maxins 800 $temp_genomes_path/$unmethyl_control_fasta -1 $read1 -2 $read2 -o $temp_path/$unmethyl_control --unmapped --nucleotide_coverage --multicore $cores
        else
    	echo Skipping Bismark mapping for $read1
        fi
        
    done
    
    
    unmapped_path=$temp_path/$unmethyl_control
    cd $unmapped_path
    
    for read1 in *_unmapped_reads_1.fq.gz
    do
        read2=$(echo $read1| sed 's/R1/R2/' | sed 's/_unmapped_reads_1/_unmapped_reads_2/')
        file_name=$(echo $read1| sed 's/.fastq.gz//')
        file_name_check=$file_name\_bismark_bt2_PE_report.txt
        if [ ! -f $temp_path/$unmapped_path/$hydroxymethyl_control/$file_name_check ]; then
    	echo $file_name_check does not exist yet
    	echo Run bismark on $read1 and $read2
    	bismark --bam --maxins 800 $temp_genomes_path/$hydroxymethyl_control_fasta -1 $read1 -2 $read2 -o $unmapped_path/$hydroxymethyl_control --unmapped --nucleotide_coverage --multicore $cores
    #--parallel ${SLURM_CPUS_PER_TASK}
        else
    	echo Skipping Bismark mapping for $read1
        fi
    done
    
    rsync -vur $temp_path/ $seq_path
    rsync -vur $temp_genomes_path/ "$data_path/genomes/"

else
    echo "unmethyl directory found and already created"
fi

echo "map_to_control_seqs.sh complete"
#####################previous extract_methylation_controls.sh###################
################################################################################
hydroxymethyl_control_count=$(find "$data_path/fastq/$unmethyl_control" -type d | grep "$hydroxymethyl_control" | wc -l)

if [ $hydroxymethyl_control_count -lt 1 ]; then  
    #temp_genomes_path=/tmp/genomes
    #WHAT IS RIGHT HERE?
    #rsync -vur "$genome_path/$unmethyl_control_fasta" $temp_genomes_path
    #rsync -vur "$genome_path/$hydroxymethyl_control_fasta" $temp_genomes_path
    
    cd $temp_path/$unmethyl_control
    
    for sample in *.trimmed_bismark_bt2_pe.bam
    do
            bismark_methylation_extractor --gzip --cytosine_report --bedGraph --genome_folder $temp_genomes_path/$unmethyl_control_fasta $sample --multicore $cores
    done
    
    cd $temp_path/$unmethyl_control/$hydroxymethyl_control
    
    for sample in *.trimmed.fastq.gz_unmapped_reads_1_bismark_bt2_pe.bam
    do
            bismark_methylation_extractor --gzip --cytosine_report --bedGraph --genome_folder $temp_genomes_path/$hydroxymethyl_control_fasta $sample --multicore $cores
    done
    
    rsync -vur $temp_path/ $seq_path
else
    echo "hydroxymethyl directory found and already created"
fi

echo "extract_methylation_controls complete"
#####################previous report_controls.sh################################
################################################################################
report_controls_count=$(find "$data_path/fastq/$unmethyl_control/$hydroxymethyl_control" -type f | grep "bismark_summary_report" | wc -l)

if [ $report_controls_count -le 1 ]; then 
    cd $temp_path/$unmethyl_control
    
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
    
    cd $temp_path/$unmethyl_control/$hydroxymethyl_control
    
    bismark2report
    bismark2summary
    #need to specify a nucleotide coverage report file in the above command! ^
    
    rsync -vur $temp_path/ $seq_path
else
    echo "bismark summary found and already created"
fi

echo "report_controls complete"
#####################previous map_to_genome_seqs.sh#############################
################################################################################
genome_alignment_count=$(find "$data_path/fastq/$unmethyl_control/$hydroxymethyl_control" -type d | grep "genome_alignment" | wc -l)
if [ $genome_alignment_count -lt 1 ]; then
    temp_genomes_path=/tmp/genomes
    genome_name=$(basename $genome_path)
    
    rsync -vur $genome_path $temp_genomes_path
    
    cd $temp_path/$unmethyl_control/$hydroxymethyl_control
    
    for read1 in *_unmapped_reads_1.fq.gz_unmapped_reads_1.fq.gz
        do
            read2=$(echo $read1| sed 's/R1/R2/' | sed 's/_unmapped_reads_1.fq.gz_unmapped_reads_1/_unmapped_reads_2.fq.gz_unmapped_reads_2/')
        #    read2=$(echo $read1| sed 's/R1/R2/') 
            bismark --bam --maxins 1000 $temp_genomes_path/$genome_name -1 $read1 -2 $read2 -o $temp_path/$unmethyl_control/$hydroxymethyl_control/genome_alignment --unmapped --nucleotide_coverage --multicore $cores
        #--parallel ${SLURM_CPUS_PER_TASK}
    done
    
    #copy all results back to #seq_path and updates to the genome
    rsync -vur $temp_path/$unmethyl_control/$hydroxymethyl_control/ $seq_path/$unmethyl_control/$hydroxymethyl_control/
    rsync -vur $temp_genomes_path/$genome_name $genome_path/..
    
else
    echo "genome alignment found and already completed"
fi

echo "map_to_genome_seqs complete"
#####################previous remove_duplicates.sh##############################
################################################################################
genome_clean_count=$(find "$data_path/fastq/$unmethyl_control/$hydroxymethyl_control/genome_alignment" -type f | grep "gz_unmapped_reads_1_bismark_bt2_pe.bam" | wc -l)
if [ $genome_clean_count -le 1 ]; then
    cd $temp_path/$unmethyl_control/$hydroxymethyl_control/genome_alignment
        
    for sample in *gz_unmapped_reads_1_bismark_bt2_pe.bam
    do
        #java -Xmx4g -XX:+UseG1GC -XX:ParallelGCThreads=8 -jar ~/picard.jar CollectInsertSizeMetrics INPUT=$sample OUTPUT={$sample}_picard_insert_size_metrics.txt HISTOGRAM_FILE={$sample}_picard_insert_size_plot.pdf METRIC_ACCUMULATION_LEVEL=ALL_READS
        #picard CollectInsertSizeMetrics INPUT=$sample OUTPUT={$sample}_picard_insert_size_metrics.txt HISTOGRAM_FILE={$sample}_picard_insert_size_plot.pdf METRIC_ACCUMULATION_LEVEL=ALL_READS
        deduplicate_bismark -p --bam $sample
    done
        
    #copy files back to seq_path directory
    rsync -vur $temp_path/ $seq_path
    
echo "remove_duplicates complete"

#####################previous extract_methylation.sh############################
################################################################################
    temp_genomes_path=/tmp/genomes
    rsync -vur $genome_path $temp_genomes_path
    genome_name=$(basename $genome_path)
        
    cd $temp_path/$unmethyl_control/$hydroxymethyl_control/genome_alignment/ 
    for sample in *.bam
    do
        bismark_methylation_extractor --buffer_size 20G --gzip --cytosine_report --bedGraph --genome_folder $temp_genomes_path/$genome_name $sample --multicore $cores
    done
        
    #copy files back to seq_path directory
    rsync -vur $temp_path/ $seq_path
    
    echo "extract_methylation complete"
#####################previous insert_size_analysis.sh###########################
################################################################################
    module load R
    module load picard/2.9.5
        
    cd $temp_path/$unmethyl_control/$hydroxymethyl_control/genome_alignment
        
        
    for sample in *gz_unmapped_reads_1_bismark_bt2_pe.bam
    do
        #java -Xmx4g -XX:+UseG1GC -XX:ParallelGCThreads=8 -jar ~/picard.jar CollectInsertSizeMetrics INPUT=$sample OUTPUT={$sample}_picard_insert_size_metrics.txt HISTOGRAM_FILE={$sample}_picard_insert_size_plot.pdf METRIC_ACCUMULATION_LEVEL=ALL_READS
        picard CollectInsertSizeMetrics INPUT=$sample OUTPUT=$sample\_picard_insert_size_metrics.txt HISTOGRAM_FILE=$sample\_picard_insert_size_plot.pdf METRIC_ACCUMULATION_LEVEL=ALL_READS
        
    done
        
    #copy files back to seq_path directory
    rsync -vur $temp_path/ $seq_path
else
    echo "genome alignment cleaning already complete"
fi

echo "insert_size_analysis complete"

#####################previous report.sh#########################################
################################################################################
report_count=$(find "$data_path/fastq/$unmethyl_control/$hydroxymethyl_control/genome_alignment" -type f | grep "bismark_summary_report" | wc -l)

if [ $report_count -le 1 ]; then
    #https://rawgit.com/FelixKrueger/Bismark/master/Docs/Bismark_User_Guide.html
    
    cd $temp_path/$unmethyl_control/$hydroxymethyl_control/genome_alignment
    
    bismark2report
    bismark2summary
    #need to specify a nucleotide coverage report file in the above command! ^
    
    #cd /oak/stanford/groups/sjaiswal/kameronr/sjaiswal_old/bisulfiteseq/bisulfite_run180925/
    
    #for sample in *.trimmed_bismark_bt2_pe.bam
    #do
    #    bismark_methylation_extractor --gzip --cytosine_report --bedGraph --genome_folder /home/kameronr/sjaiswal/kameronr/JG97/methylseq/miseq/fasta_unmethylated $sample
    #done
    
    #copy files back to seq_path directory
    rsync -vur $temp_path/ $seq_path
else
    echo "bismark summary found and already created"
fi

echo "report complete"
