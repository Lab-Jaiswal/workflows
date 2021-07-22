#!/bin/bash

#SBATCH --time=3:00:00
#SBATCH --account=sjaiswal
#SBATCH --cpus-per-task=8
#SBATCH --mem=32GB
#SBATCH --job-name=variant_call

parent_directory=$1
output_directory="$2"
min_coverage="$3"
min_var_freq="$4"
p_value="$5"

line_number=$SLURM_ARRAY_TASK_ID #get index of which file to process from $SGE_TASK_ID provided by gridengine
fastq_file="${parent_directory}/fastq_files" #provide path to file containing list of fastq files
fastq_prefix="$(sed "${line_number}q; d" "${fastq_file}")" #extract only the line number corresponding to $SGE_TASK_ID

FILENAME=$(basename "$fastq_prefix")
PREFIX=$FILENAME

R1="${fastq_prefix}_R1_001.fastq.gz"
R2="${fastq_prefix}_R2_001.fastq.gz"
READGROUP="@RG\tID:${PREFIX}\tLB:${PREFIX}\tPL:illumina\tSM:${PREFIX}"

BWA_GREF="/oak/stanford/groups/sjaiswal/Herra/CHIP_Panel_AmpliSeq/GRCh38.p12.genome.u2af1l5_mask.fa"
ASSEMBLY="GRCh38"

#Run Script and save files in the same location of the fastq files
cd "${output_directory}" || exit

#align reads and index
if [ ! -f "${PREFIX}_${ASSEMBLY}.bam" ]; then
    echo "Aligning to reference genome..."
    module load samtools
    module load bwa
    bwa mem -R "${READGROUP}" "${BWA_GREF}" "${R1}" "${R2}" | samtools view -b - | samtools sort - -o "${PREFIX}_${ASSEMBLY}.bam"
    echo "...alignment complete"
else
    echo "Alignment already completed"
fi

if [ ! -f "${PREFIX}_${ASSEMBLY}.bam.bai" ]; then
    echo "Indexing BAM file..."
    module load samtools
    samtools index "${PREFIX}_${ASSEMBLY}.bam" "${PREFIX}_${ASSEMBLY}.bam.bai"
    echo "...indexing complete"
else
    echo "Indexing of BAM already completed"
fi

if [ ! -f "${PREFIX}_${ASSEMBLY}.pileup" ]; then
    module load samtools
    echo "Generating pileup from BAM..."
    samtools mpileup \
        -A \
        --max-depth 0 \
        -C50 \
        -f "${BWA_GREF}" \
        "${PREFIX}_${ASSEMBLY}.bam" > "${PREFIX}_${ASSEMBLY}.pileup"
    echo "...pileup generated"
else
    echo "Pileup already generated"
fi

if [ ! -f "${PREFIX}_${ASSEMBLY}_varscan2.vcf" ]; then
    echo "Calling variants from pileup..."
    module load varscan
    varscan mpileup2cns \
        "${PREFIX}_${ASSEMBLY}.pileup" \
        --min-coverage "${min_coverage}" \
        --min-var-freq "${min_var_freq}" \
        --p-value "${p_value}" \
        --output-vcf 1 > "${PREFIX}_${ASSEMBLY}_varscan2.vcf"
    echo "...variants called"
else
    echo "Variants already called"
fi

if [ ! -f "${PREFIX}_${ASSEMBLY}_varscan2_filter.vcf" ]; then
    echo "Filtering variants in VCF..."
    varscan filter \
        "${PREFIX}_${ASSEMBLY}_varscan2.vcf" \
        --output-file "${PREFIX}_${ASSEMBLY}_varscan2_filter.vcf" \
        --min-coverage "${min_coverage}" \
        --min-var-freq "${min_var_freq}" \
        --p-value "${p_value}"
    echo "...variants filtered"
else
    echo "Variants already filtered"
fi


if [ ! -f "${PREFIX}_${ASSEMBLY}_varscan2_filter_annovar.hg38_multianno.vcf" ]; then
    echo "Annotating VCF with Annovar..."
    ASSEMBLY_REFGENE="hg38"
    ANNOVARROOT="/labs/sjaiswal/tools/annovar"
    "${ANNOVARROOT}/table_annovar.pl" \
        "${PREFIX}_${ASSEMBLY}_varscan2_filter.vcf" \
        "${ANNOVARROOT}/humandb" \
        --buildver "${ASSEMBLY_REFGENE}" \
        --remove \
        --outfile "${PREFIX}_${ASSEMBLY}_varscan2_filter_annovar" \
        --protocol ensGene \
        --operation g \
        --nastring '.' \
        --vcfinput \
        --thread 1
    echo "...VCF annotated"
else
    echo "VCF already annotated"
fi
