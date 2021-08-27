#!/bin/bash

#SBATCH --time=3:00:00
#SBATCH --account=sjaiswal
#SBATCH --cpus-per-task=8
#SBATCH --mem=32GB
#SBATCH --job-name=CHIP_variant_call

parent_directory=$1
output_directory="$2"
min_coverage="$3"
min_var_freq="$4"
p_value="$5"

line_number=$SLURM_ARRAY_TASK_ID #get index of which file to process from $SLURM_ARRAY_TASK_ID provided by SLURM
fastq_file="${parent_directory}/fastq_files" #provide path to file containing list of fastq files
fastq_prefix="$(sed "${line_number}q; d" "${fastq_file}")" #extract only the line number corresponding to $SLURM_ARRAY_TASK_ID

FILENAME=$(basename "${fastq_prefix}")
PREFIX=$FILENAME

R1="${fastq_prefix}_R1_001.fastq.gz"
R2="${fastq_prefix}_R2_001.fastq.gz"
READGROUP="@RG\tID:${PREFIX}\tLB:${PREFIX}\tPL:illumina\tSM:${PREFIX}"

BWA_GREF="/oak/stanford/groups/sjaiswal/Herra/CHIP_Panel_AmpliSeq/GRCh38.p12.genome.u2af1l5_mask.fa" #reference genome
TWIST_SNPS="/labs/sjaiswal/workflows/BWA_mutect_twist/twist_snps.bed" #SNPs for germline calling
ASSEMBLY="GRCh38" #Genome version

#Run script and save files in the same location of the fastq files
cd "${output_directory}" || exit

##align reads and index
#if [ ! -f ${PREFIX}_${ASSEMBLY}_R1_noadapt.fastq ] | [ ! -f ${PREFIX}_${ASSEMBLY}_R2_noadapt.fastq ]; then
    #echo "Trimming adapters with cutadapt"
    #module load cutadapt/2.4
    #cutadapt \
        #-b AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        #-b AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        #-b TGACTGGAGTTCAGACGTGTGCTCTTCCGATCT \
        #-b ACACTCTTTCCCTACACGACGCTCTTCCGATCT \
        #-B AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        #-B AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        #-B TGACTGGAGTTCAGACGTGTGCTCTTCCGATCT \
        #-B ACACTCTTTCCCTACACGACGCTCTTCCGATCT \
        #--too-short-paired-output ${PREFIX}_${ASSEMBLY}_short_R2_output.fastq \
        #--too-short-output ${PREFIX}_${ASSEMBLY}_short_R1_output.fastq \
        #--too-long-paired-output ${PREFIX}_${ASSEMBLY}_long_R2_output.fastq \
        #--too-long-output ${PREFIX}_${ASSEMBLY}_long_R1_output.fastq \
        #--report=minimal \
        #-n 3 \
        #-m 100 \
        #-o ${PREFIX}_${ASSEMBLY}_R1_noadapt.fastq \
        #-p ${PREFIX}_${ASSEMBLY}_R2_noadapt.fastq \
        #$R1 \
        #$R2 
    #echo "...adapters trimmed"
#else
    #echo "Adapters already trimmed"
#fi

if [ ! -f "${PREFIX}_${ASSEMBLY}.bam" ]; then
    echo "Aligning to reference genome..."
    module load bwa
    module load samtools
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

#if [ ! -f "${PREFIX}_${ASSEMBLY}_depth" ]; then
    #echo "Computing base pair depth..."
    #module load samtools
    #samtools depth "${PREFIX}_${ASSEMBLY}.bam" > "${PREFIX}_${ASSEMBLY}_depth"
    #echo "...base pair depth generated"
#else
    #echo "Base pair depth already generated"
#fi

#if [ ! -f "${PREFIX}_${ASSEMBLY}_depth_amplicons.txt" ]; then
    #echo "Computing average amplicon depth..."
    #module load R
    #SCRIPT_PATH=$(scontrol show job "${SLURM_JOBID}" | awk -F= '/Command=/{print $2}' | cut -d ' ' -f 1)
    #Rscript "$(dirname "${SCRIPT_PATH}")/get_amplicon_depth.R" "${PREFIX}_${ASSEMBLY}_depth" "$(dirname "${SCRIPT_PATH}")/chip_submitted_targets_Twist.xls"
    #echo "...average amplicon depth generated"
#else
    #echo "Average amplicon depth already generated"
#fi

if [ ! -f "${PREFIX}_${ASSEMBLY}_mutect2.vcf" ]; then
    echo "Calling somatic variants with Mutect2..."
    module load gatk4
    gatk Mutect2 \
        --input "${PREFIX}_${ASSEMBLY}.bam" \
        -tumor "${PREFIX}" \
        --output "${PREFIX}_${ASSEMBLY}_mutect2.vcf" \
        --reference "${BWA_GREF}" \
        --dont-use-soft-clipped-bases \
        --bamout "${PREFIX}_${ASSEMBLY}_mutect2.bam"

    module load samtools
    samtools index "${PREFIX}_${ASSEMBLY}_mutect2.bam" "${PREFIX}_${ASSEMBLY}_mutect2.bam.bai"

    echo "...somatic variants called."
else
    echo "Mutect2 somatic variants already called"
fi

if [ ! -f "${PREFIX}_${ASSEMBLY}_mutect2_filter.vcf" ]; then
    echo "Filtering somatic variants with FilterMutectCalls..."
    module load gatk4
    gatk FilterMutectCalls \
        --variant "${PREFIX}_${ASSEMBLY}_mutect2.vcf" \
        --output "${PREFIX}_${ASSEMBLY}_mutect2_filter.vcf" \
        --reference "${BWA_GREF}"
    echo "...somatic variants filtered."
else
    echo "Mutect2 somatic variants already filtered"
fi

if [ ! -f "${PREFIX}_${ASSEMBLY}_mutect2_filter_funcotator.vcf" ]; then
    TRANSCRIPT_LIST="/oak/stanford/groups/sjaiswal/Herra/CHIP_TWIST-PANEL_ATHEROMA/chip_transcript_list.txt" #Transcript list for Mutect
    FUNCOTATOR_SOURCES="/labs/sjaiswal/tools/funcotator/funcotator_dataSources.v1.6.20190124s" #Reference for Funcotator

    echo "Annotating Mutect2 VCF with Funcotator..."
    module load gatk4
    gatk Funcotator \
         --variant "${PREFIX}_${ASSEMBLY}_mutect2_filter.vcf" \
         --reference "${BWA_GREF}" \
         --ref-version hg38 \
         --data-sources-path "${FUNCOTATOR_SOURCES}" \
         --transcript-list "${TRANSCRIPT_LIST}" \
         --output "${PREFIX}_${ASSEMBLY}_mutect2_filter_funcotator.vcf" \
         --output-file-format VCF 
    echo "...VCF annotated."
else
    echo "Mutect2 VCF already annotated"
fi

if [ ! -f "${PREFIX}_${ASSEMBLY}_mutect2_filter_funcotator.maf" ]; then
    echo "Annotating VCF with Funcotator (MAF output)..."
    module load gatk4
    gatk Funcotator \
         --variant "${PREFIX}_${ASSEMBLY}_mutect2_filter.vcf" \
         --reference "${BWA_GREF}" \
         --ref-version hg38 \
         --data-sources-path "${FUNCOTATOR_SOURCES}" \
         --transcript-list "${TRANSCRIPT_LIST}" \
         --output "${PREFIX}_${ASSEMBLY}_mutect2_filter_funcotator.maf" \
         --output-file-format MAF 
         echo "...VCF annotated (MAF ouput)."
else
    echo "Mutect2 VCF already annotated (MAF output)"
fi

if [ ! -f "${PREFIX}_${ASSEMBLY}_haplotypecaller.gvcf" ]; then
    echo "Calling germline variants with HaplotypeCaller..."
    module load gatk4
    gatk HaplotypeCaller \
        --input "${PREFIX}_${ASSEMBLY}.bam" \
        --output "${PREFIX}_${ASSEMBLY}_haplotypecaller.gvcf" \
        --reference "${BWA_GREF}" \
        --intervals "${TWIST_SNPS}" \
        --bamout "${PREFIX}_${ASSEMBLY}_haplotypecaller.bam" \
        --emit-ref-confidence GVCF

    module load samtools
    samtools index "${PREFIX}_${ASSEMBLY}_haplotypecaller.bam" "${PREFIX}_${ASSEMBLY}_haplotypecaller.bam.bai"

    echo "...germline variants called."
else
    echo "HaplotypeCaller gVCF already exists"
fi

if [ ! -f "${PREFIX}_${ASSEMBLY}_haplotypecaller_genotypes.vcf" ]; then
    echo "Genotyping germline variants in gVCF with GenotypeGVCF..."
    module load gatk4
    gatk GenotypeGVCFs \
        --variant "${PREFIX}_${ASSEMBLY}_haplotypecaller.gvcf" \
        --output "${PREFIX}_${ASSEMBLY}_haplotypecaller_genotypes.vcf" \
        --reference "${BWA_GREF}" \
        --intervals "${TWIST_SNPS}" \
        --include-non-variant-sites

    echo "...germline variants genotyped."
else
    echo "Germline variants already genotyped with GenotypeGVCFs"
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
