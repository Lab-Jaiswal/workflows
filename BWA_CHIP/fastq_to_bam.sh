#!/bin/bash

echo "entering fastq_to_bam.sh"

SAMPLE_NAME=$1
READGROUP=$3
BWA_GREF=$4
R1=$5
R2=$6

echo "fastq_to_bam command used the following parameters:
$0 $1 $2 $3 $4 $5 $6"

if [ $SLURM_ARRAY_TASK_ID -eq 1 ]; then
  echo "arguments used for the fastq_to_bam.sh script:
        SAMPLE_NAME=$1
        READGROUP=$3
        BWA_GREF=$4
        R1=$5
        R2=$6
        " >> $parameter_file
fi

if [ ! -f "${SAMPLE_NAME}.bam" ]; then
    echo "Aligning to reference genome..."
    module load bwa
    module load samtools
    bwa mem -R "${READGROUP}" "${BWA_GREF}" "${R1}" "${R2}" | samtools view -b - | samtools sort - -o "${SAMPLE_NAME}.bam"
    echo "...alignment complete"
else
    echo "Alignment already completed"
fi

if [ ! -f "${SAMPLE_NAME}.bam.bai" ]; then
    echo "Indexing BAM file..."
    module load samtools
    samtools index "${SAMPLE_NAME}.bam" "${SAMPLE_NAM}.bam.bai"
    echo "...indexing complete"
else
    echo "Indexing of BAM already completed"
fi

echo "fastq_to_bam.sh complete"
