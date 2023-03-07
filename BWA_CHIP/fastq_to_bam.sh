#!/usr/bin/env bash

# Run scripts to enable activating conda environments
. $HOME/micromamba/etc/profile.d/conda.sh
. $HOME/micromamba/etc/profile.d/mamba.sh
mamba activate base

# Set bash options for verbose output and to fail immediately on errors or if variables are undefined.
set -o xtrace -o nounset -o pipefail -o errexit

# Parse command line arguments with getopt
arguments=$(getopt --options a --longoptions sample_name:,read_group:,reference_genome:,fastq_read1:,fastq_read2:,output_directory: --name 'fastq_to_bam' -- "$@")
eval set -- "$arguments"

while true; do
    case "$1" in 
        --sample_name )
            sample_name=$2 ; shift 2 ;;
        --read_group )
            read_group=$2 ; shift 2 ;;
        --reference_genome )
            reference_genome=$2 ; shift 2 ;;
        --fastq_read1 )
            fastq_read1=$2 ; shift 2 ;;
        --fastq_read2 )
            fastq_read2=$2 ; shift 2 ;;
        --output_directory )
            output_directory=$2 ; shift 2 ;;
        -- )
            shift; break ;;
        * )
            echo "Invalid argument ${1} ${2}" >&2
            exit 1
    esac
done

bam_name="${output_directory}/${sample_name}.bam"
mkdir -p "${output_directory}"

if [ ! -f "${bam_name}" ]; then
    echo "Aligning to reference genome with bwa-mem2..."
    set +o xtrace +o nounset; mamba activate --stack bwa; set -o xtrace -o nounset
    bwa-mem2 mem -R "${read_group}" "${reference_genome}" "${fastq_read1}" "${fastq_read2}" | samtools view -b - | samtools sort - -o "${bam_name}"
    set +o xtrace +o nounset; mamba deactivate; set -o xtrace -o nounset
    echo "...alignment complete"
else
    echo "Alignment already completed"
fi

if [ ! -f "${bam_name}.bai" ]; then
    echo "Indexing BAM file with samtools..."
    mamba run -n samtools samtools index "${bam_name}" "${bam_name}.bai"
    echo "...indexing complete"
else
    echo "Indexing of BAM already completed"
fi
