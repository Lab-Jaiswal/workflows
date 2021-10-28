#!/bin/bash

#star_index="/ifs/labs/sjaiswal/genomes/mm10/STAR_index"
#reference_genome="/ifs/labs/sjaiswal/genomes/mm10/GRCm38.primary_assembly.genome.fa"
#annotation_gtf="/ifs/labs/sjaiswal/genomes/mm10/gencode.vM20.annotation.gtf"

source "/labs/sjaiswal/workflows/RNASeq_hsapiens/workflow_config.sh"

module load star

STAR --runThreadN 8 `#Use 8 threads` \
    --runMode genomeGenerate `#Tell STAR you are generating an index` \
    --genomeDir ${star_index} `#Give location where index will be built`\
    --genomeFastaFiles ${reference_genome} `#Give location of reference genome FASTA`\
    --sjdbGTFfile ${annotation_gtf} `#Give location of gene annotation GTF`

