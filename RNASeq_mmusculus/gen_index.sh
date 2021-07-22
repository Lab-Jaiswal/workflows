#!/bin/bash

source "/labs/sjaiswal/workflows/RNASeq_mmusculus/workflow_config.sh"

module load star

STAR --runThreadN 16 `#Use 8 threads` \
    --runMode genomeGenerate `#Tell STAR you are generating an index` \
    --genomeDir ${star_index} `#Give location where index will be built`\
    --genomeFastaFiles ${reference_genome} `#Give location of reference genome FASTA`\
    --sjdbGTFfile ${annotation_gtf} `#Give location of gene annotation GTF`

