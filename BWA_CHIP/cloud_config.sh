#!/bin/bash
#always mandatory
    reference_genome="project-G5B07V8JPkg740v9GjfF9PzV:/References/GRCh38.p12.genome.u2af1l5_mask.fa"
    gnomad_genomes="project-G5B07V8JPkg740v9GjfF9PzV:/References/gnomad.genomes.v3.1.2.sites.maf05.vcf.bgz"

#mandatory if using funcotator
    funcotator_sources="project-G5B07V8JPkg740v9GjfF9PzV:/References/funcotator_dataSources.v1.6.20190124s"

#mandatory if using split_by_chr
    sequence_dictionary="project-G5B07V8JPkg740v9GjfF9PzV:/References/GRCh38.p12.genome.u2af1l5_mask.fa.dict"

#Optional
    chr_intervals="project-G5B07V8JPkg740v9GjfF9PzV:/References/whole_genome_intervals.interval_list"
    intervals="project-G5B07V8JPkg740v9GjfF9PzV:/References/CHIP_exons.interval_list"
    transcript_list="project-G5B07V8JPkg740v9GjfF9PzV:/References/chip_transcript_list.txt"

    #
