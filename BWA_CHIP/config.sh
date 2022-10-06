#!/bin/bash
#always mandatory
    reference_genome="/oak/stanford/groups/sjaiswal/Herra/CHIP_Panel_AmpliSeq/GRCh38.p12.genome.u2af1l5_mask.fa"
    gnomad_genomes="/oak/stanford/groups/smontgom/dnachun/workflows/gnomad.genomes.v3.1.2.sites.maf05.vcf.bgz"

#mandatory if using funcotator
    funcotator_sources="/labs/sjaiswal/tools/funcotator/funcotator_dataSources.v1.6.20190124s"

#mandatory if using split_by_chr
    sequence_dictionary="/oak/stanford/groups/sjaiswal/Herra/CHIP_Panel_AmpliSeq/GRCh38.p12.genome.u2af1l5_mask.fa.dict"

#Optional
    chr_intervals="${code_directory}/whole_genome_intervals.interval_list"
    intervals="${code_directory}/CHIP_exons.interval_list"
    transcript_list="/oak/stanford/groups/sjaiswal/Herra/CHIP_TWIST-PANEL_ATHEROMA/chip_transcript_list.txt"
