#!/bin/bash
#always mandatory
    reference_genome=~/References/GRCh38_full_analysis_set_plus_decoy_hla.fa
    gnomad_genomes=~/References/gnomad.genomes.v3.1.2.sites.maf05.vcf.bgz

#mandatory if using funcotator
    funcotator_sources=~/References/funcotator_dataSources.v1.6.20190124s

#mandatory if using split_by_chr
    sequence_dictionary=~/References/GRCh38.p12.genome.u2af1l5_mask.fa.dict

#Optional
    chr_intervals=~/References/whole_genome_intervals.interval_list
    intervals=~/References/CHIP_exons.interval_list
    transcript_list=~/References/chip_transcript_list.txt

    #
