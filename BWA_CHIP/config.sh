#!/bin/bash
#always mandatory
    reference_genome="/oak/stanford/groups/smontgom/maurertm/ADRC/Cloud_Testing_Folder/Params/GRCh38.p12.genome.u2af1l5_mask.fa"
    gnomad_genomes="/oak/stanford/groups/smontgom/maurertm/ADRC/Cloud_Testing_Folder/Params/gnomad.genomes.v3.1.2.sites.maf05.vcf.bgz"

#mandatory if using funcotator
    funcotator_sources="/oak/stanford/groups/smontgom/maurertm/ADRC/Cloud_Testing_Folder/Params/funcotator_dataSources.v1.6.20190124s"

#mandatory if using split_by_chr
    sequence_dictionary="/oak/stanford/groups/smontgom/maurertm/ADRC/Cloud_Testing_Folder/Params/GRCh38.p12.genome.u2af1l5_mask.fa.dict"

#Optional
    chr_intervals="${code_directory}/whole_genome_intervals.interval_list"
    intervals="/oak/stanford/groups/smontgom/maurertm/ADRC/Cloud_Testing_Folder/Params/CHIP_exons.interval_list"
    transcript_list="/oak/stanford/groups/smontgom/maurertm/ADRC/Cloud_Testing_Folder/Params/chip_transcript_list.txt"

    #
