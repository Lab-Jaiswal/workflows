#!/bin/bash

star_index="/labs/sjaiswal/genomes/GRCh38/STAR_genome" #path to STAR genome index

annotation_gtf="/labs/sjaiswal/genomes/GRCh38/gencode.v35.annotation.gtf" #path to annotation GTF

reference_genome="/labs/sjaiswal/genomes/GRCh38/GRCh38.p13.genome.fa" #location of reference genome used in alignment

annotation_refFlat="/labs/sjaiswal/genomes/GRCh38/gencode.v35.annotation.fixed.refFlat.txt" ## Created from ensembl GTF using gtfToGenePred -genePredExt -geneNameAsName2 from bioconda

gene_key="/labs/sjaiswal/genomes/mm10/gencode_vM20/gencode.vM20.primary_assembly.annotation.genekey.txt" #path to annotation GTF
