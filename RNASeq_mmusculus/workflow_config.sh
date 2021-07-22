#!/bin/bash

star_index="/labs/sjaiswal/genomes/mm10/STAR_index/" #path to STAR genome index

annotation_gtf="/labs/sjaiswal/genomes/mm10/gencode.vM25.annotation.gtf" #path to annotation GTF

reference_genome="/labs/sjaiswal/genomes/mm10/GRCm38.primary_assembly.genome.fa" #location of reference genome used in alignment

annotation_refFlat="/labs/sjaiswal/genomes/mm10/gencode.vM25.annotation.fixed.refFlat.txt" ## Created from ensembl GTF using gtfToGenePred -genePredExt -geneNameAsName2 from bioconda

gene_key="/labs/sjaiswal/genomes/mm10/gencode.vM25.annotation.genekey.txt" #path to annotation GTF
