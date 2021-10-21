library(rtracklayer)
library(magrittr)
library(tidyverse)

gtf_ranges <- readGFF("/labs/sjaiswal/genomes/mm10/gencode.vM25.annotation.gtf")
gtf_genes_rrna <- filter(gtf_ranges, type == "gene" & gene_type == "rRNA")
gtf_genes_rrna_select <- select(gtf_genes_rrna, seqid, start, end, strand, gene_id)


