library(magrittr)
library(rtracklayer)
library(tidyverse)

command_args <- commandArgs(trailingOnly = TRUE) #Read in arguments from command line
gtf_file <- command_args[1] #Assign first argument to gtf_file

read_gtf <- readGFF(gtf_file) %>% filter(type == "gene") #Read in GTF file
gtf_names_df <- select(read_gtf, gene_id, gene_name) %>% set_colnames(c("Ensembl_ID", "Symbol")) %>% distinct
gtf_basename <- str_replace(gtf_file, ".gtf", "")
genekey_filename <- str_c(gtf_basename, ".genekey.txt")
write_tsv(gtf_names_df, genekey_filename) #Write gene key tibble to disk as a text file
