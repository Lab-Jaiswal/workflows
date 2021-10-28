library(magrittr)
library(tidyverse)

command_args <- commandArgs(trailingOnly = TRUE) #Read in arguments from command line
gtf_file <- command_args[1] #Assign first argument to gtf_file

read_gtf <- read_tsv(gtf_file, col_names = FALSE, skip = 5) #Read in GTF file
gtf_annot_split <- str_split(read_gtf$X9, ";") #Split GTF annotation column by ; character
gtf_ids <- map_chr(gtf_annot_split, magrittr::extract, 1) %>% str_replace_all("gene_id ", "") %>% str_replace_all('\\\"', '') #Extract first field which will always be gene_id and strip out characters that are not part of Ensembl ID
gtf_name_locs <- map(gtf_annot_split, str_detect, "gene_name") %>% map_int(which) #Find out which field is gene name - this location will be variable for each line 
gtf_names <- map2(gtf_annot_split, gtf_name_locs, magrittr::extract) %>% str_replace_all("gene_name ", "") %>% str_replace_all('\\\"', '') %>% str_replace_all(" ", "") #Extract field corresponding to gene_name and strip characters which are not part of gene symbol

gtf_names_df <- tibble(Ensembl_ID = gtf_ids, Symbol = gtf_names) %>% distinct #Combined gene IDs and symbol columns into a tibble and only keep unique rows
gtf_basename <- str_replace(gtf_file, ".gtf", "")
genekey_filename <- str_c(gtf_basename, ".genekey.txt")
write_tsv(gtf_names_df, genekey_filename) #Write gene key tibble to disk as a text file
