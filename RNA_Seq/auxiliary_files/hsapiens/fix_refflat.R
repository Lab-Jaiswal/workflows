library(magrittr)
library(tidyverse)

command_args <- commandArgs(trailingOnly = TRUE) #Read in arguments from command line
refFlat_name <- command_args[1] #Assign first argument to refFlat_name

gtf_read <- read.delim(refFlat_name, header = FALSE, stringsAsFactors = FALSE) %>% #Read in refFlat file created by gtfToGenePred
    select(V1, V12, V2:10) #Reorder columns to comply with genePred specification

file_basename <- str_replace(refFlat_name, "\\.refFlat\\.txt", "") #Delete file extension
fixed_filename <- str_c(file_basename, ".fixed.refFlat.txt") #Add fixed to name, as well as file extension
write_tsv(gtf_read, fixed_filename, col_names = FALSE) #Write reordered tibble back to disk as text file
