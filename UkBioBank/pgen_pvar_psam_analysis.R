library(maftools)
library(GenomicRanges)
library(readxl)
library(magrittr)
library(tidyverse)
#devtools::install_github('chrchang/plink-ng', subdir='2.0/pgenlibr')
library(pgenlibr)
devtools::install_github("Lab-Jaiswal/tools4ukbb", force=TRUE)
library(tools4ukbb)
library(ukbtools)
library(lubridate)

sequnce <- seq(1:22) %>% append(c("X", "Y"))
chr_list<- sprintf("c%s", sequnce) 

pfile <- '/oak/stanford/groups/sjaiswal/maurertm/22828_pgens/'
directory <- '/oak/stanford/groups/sjaiswal/maurertm/22828_pgens/'

make_pvar <- function(chr){
  chrs_pvar <- paste0(chr, "_b0_v3.pvar", sep="")
  pvar_file_name <- list.files(directory, pattern = chrs_pvar, full.names = TRUE)
  pvar <- NewPvar(pvar_file_name)
}

make_pgen <- function(chr){
  chrs_pgen <- paste0(chr, "_b0_v3.pgen", sep="")
  chrs_pvar <- paste0(chr, "_b0_v3.pvar", sep="")
   
  pgen_file_name <- list.files(directory, pattern = chrs_pgen, full.names = TRUE)
  pvar_file_name <- list.files(directory, pattern = chrs_pvar, full.names = TRUE)
  pvar <- NewPvar(pvar_file_name)
  pgen <- NewPgen(pgen_file_name, pvar=pvar)
}

make_psam <- function(chr){
  chrs_psam <- paste0(chr, "_b0_v3.psam", sep="")
  psam_file_name <- list.files(directory, pattern = chrs_psam, full.names = TRUE)
  psam <- read_tsv(psam_file_name)
}


get_variant <- function(id_num){
  GetVariantId(pvars_current, id_num)
}

#example- make a pvar, psam, and pgen of chromosome 19:
  #pvars_current <- make_pvar("c19")
  #psams_current <- make_psam("c19")
  #pgen_current <- make_pgen("c19")

#example- make a pvar, psam, and pgen of every chromosome
  #for (i in chr_list){
   # pvar <- paste0(i, "_pvar", sep = "")
   # pgen <- paste0(i, "_pgen", sep = "")
   # psam <- paste0(i, "_psam", sep = "")
   # assign(pvar, make_pvar(i))
   # assign(pgen, make_pgen(i))
   # assign(psam, make_psam(i))
  #}    

#example- get the APOE2 and APOE4 variants
#variants_apoe <- filter(variants_df, is_in(value, c("rs429358", "rs7412")))
#output:
   # 1 rs429358 1569023
   # 2 rs7412   1569031

#example- get dataframe with eid, and APOE2 and APOE4 genotypes:
  #geno_mat_orig <- ReadList(pgen_current, variants_apoe$index, meanimpute=F)  
  #geno_mat_flip <- round(2 - geno_mat_orig)
  #geno_mat <- as_tibble(geno_mat_flip)
  #geno_mat$eid <- psams_current$IID
  #colnames(geno_mat) <- c("APOE4", "APOE2", "eid")
