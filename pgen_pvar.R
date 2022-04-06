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


pfile <- '/oak/stanford/projects/ukbb/genotypes/pgen_app13721_v3/'
directory <- "/oak/stanford/projects/ukbb/genotypes/pgen_app13721_v3"

#sample_names <- list.files(directory, pattern = "*.mac1.pvar$") %>% str_remove_all(".pvar*$")
sequnce <- seq(1:22) %>% append(c("X", "Y"))
chr_list<- sprintf("chr%s", sequnce) 

make_pvar_old <- function(chr){
  chrs_pvar <- paste0(chr, "_v3.mac1.pvar", sep="")
  #chrs_pvar <- paste0("ukb22828_", "chr1", "_b0_v3.pvar", sep="")
  #chrs_pvar <- paste0(chr, "_b0_v3.pvar", sep="")
  
  
  pvar_file_name <- list.files(directory, pattern = chrs_pvar, full.names = TRUE)
  pvar <- NewPvar(pvar_file_name)
}

make_pgen_old <- function(chr){
  #chrs_pgen <- paste0(chr, "_b0_v3.pgen", sep="")
  #chrs_pvar <- paste0(chr, "_b0_v3.pvar", sep="")
  chrs_pgen <- paste0(chr, "_v3.mac1.pgen", sep="")
  chrs_pvar <- paste0(chr, "_v3.mac1.pvar", sep="")
  
  pgen_file_name <- list.files(directory, pattern = chrs_pgen, full.names = TRUE)
  pvar_file_name <- list.files(directory, pattern = chrs_pvar, full.names = TRUE)
  pvar <- NewPvar(pvar_file_name)
  pgen <- NewPgen(pgen_file_name, pvar=pvar)
}

make_psam_old <- function(chr){
  chrs_psam <- paste0(chr, "_v3.mac1.psam", sep="")
  #chrs_psam <- paste0(chr, "_b0_v3.psam", sep="")
  psam_file_name <- list.files(directory, pattern = chrs_psam, full.names = TRUE)
  psam <- read_tsv(psam_file_name)
}

pvars_old <- make_pvar_old("chr19")
psams_old <- make_psam_old("chr19")
pgen_old <- make_pgen_old("chr19")



sequnce <- seq(1:22) %>% append(c("X", "Y"))
chr_list<- sprintf("c%s", sequnce) 

pfile <- '/oak/stanford/groups/sjaiswal/maurertm/22828_pgens/'
#directory <- "/oak/stanford/projects/ukbb/genotypes/exome.v1/vcf/ukb_spb_exm_vcf/59"
directory <- '/oak/stanford/groups/sjaiswal/maurertm/22828_pgens/'

make_pvar <- function(chr){
  #chrs_pvar <- paste0(chr, "_v3.mac1.pvar", sep="")
  #chrs_pvar <- paste0("ukb22828_", "chr1", "_b0_v3.pvar", sep="")
  chrs_pvar <- paste0(chr, "_b0_v3.pvar", sep="")
  
  
  pvar_file_name <- list.files(directory, pattern = chrs_pvar, full.names = TRUE)
  pvar <- NewPvar(pvar_file_name)
}

make_pgen <- function(chr){
  chrs_pgen <- paste0(chr, "_b0_v3.pgen", sep="")
  chrs_pvar <- paste0(chr, "_b0_v3.pvar", sep="")
  #chrs_pgen <- paste0(chr, "_v3.mac1.pgen", sep="")
  #chrs_pvar <- paste0(chr, "_v3.mac1.pvar", sep="")
  
  pgen_file_name <- list.files(directory, pattern = chrs_pgen, full.names = TRUE)
  pvar_file_name <- list.files(directory, pattern = chrs_pvar, full.names = TRUE)
  pvar <- NewPvar(pvar_file_name)
  pgen <- NewPgen(pgen_file_name, pvar=pvar)
}

make_psam <- function(chr){
  #chrs_psam <- paste0(chr, "_v3.mac1.psam", sep="")
  chrs_psam <- paste0(chr, "_b0_v3.psam", sep="")
  psam_file_name <- list.files(directory, pattern = chrs_psam, full.names = TRUE)
  psam <- read_tsv(psam_file_name)
}

pvars_current <- make_pvar("c19")
psams_current <- make_psam("c19")
pgen_current <- make_pgen("c19")



chr_list_edited <- chr_list[c(1:5, 7:20, 22)]
#no X, Y, 6, 

for (i in chr_list_edited){
  pvar <- paste0(i, "_pvar", sep = "")
  pgen <- paste0(i, "_pgen", sep = "")
  psam <- paste0(i, "_psam", sep = "")
  assign(pvar, make_pvar(i))
  assign(pgen, make_pgen(i))
  assign(psam, make_psam(i))
}    


geno_mat_old <- ReadList(pgen_old, 1:10, meanimpute=F)
geno_mat_current <- ReadList(pgen_current, 1:10, meanimpute=F)

get_variant <- function(id_num){
  GetVariantId(pvars_current, id_num)
}

variants <- map(1:2074989, get_variant) %>% unlist
variants_df <- as_tibble(variants)
variants_df$index <- rownames(variants_df) %>% as.integer
#variants_df[grep("^rs7412$", variants_df$value),]
#variants_df[grep("^19:45412079", variants_df$value),]
variants_apoe <- filter(variants_df, is_in(value, c("rs429358", "rs7412")))
#1 rs429358 1569023
#2 rs7412   1569031

for(i in 1:length(variants)){
  vector[i] <- grepl(variants[i],"^19:45411941")
}





#GetVariantCt(pvar)
#GetRawSampleCt(pgen)
geno_mat_orig <- ReadList(pgen_current, variants_apoe$index, meanimpute=F)  
geno_mat_flip <- round(2 - geno_mat_orig)
geno_mat <- as_tibble(geno_mat_flip)
geno_mat$eid <- psams_current$IID
colnames(geno_mat) <- c("APOE4", "APOE2", "eid")


ukb_data  <- read_rds("/oak/stanford/groups/sjaiswal/maurertm/processing/ukb_pruned.rds")

#Alzheimer's
#ICD10 and ICD9
AD10 <- ukb_icd_keyword("alzheim", icd.version = 10) %>% as_tibble()
AD9 <-  ukb_icd_keyword("alzheim", icd.version = 9) %>% as_tibble()
AD_icd_10 <- ukb_icd_keyword("alzheim", icd.version = 10) %>% as_tibble() %>% select(code) %>% as.list()  
AD_icd_9 <- c("3310")
#Combined ICDs
AD_UF <- rbind(AD10, AD9)
AD_icd <- c(AD_icd_10, AD_icd_9) %>% flatten() %>% unlist() 
#Self Reported
#filter(self_reported_codes, grepl("alz", meaning)) %>% as_tibble() #Nothing applicable here   
#Collect the number of people with icd9 and icd10 codes consistent with AD
AD <- individuals_with_disease(AD_icd, ukb_data)
#Column Specific Information
alzheimers_subset<- filter(ukb_data, !is.na(date_of_alzheimers_disease_report_f42020_0_0)) %>% select(eid, date_of_alzheimers_disease_report_f42020_0_0)
eid_AD_subset <- alzheimers_subset$eid

#Combine column specific and ICD information
combined <- append(AD$icd9_and_icd10, eid_AD_subset) %>% unique()

#add Hx column
AD_df <- ukb_data[ukb_data$eid %in% combined, ]
nonAD_df <- ukb_data[!ukb_data$eid %in% combined, ]
AD_df$AD <- 1
nonAD_df$AD <- 0
AD_combined <- rbind(AD_df, nonAD_df)

#Remove those in the ukb_data bank not in the APOE bank
ukb_eid <- ukb_data$eid
apoe_eid <- geno_mat$eid
genotyped <- intersect(apoe_eid, ukb_eid)
AD_filtered <- filter(AD_combined, eid %in% genotyped)

#Get DOB and age
DOB <- date_of_birth(AD_filtered)
DOB$Current_age = as.numeric(difftime(Sys.Date(),DOB$DOB, units = "weeks"))/52.25

#Merge AD_filtered and DOB on 
AD_select_AD <- AD_filtered %>% select(AD, eid)
AD_AD_eid <- AD_select_AD %>% merge(DOB, by = 'eid')

#Merge AD_AD_eid and APOE results
Final_APOE_df <- geno_mat %>% merge(AD_AD_eid, by = 'eid') %>% select(-DOB)
colnames(Final_APOE_df) <- c("eid", "APOE4", "APOE2", "AD", "Age")

#glm analysis
ad_glm_APOE4 <- glm(AD ~ APOE4 + Age, family = "binomial", data = Final_APOE_df) # Should be a positive effect
#Coefficients:
#  (Intercept)        APOE4          Age
#-21.8643       1.2983       0.2133

#Degrees of Freedom: 487202 Total (i.e. Null);  487200 Residual
#Null Deviance:      34370
#Residual Deviance: 28670        AIC: 28680
ad_glm_APOE2 <- glm(AD ~ APOE2 + Age, family = "binomial", data = Final_APOE_df) # Should be a negative effect
#Coefficients:
#(Intercept)        APOE2          Age
#-20.9002      -0.8318       0.2107

#Degrees of Freedom: 487202 Total (i.e. Null);  487200 Residual
#Null Deviance:      34370
#Residual Deviance: 30360        AIC: 30360