#install.packages("devtools")
#library(devtools)
devtools::install_github("Lab-Jaiswal/tidyUkBioBank", force=TRUE)
library(readxl)
library(magrittr)
library(tidyverse)
#devtools::install_github('chrchang/plink-ng', subdir='2.0/pgenlibr')
library(tidyUkBioBank)
library(ukbtools)
library(lubridate)

#Hanscombe KB, Coleman J, Traylor M, Lewis CM (e-print 158113). “ukbtoo
#UK Biobank data.” _bioRxiv_. <URL: https://doi.org/10.1101/158113>.

#read in the pruned ukbiobank dataframe (unnecessary columns removed)
ukb_data  <- read_rds("/home/maurertm/groups/maurertm/outputs/R_objects/ukb_pruned.rds")
#filter the dataframe to only contain eids that have a cram or crai file on DNANexus
eid_list <- read.table("/oak/stanford/groups/sjaiswal/maurertm/R_scripts/eids_with_cram_files", header=FALSE)
eids_with_cram_files <- eid_list$V1
ukb_data_with_cram_files <- filter(ukb_data, is_in(eid, eids_with_cram_files))

#get 5000 youngest 
ukb_data_eids_dates <- ukb_data_with_cram_files %>% select(eid, age_when_attended_assessment_centre_f21003_0_0)
ranked_ages <- ukb_data_eids_dates[with(ukb_data_eids_dates, order(as.numeric(as.character(ukb_data_eids_dates$age_when_attended_assessment_centre_f21003_0_0)))), ]
top_5000 <- head(ranked_ages, 5000) %>% pull(eid)
ukb_data_youngest <- filter(ukb_data, is_in(eid, top_5000))

#filter out any one who reported having cancer
cancer_colnames <- ukb_data %>%  select(contains("cancer", ignore.case = TRUE)) %>% colnames
cancer_colnames <- Filter(function(x) !any(grepl("noncancer", x)), cancer_colnames)

#get ICD codes for malignant neoplasms, autoimmune disorders (lupus, ra, sjrogens, crohns, uc, ibd, and autoimmune thyroiditis), and hx of DMARD rx
no_cancer <- filter(ukb_data_youngest, is.na(reported_occurrences_of_cancer_f40009_0_0))
ukb_data_youngest$reported_occurrences_of_cancer_f40009_0_0 %>% unique
#cancer is not included in self reported chart
cancer_10 <- ukb_icd_keyword("malignant", icd.version = 10) %>% as_tibble
cancer_10 %>% print(n=200)
cancer_10_modified <- cancer_10[c(1:135, 139:145, 158:167),] 
cancer_10_codes <- cancer_10_modified %>% pull(code)
cancer_9 <- ukb_icd_keyword("malignant", icd.version = 9) %>% as_tibble
cancer_9 %>% print(n=500)
cancer_9_modified <- cancer_9[c(1, 3:380, 382:397),] 
cancer_9_codes <- cancer_9_modified %>% pull(code)

lupus_10 <- ukb_icd_keyword("lupus", icd.version = 10) %>% as_tibble
lupus_10_modified <- lupus_10[c(1, 5, 7:9),]
lupus_10_codes <- lupus_10_modified %>% pull(code)
lupus_9 <- ukb_icd_keyword("lupus", icd.version = 9) %>% as_tibble
lupus_9_codes <- lupus_9 %>% pull(code)

RA_10 <- ukb_icd_keyword("rheum", icd.version = 10) %>% as_tibble
RA_10_modified <- RA_10[c(64:108, 131:163),]
RA_10_codes <- RA_10_modified %>% pull(code)
RA_9 <- ukb_icd_keyword("rheum", icd.version = 9) %>% as_tibble
RA_9_modified <- RA_9[c(6:19, 131:163),]
RA_9_codes <- RA_9_modified %>% pull(code)

antirheum_10 <- ukb_icd_keyword("antirheum", icd.version = 10) %>% as_tibble
antirheum_10_modified <- antirheum_10[c(2, 7),]
antirheum_10_codes <- antirheum_10_modified %>% pull(code)
antirheum_9 <- ukb_icd_keyword("antirheum", icd.version = 9) %>% as_tibble
antirheum_9_modified <- antirheum_9[c(18),]
antirheum_9_codes <- antirheum_9_modified %>% pull(code)

SJ_10 <- ukb_icd_keyword("sicca", icd.version = 10) %>% as_tibble() #code: M350
SJ_10_codes <- SJ_10 %>% pull(code)
SJ_9 <- ukb_icd_keyword("sicca", icd.version = 9) %>% as_tibble()  #code: 7102
SJ_9_codes <- SJ_9 %>% pull(code)

Crohn_10 <- ukb_icd_keyword("Crohn", icd.version = 10) %>% as_tibble()
Crohn_10_codes <- Crohn_10 %>% pull(code)
Crohn_9 <- ukb_icd_keyword("enteritis", icd.version = 9) %>% as_tibble()      
Crohn_9_modified <- Crohn_9[c(7:16),]
Crohn_9_codes <- Crohn_9_modified %>% pull(code)

UC_10 <- ukb_icd_keyword("Colitis", icd.version = 10) %>% as_tibble() 
UC_10_modified <- UC_10[c(6:8, 10:11, 19:40), ]
UC_10_codes <- UC_10_modified %>% pull(code)

ukb_icd_keyword("ulerative", icd.version = 9) %>% as_tibble() %>% print(n=100)
ukb_icd_keyword("colitis", icd.version = 9) %>% as_tibble() %>% print(n=100)
#Could not find anything. Did this diagnosis have another name in the icd9? Both colitis and ulerative do not bring up UC.

Thy_10 <- ukb_icd_keyword("thyroiditis", icd.version = 10) %>% as_tibble()
Thy_10_modified <- Thy_10[c(1:2, 4:5, 7:8), "code"] 
Thy_10_codes <- Thy_10_modified %>% pull(code)
Thy_9 <- ukb_icd_keyword("thyroiditis", icd.version = 9) %>% as_tibble() 
Thy_9_modified <- Thy_9[c(1, 8, 2), "code"] 
Thy_9_codes <- Thy_9_modified %>% pull(code)


excluded_icds <- cancer_10_codes %>% append(cancer_9_codes) %>% append(lupus_10_codes) %>% append(lupus_9_codes) %>%
                append(RA_10_codes) %>% append(RA_9_codes) %>% append(antirheum_10_codes) %>% append(antirheum_9_codes) %>%
                append(SJ_10_codes) %>% append(SJ_9_codes) %>% append(Crohn_10_codes) %>% append(Crohn_9_codes)  %>% append(UC_10_codes) %>%
                append(Thy_10_codes) %>% append(Thy_9_codes)
excluded_icds <- excluded_icds[!is.na(excluded_icds)]

#in the 5000 youngest with no self reported cancer dx, filter out anyone with the above icd codes (in excluded_icds)
indiv_with_excluded_icds <- individuals_with_disease(excluded_icds, no_cancer) 
excluded_icd_eids <- indiv_with_excluded_icds$combined

top_5000_without_predis <- filter(no_cancer, !is_in(eid, excluded_icd_eids))
                          
#get 1000 youngest
top_1000_without_predis <- head(top_5000_without_predis, 1000) %>% pull(eid)

#write output
write_rds(top_1000_without_predis, "/oak/stanford/groups/sjaiswal/maurertm/R_scripts/top_1000_without_predis_with_crams.rds")
write_lines(top_1000_without_predis, "/oak/stanford/groups/sjaiswal/maurertm/R_scripts/top_1000_without_predis_with_crams")
