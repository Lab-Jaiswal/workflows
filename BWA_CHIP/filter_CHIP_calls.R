library(maftools)
library(GenomicRanges)
library(readxl)
library(magrittr)
library(tidyverse)

#specificiations
#tumor_f 0.02-0.35
#DP -ge 20
#filter pass, weak
#t_alt_count -ge 5
#as_sb_table can't contain 0 for any 4 values
# for ^|0 or , 0 or [0 or
#need whitelist table


# get benchmark at this point

#get rid of vcf_f2r1_alternate2 if not NA
#get rid of maf_mmq_alternate2 if not NA (or any col with alternate 2 if not na)
#remove if n of 4 is the same
#paste hugo and protein chaige and if 4+ have it remove it
#and have version of ^ with cutoff -ge 2 (or -lt 2)
#drop SF3A1, GATA1, GATA3, PTEN, SF1, STAG1, IKZF2, IKZF3, PDSS2, LUC7L2, JAK1, JAK3, GNA13, KMT2A, KMT2D, and CSF1R


black_list <- c("SF3A1", "GATA1", "GATA3", "PTEN", "SF1", "STAG1", "IKZF2", "IKZF3", "PDSS2", "LUC7L2", "JAK1", "JAK3", "GNA13", "KMT2A", "KMT2D", "CSF1R")
chip <- "/Users/maurertm/Downloads/whitelisted.csv"
chip <- read_csv(chip) 
chip$hugo_protein <- str_c(chip$Hugo_Symbol, "_", chip$Protein_Change)
chip_filtered <- chip %>% filter(tumor_f <= 0.35) %>% filter(tumor_f >= 0.02) %>%
  filter(DP >= 20) %>% filter(t_alt_count >= 5) %>% filter(vcf_FILTER == "PASS" | vcf_FILTER == "weak_evidence") %>%
  #filter(!is.na(vcf_f1r2_alternate2)) %>% filter(!is.na(vcf_f2r1_alternate2)) %>% filter(!is.na(maf_mmq_alternate2)) 
  #^commented out filters out multi-alelic variants
  duplicates <- chip_filtered %>% group_by(hugo_protein) %>% tally() %>% arrange(desc(n))
duplicates <- filter(duplicates, n >= 2)
duplicate_list <- duplicates$hugo_protein

no_duplicates <- dplyr::filter(chip_filtered, !magrittr::is_in(hugo_protein, duplicate_list))
no_blacklisted <- dplyr::filter(no_duplicates, !magrittr::is_in(Hugo_Symbol, black_list))
final <- dplyr::filter(no_blacklisted, !str_detect(AS_SB_TABLE, "\\[0|\\|0|, 0"))
