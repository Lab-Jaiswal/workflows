library(maftools)
library(readxl)
library(magrittr)
library(tidyverse)
library(dplyr)
library(purrr)
library(openxlsx)

string_split_slash <- function(df, split_number){
  slash_split <- str_split_fixed(df$Variant_calling_mutations, "/", split_number) %>% as_tibble()
  slash_split_joined <- bind_cols(df, slash_split) %>% select(-Variant_calling_mutations)
  slash_split_long <- pivot_longer(slash_split_joined, -Gene_name, names_to = "variants", values_to = "Variant_types") %>% select(-variants)
  slash_split_long[slash_split_long==""] <-NA
  filter(slash_split_long, !is.na(Variant_types))
}

set_columns <- function(df, start_position, end_position, WL, MR, note){
  df$start <- start_position
  df$end <- end_position
  df$whitelist <- WL
  df$manual_review <- MR
  df$notes <- note
  df
}

string_split_comma <- function(df, split_number) {
  split_comma <- str_split_fixed(df$Variant_types, ",", split_number) %>% as_tibble()
  split_comma_joined <- bind_cols(df, split_comma) %>% select(-Variant_types)
  variant_type_long <- pivot_longer(split_comma_joined, -Gene_name, names_to = "variants", values_to = "Variant_types") %>% select(-variants)
  variant_type_long[variant_type_long==""] <-NA
  variant_cleaned <- filter(variant_type_long, !is.na(Variant_types))
  variant_no_qoutes <- as.data.frame(sapply(variant_cleaned, function(x) gsub("\"", "", x)))
  types <- c("Frameshift", "nonsense", "splice-site")
    if (nrow(filter(variant_no_qoutes, is_in(Variant_types, types))) > 1) {
      non_missense <- filter(variant_no_qoutes, is_in(Variant_types, types))
      non_missense$missense_type <- NA
      missense <- filter(variant_no_qoutes, !is_in(Variant_types, types))
      missense$missense_type <- missense$Variant_types
      missense$Variant_types <- "missense"
      variants <- bind_rows(missense, non_missense)
    } else {
      variant_no_qoutes$missense_type <- variant_no_qoutes$Variant_types
      variant_no_qoutes$Variant_types <-  "missense"
      variants <- variant_no_qoutes
    }
  
  
}

to_char <- function(mydata){
  mydata <- mydata %>% mutate_all(as.character)
  
}

command_args <- commandArgs(trailingOnly = TRUE)
variantWL_location <- command_args[1]
output <- command_args[2]
variantWL <- read_excel(variantWL_location, col_names = T)

variantWL_rmAccession <- variantWL %>% select(-Accession, -...4, -...5, -...6, -...7, -...8)
names(variantWL_rmAccession) <- c("Gene_name", "Variant_calling_mutations")
contains_pdot <- c("ASXL1", "ASXL2", "CBL","CBLB", "NPM1", "TET2")
ZNF_domain <- c("GATA3")
contains_parentheses <- c("GNAS", "PRPF8")
ins_del <- c("JAK2", "KIT", "KDM6A", "MPL", "CREBBP", "DNMT3A")
unique <- c("PPM1D", "TP53", "EP300", "FLT3")

complicated <- list(contains_pdot, ZNF_domain, contains_parentheses, ins_del, unique) %>% flatten() %>% unlist()

WL_filtered <- filter(variantWL_rmAccession, !is_in(Gene_name, complicated))
filtered_split_slash <- string_split_slash(WL_filtered, 4)
filtered_split_comma <- string_split_comma(filtered_split_slash, 600) 

filtered_expanded <- set_columns(filtered_split_comma, NA, NA, 1, 0, NA)
filtered_final <- filtered_expanded

WL_complicated <- filter(variantWL_rmAccession, is_in(Gene_name, complicated))

#################Contains positional information###############################
############"ASXL1", "ASXL2", "CBL","CBLB", "NPM1", "TET2")#####################
contains_position <- WL_complicated[c(1:4, 15, 18),]
contains_position_split <- string_split_slash(contains_position, 5)

#ASXL1/ASXL2
ASXL_358 <- contains_position_split[1:3, ]
ASXL_380 <- contains_position_split[4:6, ]
ASXL_358$missense_type <- NA
ASXL_380$missense_type <- NA
ASXL_358_expanded <- set_columns(ASXL_358, "358", "1541", 1, 0, NA)
ASXL_380_expanded <- set_columns(ASXL_380, "380", "1435", 1, 0, NA)
ASXL_expanded <- rbind(ASXL_358_expanded, ASXL_380_expanded)
ASXL_expanded$Variant_types <- c("Frameshift", "nonsense", "splice-site", "Frameshift", "nonsense", "splice-site")
ASXL_final <- ASXL_expanded

#CBL/CBLB
CBL_CBLB <- contains_position_split[7:8, ]
CBL_CBLB$missense_type <- NA
CBL_CBLB_expanded <- set_columns(CBL_CBLB, c("381", "372"), c("421", "412"), 0, 1, c("RING finger missense p.381-421", "RING finger missense p.372-412"))
CBL_CBLB_expanded$Variant_types <- "missense"
CBL_CBLB_final <- CBL_CBLB_expanded

#NMP1
NPM1 <- contains_position_split[c(9, 9, 9, 9), ]
NPM1$missense_type <- NA
NPM1_expanded <- set_columns(NPM1, c("859", "860", "862", "863"), c("860", "861", "863","864"), 0, 1, "Frameshift p.W288fs (insertion at c.859_860, 860_861, 862_863, 863_864)")
NPM1_expanded$Variant_types <- "Frameshift"
NPM1_final <- NPM1_expanded

#TET2
TET2 <- contains_position_split[10:12, ]
TET2_long <- string_split_comma(TET2, 600)
TET2_expanded <- set_columns(TET2_long, NA, NA, 1, 0, NA)

TET2_manual_1104 <- TET2_expanded[c(1, 268:270),]
TET2_manual_1843 <- TET2_expanded[c(1, 268:270),]

TET2_expanded <- TET2_expanded[-c(1, 268:270),]

TET2_manual_1104_expanded <- set_columns(TET2_manual_1104, "1104", "1843", 0, 1, "in catalytic domains (p.1104-1481 and 1843-2002)")
TET2_manual_1104_expanded[1,3] <- NA
TET2_manual_1843_expanded <- set_columns(TET2_manual_1843, "1843", "2002", 0, 1, "in catalytic domains (p.1104-1481 and 1843-2002)")
TET2_manual_1843_expanded[1,3] <- NA

TET2_final <- rbind(TET2_manual_1104_expanded, TET2_manual_1843_expanded, TET2_expanded)

############################Others##############################################
###############################################################################

#CREBBP
CREBBP <- WL_complicated %>% filter(Gene_name == "CREBBP")
CREBBP_split_slash <- string_split_slash(CREBBP, 3)
CREBBP_split_comma <- string_split_comma(CREBBP_split_slash, 15) 
CREBBP_expanded <- set_columns(CREBBP_split_comma, NA, NA, 1, 0, NA)
CREBBP_expanded[13,] <- c("CREBBP", "missense", "S1680", NA, NA, 0, 1, "deletion")
CREBBP_final <- CREBBP_expanded


#DNMT3A
DNMT3A <- WL_complicated %>% filter(Gene_name == "DNMT3A")
DNMT3A_split_slash <- string_split_slash(DNMT3A, 3)
DNMT3A_split_comma <- string_split_comma(DNMT3A_split_slash, 590)
DNMT3A_expanded <- set_columns(DNMT3A_split_comma, NA, NA, 1, 0, NA)

DNMT3A_long_259 <- DNMT3A_expanded[c(259, 259),]
DNMT3A_long_259[1,3] <- "D702G"
DNMT3A_long_259[2, 3] <-"D702E"

DNMT3A_long_247 <- DNMT3A_expanded[c(247, 247),]
DNMT3A_long_247[1,3] <- "W698S"
DNMT3A_long_247[2,3] <- "G699R"

DNMT3A_long_243 <- DNMT3A_expanded[c(243, 243),]
DNMT3A_long_243[1,3] <-  "I695F"
DNMT3A_long_243[2,3] <- "I695T"
DNMT3A_long_doubles <- rbind(DNMT3A_long_243, DNMT3A_long_247, DNMT3A_long_259)

DNMT3A_long_del <- DNMT3A_expanded[c(313, 363),]
DNMT3A_long_del[1, 8] <- "deletion"
DNMT3A_long_del[2,8] <- "deletion"

DNMT3A_expanded <- DNMT3A_expanded[-259,]
DNMT3A_expanded <- DNMT3A_expanded[-247,]
DNMT3A_expanded <- DNMT3A_expanded[-243,]
DNMT3A_expanded <- DNMT3A_expanded[-313,]
DNMT3A_expanded <- DNMT3A_expanded[-363,]

DNMT3A_final <- rbind(DNMT3A_expanded, DNMT3A_long_doubles, DNMT3A_long_del)

#EP300
EP300 <- WL_complicated %>% filter(Gene_name == "EP300")
EP300_split_slash <- string_split_slash(EP300, 3)
EP300_split_comma <- string_split_comma(EP300_split_slash, 10)

EP300_expanded <- set_columns(EP300_split_comma, NA, NA, 1, 0, NA)

EP300_manual <- EP300_expanded[c(1,1),]
EP300_manual[1,] <- c("EP300", "missense", "1148", NA, NA, 0, 1, "VF1148_1149del")
EP300_manual[2,] <- c("EP300", "missense", "1149", NA, NA, 0, 1, "VF1148_1149del")

EP300_final <- rbind(EP300_manual, EP300_expanded)
EP300_final <- EP300_final[-3,]

#FLT3
#fix FY590-591GD
FLT3 <- WL_complicated %>% filter(Gene_name == "FLT3")
FLT3_split_slash <- string_split_slash(FLT3, 3)
FLT3_split_comma <- string_split_comma(FLT3_split_slash, 10)

FLT3_expanded <- set_columns(FLT3_split_comma, NA, NA, 1, 0, NA)

FLT3_manual <- FLT3_expanded[c(4,8),]
FLT3_manual[1,] <- c("FLT3", "missense", "590", NA, NA, 0, 1, "FY590-591GD")
FLT3_manual[2,] <- c("FLT3", "missense", "591", NA, NA, 0, 1, "FY590-591GD")
FLT3_manual[3,] <- c("FLT3", "missense", "835", NA, NA, 0, 1, "deletion")

FLT3_expanded <- FLT3_expanded[-c(4,8),]

FLT3_final <- rbind(FLT3_manual, FLT3_expanded)

#GATA3
GATA3 <- WL_complicated %>% filter(Gene_name == "GATA3")
GATA3_split_slash <- string_split_slash(GATA3, 4)
GATA3_split_comma <- string_split_comma(GATA3_split_slash, 10)
GATA3_expanded <- set_columns(GATA3_split_comma, NA, NA, c(1, 0, 0, 0, 0, 1, 1), c(0, 1, 1, 1, 1, 0, 0), c("ZNF_domain", NA, NA, NA, NA, "ZNF_domain", "ZNF_domain"))
GATA3_expanded[1, 2:3] <- c("splice-site", NA)
GATA3_expanded[1, 3] <- NA
GATA3_final <- GATA3_expanded

#GNAS
#ask sidd about this
GNAS <- WL_complicated %>% filter(Gene_name == "GNAS")
GNAS$Variant_calling_mutations <- c("R201G, R201S, R201C, R201H, R201L, Q227K, Q227R, Q227L, Q227H, R374C")
GNAS_split_slash <- string_split_slash(GNAS, 4)
GNAS_split_comma <- string_split_comma(GNAS_split_slash, 10)
GNAS_expanded <- set_columns(GNAS_split_comma, NA, NA, 1, 0, c("R201(844)G", "R201(844)S",  "R201(844)C",  "R201(844)H",  "R201(844)L",  "Q227(870)K",  "Q227(870)R",  "Q227(870)L",  "Q227(870)H",  "R374(1017)C"))
GNAS_final <- GNAS_expanded

#JAK2
#ask sidd
JAK2 <- WL_complicated %>% filter(Gene_name == "JAK2")
JAK2_split_slash <- string_split_slash(JAK2, 4)
JAK2_split_comma <- string_split_comma(JAK2_split_slash, 12)



JAK2_manual <- JAK2_split_comma[c(13:22),]
JAK2_manual$missense_type <- c("537", "538", "539", "540", "541", "542", "543", "544", "546", "547")
JAK2_manual_expanded <- set_columns(JAK2_manual, NA, NA, 0, 1, c("del/ins537-539L", "del/ins537-539L_del/ins538-539L", "del/ins537-539L_del/ins538-539L", "del/ins540-543MK_del/ins540-544MK","del/ins540-543MK_del/ins540-544MK_del/ins541-543K", #41
                                                                 "del/ins540-543MK_del/ins540-544MK_del/ins541-543K_del542-543", "del/ins540-543MK_del/ins540-544MK_del/ins541-543K_del542-543_del543-544", "del/ins540-544MK_del543-544", #44
                                                                 "ins11546-547", "ins11546-547"))

JAK2_long <- JAK2_split_comma[-c(12:22),]
JAK2_long_expanded <- set_columns(JAK2_long, NA, NA, 1, 0, NA)

JAK2_final <- rbind(JAK2_long_expanded, JAK2_manual_expanded)

#KDM6A
KDM6A <- WL_complicated %>% filter(Gene_name == "KDM6A")
KDM6A_split_slash <- string_split_slash(KDM6A, 4)
KDM6A_split_comma <- string_split_comma(KDM6A_split_slash, 12)
KDM6A_expanded <- set_columns(KDM6A_split_comma, NA, NA, 1, 0, NA)
KDM6A_expanded[1,] <- c("KDM6A", "missense", "419", NA, NA, 0, 1, "del419")

KDM6A_final <- KDM6A_expanded

#KIT
KIT <- WL_complicated %>% filter(Gene_name == "KIT")
KIT_split_slash <- string_split_slash(KIT, 4)
KIT_split_comma <- string_split_comma(KIT_split_slash, 30)

KIT_expanded <- set_columns(KIT_split_comma, NA, NA, 1, 0, NA)

KIT_expanded[1, ] <- c("KIT", "missense", "503", NA, NA, 0, 1, "insertion") 
KIT_expanded[10, ] <- c("KIT", "missense", "560", NA, NA, 0, 1, "deletion")
KIT_expanded[12, ] <- c("KIT", "missense", "579", NA, NA, 0, 1, "deletion") 
KIT_expanded[28, ] <- c("KIT", "missense", "551", NA, NA, 0, 1, "del551-559") 
KIT_expanded[29, ] <- c("KIT", "missense", "552", NA, NA, 0, 1, "del551-559") 
KIT_expanded[30, ] <- c("KIT", "missense", "553", NA, NA, 0, 1, "del551-559") 
KIT_expanded[31, ] <- c("KIT", "missense", "554", NA, NA, 0, 1, "del551-559") 
KIT_expanded[32, ] <- c("KIT", "missense", "555", NA, NA, 0, 1, "del551-559") 
KIT_expanded[33, ] <- c("KIT", "missense", "556", NA, NA, 0, 1, "del551-559") 
KIT_expanded[34, ] <- c("KIT", "missense", "557", NA, NA, 0, 1, "del551-559") 
KIT_expanded[35, ] <- c("KIT", "missense", "558", NA, NA, 0, 1, "del551-559") 
KIT_expanded[36, ] <- c("KIT", "missense", "559", NA, NA, 0, 1, "del551-559") 

KIT_final <- KIT_expanded

#MPL
#clarify W518KT (JAK, KIT, MPL)
MPL <- WL_complicated %>% filter(Gene_name == "MPL")
MPL_split_slash <- string_split_slash(MPL, 4)
MPL_split_comma <- string_split_comma(MPL_split_slash, 30)

MPL_expanded <- set_columns(MPL_split_comma, NA, NA, 1, 0, NA)



MPL_expanded[5, ] <- c("MPL", "missense", "513", NA, NA, 0, 1, "deletion")
MPL_expanded[14, ] <- c("MPL", "missense", "515", NA, NA, 0, 1, "W515-518KT")
MPL_expanded[15, ] <- c("MPL", "missense", "516", NA, NA, 0, 1, "W515-518KT")
MPL_expanded[16, ] <- c("MPL", "missense", "517", NA, NA, 0, 1, "W515-518KT")
MPL_expanded[17, ] <- c("MPL", "missense", "518", NA, NA, 0, 1, "W515-518KT")

MPL_final <- MPL_expanded

#PPM1D
PPM1D <- WL_complicated %>% filter(Gene_name == "PPM1D")
PPM1D_long <- data.frame(Gene_name ="PPM1D", Variant_types = c("Frameshift", "nonsense"), missense_type= NA,  start = NA, end = NA, whitelist=0, manual_review=1, notes="exon 5 or 6")
PPM1D_final <- PPM1D_long

#PRPF8
PRPF8 <- WL_complicated %>% filter(Gene_name == "PRPF8")
PRPF8_long <- data.frame(Gene_name ="PRPF8", Variant_types = "missense", missense_type= c("M1307I", "C1594W", "D1598Y", "D1598N", "D1598V"),  start = NA, end = NA, whitelist=1, manual_review=0, notes = NA)
PRPF8_final <- PRPF8_long

#TP53
TP53 <- WL_complicated %>% filter(Gene_name == "TP53")
TP53_split_slash <- string_split_slash(TP53, 4)
TP53_split_comma <- string_split_comma(TP53_split_slash, 600)
TP53_expanded <- set_columns(TP53_split_comma, NA, NA, 1, 0, NA)
TP53_expanded[22, ] <- c("TP53", "missense", "F134V", NA, NA, 1, 0, NA)
TP53_expanded[291, ] <- c("TP53", "missense", "F134L", NA, NA, 1, 0, NA)

TP53_final <- TP53_expanded

complex_list <- list(filtered_final, ASXL_final, CBL_CBLB_final, CREBBP_final, EP300_final, FLT3_final, GATA3_final, GNAS_final, JAK2_final, KDM6A_final, KIT_final, MPL_final, NPM1_final, PPM1D_final, PRPF8_final, TET2_final, TP53_final)

complex_char <- map(complex_list, to_char)
complete<- bind_rows(complex_char)

complete_transformed <- transform(complete, start = as.numeric(start), 
                                  end = as.numeric(end))

first_protein <- do.call('rbind', strsplit(complete_transformed$missense_type, split = "(?<=[a-zA-Z])\\s*(?=[0-9])", perl = TRUE)) %>% as_tibble()

first_protein$number <- as.integer(!grepl("^[0-9]+$", first_protein$V1)) 
first_protein$index <- rownames(complete_transformed)
first_number <- filter(first_protein, number == 0)
first_number$Initial_Protein <- NA
first_number$Protein_Position <- first_number$V2

first_letter <- filter(first_protein, number == 1)
first_letter$Initial_Protein <- first_letter$V1
first_letter$Protein_Position <- first_letter$V2

first_protein <- rbind(first_letter, first_number) %>% select(Initial_Protein, Protein_Position, index)

last_protein<-  do.call('rbind', strsplit(first_protein$Protein_Position, split = "(?<=[0-9])\\s*(?=[a-zA-Z])", perl = TRUE)) %>% as_tibble()
last_protein$index <- first_protein$index
last_protein$V2 <- stringr::str_replace(last_protein$V2, '\\*', '')
last_protein$number <- as.integer(!grepl("^[0-9]+$", last_protein$V2))

last_number <- filter(last_protein, number == 0)
last_number$Protein_Position <- last_number$V2
last_number$Final_Protein <- "*"

last_letter <- filter(last_protein, number == 1)
last_letter$Protein_Position <- last_letter$V1
last_letter$Final_Protein <- last_letter$V2

last_protein_position <- rbind(last_letter, last_number) %>% select(Protein_Position, Final_Protein, index)
first_protein <- first_protein %>% select(-Protein_Position)

Protein_Positions <- merge(first_protein, last_protein_position, by="index") 

complete_transformed$index <- rownames(complete_transformed)
final <- merge(complete_transformed, Protein_Positions, by="index") %>% transform(Protein_Position = as.numeric(Protein_Position)) %>% select(-index)

final_complex_missense <- final %>% filter(Variant_types == "missense") %>% filter(!is.na(start))
CBL <- do.call("rbind", replicate(41, final_complex_missense[1,], simplify = FALSE)) %>% mutate(missense_type = c(381:421))
CBLB <- do.call("rbind", replicate(41, final_complex_missense[2,], simplify = FALSE)) %>% mutate(missense_type = c(372:412))
TET2_Missense <- do.call("rbind", replicate(899, final_complex_missense[3,], simplify = FALSE)) %>% mutate(missense_type = c(1104:2002))

final_complex_nonmissense <- final %>% filter(Variant_types != "missense") %>% filter(!is.na(start))
ASXL1 <- do.call("rbind", replicate(1184, final_complex_nonmissense[1,], simplify = FALSE)) %>% mutate(Protein_Position = c(358:1541))
ASXL2 <- do.call("rbind", replicate(1162, final_complex_nonmissense[4,], simplify = FALSE)) %>% mutate(Protein_Position = c(380:1541))
NPM1 <- do.call("rbind", replicate(6, final_complex_nonmissense[7,], simplify = FALSE)) %>% mutate(Protein_Position = c(859:864))
TET2_Frameshift <- do.call("rbind", replicate(899, final_complex_nonmissense[11,], simplify = FALSE)) %>% mutate(Protein_Position = c(1104:2002))
TET2_nonsense <- do.call("rbind", replicate(899, final_complex_nonmissense[13,], simplify = FALSE)) %>% mutate(Protein_Position = c(1104:2002))
TET2_splicesite <- do.call("rbind", replicate(899, final_complex_nonmissense[13,], simplify = FALSE)) %>% mutate(Protein_Position = c(1104:2002))

final_complex_merged <- rbind(CBL, CBLB, TET2_Missense, TET2_Frameshift, TET2_nonsense, TET2_splicesite, ASXL1, ASXL2, NPM1)

final$index <- rownames(final)
final <- final[!(final$index %in% c(904:911, 1038:1041, 1266:1272, 1073, 1090)),]
final_merged <- final %>% select(-index) %>% rbind(final_complex_merged)

whitelist_file_location <- str_c(output, "/variant_whitelist.xlsx")

if (file.exists(whitelist_file_location) == TRUE) {
    file.remove(whitelist_file_location)
}

write.xlsx(final_merged, whitelist_file_location, keepNA = TRUE)
