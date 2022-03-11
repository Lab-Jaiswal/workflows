library(readxl)
library(magrittr)
library(tidyverse)
require(data.table)

rename <- function(df, column, new){
  x <- names(df)                               #Did this to avoid typing twice
  if (is.numeric(column)) column <- x[column]  #Take numeric input by indexing
  names(df)[x %in% column] <- new              #What you're interested in
  return(df)
}


command_args <- commandArgs(trailingOnly = TRUE)
whitelist_coordinates <- "/home/maurertm/labs/variant_whitelist_real.xls"
mutect_coordinates <- "/home/maurertm/smontgom/maurertm/mutect_aggregated_noWL_Mar1_1147.tsv"
output <- "/home/maurertm/smontgom/maurertm/"

mutect<-as.data.frame(fread(mutect_coordinates))
white_list <- read_excel(whitelist_coordinates, col_names=T)

Missense_pattern <- c("MISS", "Miss", "miss")
Frameshift_pattern <- c("FRAME", "Frame", "frame")
Splice_site_pattern <- c("SPLICE", "Splice", "splice")
Nonsense_pattern <- c("NONSENSE", "Nonsense", "NonSense", "nonsense")

Missense <- (grep(paste(Missense_pattern,collapse="|"), 
                  mutect$Variant_Classification))
Frame_shift <- (grep(paste(Frameshift_pattern,collapse="|"), 
                     mutect$Variant_Classification))
Splice_site <- (grep(paste(Splice_site_pattern,collapse="|"), 
                     mutect$Variant_Classification))
Nonsense <- (grep(paste(Nonsense_pattern,collapse="|"), 
                  mutect$Variant_Classification))
for (i in Missense){
  mutect[i, "Variant_Classification"] <- 'missense'
}
for (i in Frame_shift){
  mutect[i, "Variant_Classification"]  <- 'Frameshift'
}
for (i in Splice_site){
  mutect[i, "Variant_Classification"] <- 'splice-site'
}
for (i in Nonsense){
  mutect[i, "Variant_Classification"] <- 'nonsense'
}

#Get Initial_Protein, Protein_Position, and Final_Protein column from mutect
#First find the values in Protein Change that contain multiple positions
mutect_no_p <- mutate(mutect, Protein_Change_No_P = gsub("^.*\\.","", mutect$Protein_Change)) 
normal <- dplyr::filter(mutect_no_p, !grepl("_", Protein_Change_No_P)) %>% dplyr::filter(!grepl("del", Protein_Change_No_P)) %>% dplyr::filter(!grepl("\\*", Protein_Change_No_P)) %>% dplyr::filter(Protein_Change_No_P != "")

deviant_insertion <- dplyr::filter(mutect_no_p, grepl("_", Protein_Change_No_P))
deviant_del <- dplyr::filter(mutect_no_p, grepl("del", Protein_Change_No_P))
deviant_asterisk <- dplyr::filter(mutect_no_p, grepl("\\*", Protein_Change_No_P))
deviant_empty <-  dplyr::filter(mutect_no_p, Protein_Change_No_P == "")

#split up deviant by insertion type
deviant_ins <- dplyr::filter(deviant_insertion, grepl("ins", Protein_Change_No_P))
deviant_no_ins <- dplyr::filter(deviant_insertion, !grepl("ins", Protein_Change_No_P))

#clean deviant_ins
    deviant_ins_split <- str_split_fixed(deviant_ins$Protein_Change_No_P, '_', 2) %>% as_tibble()
    deviant_ins_dataframe <- mutate(deviant_ins_split, Intermediate = gsub("[[:digit:]]","", V2)) %>% 
        mutate(Final_Protein = gsub("ins","", Intermediate)) %>% 
        mutate(Protein_Position = gsub("[^0-9.-]", "", V2)) %>% 
        mutate(Initial_Protein = NA) %>%
        select(Initial_Protein, Protein_Position, Final_Protein)
    deviant_ins <- cbind(deviant_ins, deviant_ins_dataframe)

#clean deviant_no_ins
     deviant_noins_split <- str_split_fixed(deviant_no_ins$Protein_Change_No_P, '_', 2) %>% as_tibble()
     deviant_noins_dataframe <- mutate(deviant_noins_split, Intermediate = gsub("[[:digit:]]","", V2)) %>% 
        mutate(Final_Protein = gsub("[^>]*>(.*)", "\\1", Intermediate)) %>% 
        mutate(Protein_Position = gsub("[^0-9.-]", "", V2)) %>% 
        mutate(Initial_Protein = gsub("(.*)>.*", "\\1", Intermediate) ) %>%
        select(Initial_Protein, Protein_Position, Final_Protein)
    deviant_no_ins <- cbind(deviant_no_ins, deviant_noins_dataframe)

#clean deviant_del    
     deviant_del_dataframe <- deviant_del %>% 
        mutate(Final_Protein = NA) %>% 
        mutate(Protein_Position = gsub("[^0-9.-]", "", Protein_Change_No_P)) %>% 
        mutate(Intermediate = gsub("\\d+", "", Protein_Change_No_P)) %>%
        mutate(Initial_Protein = gsub("del", "", Intermediate)) %>%
        select(Initial_Protein, Protein_Position, Final_Protein)
    deviant_del <- cbind(deviant_del, deviant_del_dataframe)

#clean asterisk 
    deviant_ast_dataframe <- deviant_asterisk %>% 
        mutate(Final_Protein = NA) %>% 
        mutate(Protein_Position = gsub("[^0-9.-]", "", Protein_Change_No_P)) %>% 
        mutate(Intermediate = gsub("\\d+", "", Protein_Change_No_P)) %>%
        mutate(Initial_Protein = gsub("\\*", "", Intermediate)) %>% 
        select(Initial_Protein, Protein_Position, Final_Protein)
    deviant_asterisk <- cbind(deviant_asterisk, deviant_ast_dataframe)

#clean deviant empty
    deviant_empty_dataframe <- deviant_empty %>% 
        mutate(Final_Protein = NA) %>% 
        mutate(Protein_Position = NA) %>% 
        mutate(Initial_Protein = NA) %>% 
        select(Initial_Protein, Protein_Position, Final_Protein)
    deviant_empty <- cbind(deviant_empty, deviant_empty_dataframe)

#clean normal
first_protein <- do.call('rbind', strsplit(normal$Protein_Change_No_P, split = "(?<=[a-zA-Z])\\s*(?=[0-9])", perl = TRUE)) %>% as_tibble()
first_protein$number <- as.integer(!grepl("^[0-9]+$", first_protein$V1)) 
first_protein$index <- rownames(normal)
first_number <- filter(first_protein, number == 0)
first_number$Initial_Protein <- NA
first_number$Protein_Position <- first_number$V2

first_letter <- filter(first_protein, number == 1)
first_letter$Initial_Protein <- first_letter$V1
first_letter$Protein_Position <- first_letter$V2

first_protein <- rbind(first_letter, first_number) %>% dplyr::select(Initial_Protein, Protein_Position, index)

last_protein <-  do.call('rbind', strsplit(first_protein$Protein_Position, split = "(?<=[0-9])\\s*(?=[a-zA-Z])", perl = TRUE)) %>% as_tibble()
last_protein$index <- first_protein$index
last_protein$V2 <- stringr::str_replace(last_protein$V2, '\\*', '')
last_protein$number <- as.integer(!grepl("^[0-9]+$", last_protein$V2))

last_number <- filter(last_protein, number == 0)
last_number$Protein_Position <- last_number$V2
last_number$Final_Protein <- NA

last_letter <- filter(last_protein, number == 1)
last_letter$Protein_Position <- last_letter$V1
last_letter$Final_Protein <- last_letter$V2

last_protein_position <- rbind(last_letter, last_number) %>% dplyr::select(Protein_Position, Final_Protein, index)
first_protein <- first_protein %>% dplyr::select(-Protein_Position)

Protein_Positions <- merge(first_protein, last_protein_position, by="index")

normal$index <- rownames(normal)
normal <- merge(normal, Protein_Positions, by="index") %>% select(-index)
   
mutect_joined <- rbind(normal, deviant_del, deviant_ins, deviant_no_ins, deviant_asterisk, deviant_empty)
mutect_annotated_joined <- transform(mutect_joined, Protein_Positions = as.numeric(Protein_Position))
mutect_annotated_joined$index <- rownames(mutect_annotated_joined)

####################NON-MISSENSE##############################
##############################################################
#Simple cases, where the p.Position does not affect if the Frame Shift, Splice Site, or Non-Sense mutation is whitelisted 
complex_nonmissense <- c("ASXL1", "ASXL2", "TET2", "NPM1")
non_missense_mutect_simple <- mutect_annotated_joined %>% filter(Variant_Classification != "missense") %>% filter(!is_in(Hugo_Symbol, complex_nonmissense)) 
non_missense_whitelist_simple <- white_list %>% filter(Variant_types!= "missense") %>% filter(!is_in(Gene_name, complex_nonmissense)) %>% dplyr::select(-missense_type, -Initial_Protein, -Protein_Position, -Final_Protein)

non_missense_mutect_simple$Hugo_Variant <- paste(non_missense_mutect_simple$Hugo_Symbol, non_missense_mutect_simple$Variant_Classification, sep = "_")
non_missense_whitelist_simple$Hugo_Variant <- paste(non_missense_whitelist_simple$Gene_name, non_missense_whitelist_simple$Variant_types, sep = "_")

WL_NonMissense_Simple_Gene_Variants <- non_missense_whitelist_simple$Hugo_Variant
#mutect <- unique(non_missense_mutect_simple$Hugo_Variant)
#wl <- unique(non_missense_whitelist_simple$Hugo_Variant)

#intersect(wl, mutect)
mutect_nonmissense_variants_simple <- non_missense_mutect_simple %>% filter(is_in(Hugo_Variant, WL_NonMissense_Simple_Gene_Variants)) %>% left_join(non_missense_whitelist_simple, by="Hugo_Variant") %>% dplyr::select(-Hugo_Variant, -Protein_Change_No_P, -Protein_Positions)

#Complex cases, where the p.Position affects if the Frame Shift, Splice Site, or Non-Sense mutation is whitelisted
non_missense_mutect_complex <- mutect_annotated_joined %>% filter(Variant_Classification != "missense") %>% filter(is_in(Hugo_Symbol, complex_nonmissense)) 
non_missense_whitelist_complex <- white_list %>% filter(Variant_types!= "missense") %>% filter(is_in(Gene_name, complex_nonmissense)) %>% dplyr::select(-missense_type, -Initial_Protein, -Final_Protein)

non_missense_mutect_complex$Hugo_Variant_Position <- paste(non_missense_mutect_complex$Hugo_Symbol, non_missense_mutect_complex$Variant_Classification, non_missense_mutect_complex$Protein_Position, sep = "_")
non_missense_mutect_complex <- non_missense_mutect_complex %>% dplyr::select(-Protein_Position)
non_missense_whitelist_complex$Hugo_Variant_Position <- paste(non_missense_whitelist_complex$Gene_name, non_missense_whitelist_complex$Variant_types, non_missense_whitelist_complex$Protein_Position, sep = "_")

WL_NonMissense_Complex_Gene_Variants <- non_missense_whitelist_complex$Hugo_Variant_Position
#non_missense_mutect_complex$Hugo_Variant_Position

mutect_nonmissense_variants_complex <- non_missense_mutect_complex %>% filter(is_in(Hugo_Variant_Position, WL_NonMissense_Complex_Gene_Variants)) %>% left_join(non_missense_whitelist_complex, by="Hugo_Variant_Position") %>% dplyr::select(-Protein_Change_No_P, -Protein_Positions, -Hugo_Variant_Position)

####################MISSENSE####################################
#################################################################
#Simple Cases, where there is both a Initial_Protein and Final_Protein value
missense_mutect <- mutect_annotated_joined %>% filter(Variant_Classification == "missense") 
missense_whitelist_simple <- white_list %>% filter(Variant_types == "missense") %>% filter(!is.na(Initial_Protein)) %>% filter(!is.na(Final_Protein))
missense_mutect$Hugo_ProteinChange <- paste(missense_mutect$Hugo_Symbol, missense_mutect$Protein_Change_No_P, sep = "_")
missense_whitelist_simple$Hugo_ProteinChange <- paste(missense_whitelist_simple$Gene_name, missense_whitelist_simple$missense_type, sep = "_")
missense_whitelist_simple <- missense_whitelist_simple %>% dplyr::select(-missense_type)

WL_Missense_Simple_Gene_Variants <- missense_whitelist_simple$Hugo_ProteinChange

mutect_missense_variants_simple <- missense_mutect %>% dplyr::select(-Initial_Protein, -Protein_Position, -Final_Protein) %>% filter(is_in(Hugo_ProteinChange, WL_Missense_Simple_Gene_Variants)) %>% left_join(missense_whitelist_simple, by="Hugo_ProteinChange") %>% dplyr::select(-Hugo_ProteinChange, -Protein_Change_No_P, -Protein_Positions)


#Account for complex cases where their is an insertion or deletion
missense_whitelist_complex_first <- white_list %>% filter(Variant_types == "missense") %>% filter(Initial_Protein == "NA" | is.na(Initial_Protein))
missense_whitelist_complex_last <- white_list %>% filter(Variant_types == "missense") %>% filter(Final_Protein == "NA" | is.na(Final_Protein))
missense_whitelist_complex <- rbind(missense_whitelist_complex_first, missense_whitelist_complex_last) %>% dplyr::select(-Initial_Protein, -Protein_Position, -Final_Protein) 

missense_whitelist_complex$Hugo_ProteinPosition <- paste(missense_whitelist_complex$Gene_name, missense_whitelist_complex$missense_type, sep = "_")
missense_whitelist_complex <- missense_whitelist_complex %>% dplyr::select(-missense_type)

missense_mutect$Hugo_ProteinPosition <- paste(missense_mutect$Hugo_Symbol, missense_mutect$Protein_Position, sep = "_")

WL_Missense_Complex_Gene_Variants <- missense_whitelist_complex$Hugo_ProteinPosition

mutect_missense_variants_complex <- missense_mutect %>% filter(is_in(Hugo_ProteinPosition, WL_Missense_Complex_Gene_Variants)) %>% left_join(missense_whitelist_complex, by="Hugo_ProteinPosition") %>% dplyr::select(-Hugo_ProteinPosition, -Hugo_ProteinChange, -Protein_Change_No_P, -Protein_Positions)

###############JOIN MISSENSE AND NON-MISSENSE####################
#################################################################
whitelisted_variants <- rbind(mutect_nonmissense_variants_simple, mutect_nonmissense_variants_complex, mutect_missense_variants_simple, mutect_missense_variants_complex) %>% unique() %>% dplyr::select(-"Gene_name", -"Variant_types", -"start", -"end")

whitelisted_index <- as.numeric(unlist(whitelisted_variants$index))
mutect_no_WL <- mutect_annotated_joined %>% filter(!is_in(index, whitelisted_index))
mutect_no_WL$whitelist <- 0
mutect_no_WL$manual_review <- 0
mutect_no_WL$notes <- NA

col_WL <- colnames(whitelisted_variants)
not_whitelisted <- dplyr::select(mutect_no_WL, all_of(col_WL))

mutect_whitelist_annotated <- rbind(whitelisted_variants, not_whitelisted) %>% arrange(Hugo_Symbol)

annotated_file_location <- str_c(output, "/mutect_annotated_whitelist_Mar1_431.tsv")

if (file.exists(annotated_file_location) == TRUE) {
  file.remove(annotated_file_location)
}

write_tsv(mutect_whitelist_annotated, annotated_file_location)
