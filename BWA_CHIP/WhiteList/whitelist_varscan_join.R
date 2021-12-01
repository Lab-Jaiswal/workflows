library(readxl)
library(magrittr)
library(tidyverse)

rename <- function(df, column, new){
    x <- names(df)                               #Did this to avoid typing twice
    if (is.numeric(column)) column <- x[column]  #Take numeric input by indexing
    names(df)[x %in% column] <- new              #What you're interested in
    return(df)
}


command_args <- commandArgs(trailingOnly = TRUE)
whitelist_coordinates <- command_args[1]
varscan_coordinates <- command_args[2]
output <- command_args[3]

white_list <- read_excel(whitelist_coordinates, col_names=T)
varscan <- read_table(varscan_coordinates)

varscan$Variant_Classification <- varscan$ExonicFunc.ensGene
varscan$Variant_Classification <- gsub('nonsynonymous_SNV', 'missense', varscan$Variant_Classification)
varscan$Variant_Classification <- gsub("frameshift_deletion", "Frameshift", varscan$Variant_Classification)
white_list$Variant_types <- gsub("splice-site", "Splicesite", white_list$Variant_types)

#Get Initial_Protein, Protein_Position, and Final_Protein column from varscan
no_p <- do.call('rbind', strsplit(varscan$Protein_change, split = "p.", perl = TRUE)) %>% as_tibble()
varscan <- bind_cols(varscan, no_p)
first_protein <- do.call('rbind', strsplit(no_p$V2, split = "(?<=[a-zA-Z])\\s*(?=[0-9])", perl = TRUE)) %>% as_tibble()

first_protein$number <- as.integer(!grepl("^[0-9]+$", first_protein$V1)) 
first_protein$index <- rownames(varscan)
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
last_number$Final_Protein <- NA

last_letter <- filter(last_protein, number == 1)
last_letter$Protein_Position <- last_letter$V1
last_letter$Final_Protein <- last_letter$V2

last_protein_position <- rbind(last_letter, last_number) %>% select(Protein_Position, Final_Protein, index)
first_protein <- first_protein %>% select(-Protein_Position)

Protein_Positions <- merge(first_protein, last_protein_position, by="index") 

varscan$index <- rownames(varscan)
varscan <- merge(varscan, Protein_Positions, by="index")
varscan <- transform(varscan, Protein_Position = as.numeric(Protein_Position))
varscan <- rename(varscan, "Gene.ensGene", "Hugo_Symbol")

####################NON-MISSENSE##############################
##############################################################
#Simple cases, where the p.Position does not affect if the Frame Shift, Splice Site, or Non-Sense mutation is whitelisted 
complex_nonmissense <- c("ASXL1", "ASXL2", "TET2", "NPM1")
non_missense_varscan_simple <- varscan %>% filter(Variant_Classification != "missense") %>% filter(!is_in(Hugo_Symbol, complex_nonmissense))
non_missense_whitelist_simple <- white_list %>% filter(Variant_types!= "missense") %>% filter(!is_in(Gene_name, complex_nonmissense)) %>% select(-missense_type, -Initial_Protein, -Protein_Position, -Final_Protein)

non_missense_varscan_simple$Hugo_Variant <- paste(non_missense_varscan_simple$Hugo_Symbol, non_missense_varscan_simple$Variant_Classification, sep = "_")
non_missense_whitelist_simple$Hugo_Variant <- paste(non_missense_whitelist_simple$Gene_name, non_missense_whitelist_simple$Variant_types, sep = "_")

WL_NonMissense_Simple_Gene_Variants <- non_missense_whitelist_simple$Hugo_Variant


varscan_nonmissense_variants_simple <- non_missense_varscan_simple %>% filter(is_in(Hugo_Variant, WL_NonMissense_Simple_Gene_Variants)) %>% left_join(non_missense_whitelist_simple, by="Hugo_Variant") %>% select(-Hugo_Variant)

#Complex cases, where the p.Position affects if the Frame Shift, Splice Site, or Non-Sense mutation is whitelisted
non_missense_varscan_complex <- varscan %>% filter(Variant_Classification != "missense") %>% filter(is_in(Hugo_Symbol, complex_nonmissense))
non_missense_whitelist_complex <- white_list %>% filter(Variant_types!= "missense") %>% filter(is_in(Gene_name, complex_nonmissense)) %>% select(-missense_type, -Initial_Protein, -Final_Protein)

non_missense_varscan_complex$Hugo_Variant_Position <- paste(non_missense_varscan_complex$Hugo_Symbol, non_missense_varscan_complex$Variant_Classification, non_missense_varscan_complex$Protein_Position, sep = "_")
non_missense_varscan_complex <- non_missense_varscan_complex %>% select(-Protein_Position)
non_missense_whitelist_complex$Hugo_Variant_Position <- paste(non_missense_whitelist_complex$Gene_name, non_missense_whitelist_complex$Variant_types, non_missense_whitelist_complex$Protein_Position, sep = "_")

WL_NonMissense_Complex_Gene_Variants <- non_missense_whitelist_complex$Hugo_Variant_Position


varscan_nonmissense_variants_complex <- non_missense_varscan_complex %>% filter(is_in(Hugo_Variant_Position, WL_NonMissense_Complex_Gene_Variants)) %>% left_join(non_missense_whitelist_complex, by="Hugo_Variant_Position") %>% select(-Hugo_Variant_Position)

####################MISSENSE####################################
#################################################################
#Simple Cases, where there is both a Initial_Protein and Final_Protein value
missense_varscan <- varscan %>% filter(Variant_Classification == "missense")
missense_whitelist_simple <- white_list %>% filter(Variant_types == "missense") %>% filter(!is.na(Initial_Protein)) %>% filter(!is.na(Final_Protein)) %>% select(-Initial_Protein, -Protein_Position, -Final_Protein)
missense_varscan$Hugo_ProteinChange <- paste(missense_varscan$Hugo_Symbol, missense_varscan$V2, sep = "_")
missense_whitelist_simple$Hugo_ProteinChange <- paste(missense_whitelist_simple$Gene_name, missense_whitelist_simple$missense_type, sep = "_")
missense_whitelist_simple <- missense_whitelist_simple %>% select(-missense_type)

WL_Missense_Simple_Gene_Variants <- missense_whitelist_simple$Hugo_ProteinChange

varscan_missense_variants_simple <- missense_varscan %>% filter(is_in(Hugo_ProteinChange, WL_Missense_Simple_Gene_Variants)) %>% left_join(missense_whitelist_simple, by="Hugo_ProteinChange") %>% select(-Hugo_ProteinChange)

#Account for complex cases where their is an insertion or deletion
missense_whitelist_complex_first <- white_list %>% filter(Variant_types == "missense") %>% filter(is.na(Initial_Protein)) 
missense_whitelist_complex_last <- white_list %>% filter(Variant_types == "missense") %>% filter(is.na(Final_Protein)) 
missense_whitelist_complex <- rbind(missense_whitelist_complex_first, missense_whitelist_complex_last) %>% select(-Initial_Protein, -Protein_Position, -Final_Protein)

missense_whitelist_complex$Hugo_ProteinPosition <- paste(missense_whitelist_complex$Gene_name, missense_whitelist_complex$missense_type, sep = "_")
missense_whitelist_complex <- missense_whitelist_complex %>% select(-missense_type)

missense_varscan$Hugo_ProteinPosition <- paste(missense_varscan$Hugo_Symbol, missense_varscan$Protein_Position, sep = "_")

WL_Missense_Complex_Gene_Variants <- missense_whitelist_complex$Hugo_ProteinPosition

varscan_missense_variants_complex <- missense_varscan %>% filter(is_in(Hugo_ProteinPosition, WL_Missense_Complex_Gene_Variants)) %>% left_join(missense_whitelist_complex, by="Hugo_ProteinPosition") %>% select(-Hugo_ProteinPosition, Hugo_ProteinChange)

###############JOIN MISSENSE AND NON-MISSENSE####################
#################################################################
whitelisted_variants <- rbind(varscan_nonmissense_variants_simple, varscan_nonmissense_variants_complex, varscan_missense_variants_simple, varscan_nonmissense_variants_complex) %>% unique() %>% select(-"Gene_name", -"Variant_types", -"start", -"end")
whitelisted_index <- as.numeric(unlist(whitelisted_variants$index))
varscan <- varscan %>% filter(!is_in(index, whitelisted_index))
varscan$whitelist <- 0
varscan$manual_review <- 0
varscan$notes <- NA

varscan_whitelist_annotated <- rbind(whitelisted_variants, varscan) %>% arrange(Hugo_Symbol)

annotated_file_location <- str_c(output, "/varscan_annotated.tsv")

if (file.exists(annotated_file_location) == TRUE) {
    file.remove(annotated_file_location)
}

write_tsv(varscan_whitelist_annotated, annotated_file_location)
