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
mutect_coordinates <- command_args[2]
output <- command_args[3]

white_list <- read_excel(whitelist_coordinates, col_names=T)
mutect <- read_excel(mutect_coordinates, col_names=T)

mutect$Variant_Classification <- gsub('MISSENSE', 'missense', mutect$Variant_Classification)
mutect$Variant_Classification <- gsub("FRAME_SHIFT_DEL", "Frameshift", mutect$Variant_Classification)
mutect$Variant_Classification <- gsub("SPLICE_SITE", "Splicesite", mutect$Variant_Classification)
mutect$Variant_Classification <- gsub("NONSENSE", "nonsense", mutect$Variant_Classification)
white_list$Variant_types <- gsub("splice-site", "Splicesite", white_list$Variant_types)

manual_review <- white_list %>% filter(manual_review == 1)
white_listed <- white_list %>% filter(whitelist == 1)
#Get Initial_Protein, Protein_Position, and Final_Protein column from mutect
no_p <- do.call('rbind', strsplit(mutect$Protein_Change, split = "p.", perl = TRUE)) %>% as_tibble()
mutect <- bind_cols(mutect, no_p)
first_protein <- do.call('rbind', strsplit(no_p$V2, split = "(?<=[a-zA-Z])\\s*(?=[0-9])", perl = TRUE)) %>% as_tibble()

first_protein$number <- as.integer(!grepl("^[0-9]+$", first_protein$V1)) 
first_protein$index <- rownames(mutect)
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

mutect$index <- rownames(mutect)
mutect <- merge(mutect, Protein_Positions, by="index")
mutect <- transform(mutect, Protein_Position = as.numeric(Protein_Position))

####################NON-MISSENSE##############################
##############################################################
#Simple cases, where the p.Position does not affect if the Frame Shift, Splice Site, or Non-Sense mutation is whitelisted 
complex_nonmissense <- c("ASXL1", "ASXL2", "TET2", "NPM1")
non_missense_mutect_simple <- mutect %>% filter(Variant_Classification != "missense") %>% filter(!is_in(Hugo_Symbol, complex_nonmissense))
non_missense_whitelist_simple <- white_list %>% filter(Variant_types!= "missense") %>% filter(!is_in(Gene_name, complex_nonmissense)) %>% select(-missense_type, -Initial_Protein, -Protein_Position, -Final_Protein)

non_missense_mutect_simple$Hugo_Variant <- paste(non_missense_mutect_simple$Hugo_Symbol, non_missense_mutect_simple$Variant_Classification, sep = "_")
non_missense_whitelist_simple$Hugo_Variant <- paste(non_missense_whitelist_simple$Gene_name, non_missense_whitelist_simple$Variant_types, sep = "_")

WL_NonMissense_Simple_Gene_Variants <- non_missense_whitelist_simple$Hugo_Variant


mutect_nonmissense_variants_simple <- non_missense_mutect_simple %>% filter(is_in(Hugo_Variant, WL_NonMissense_Simple_Gene_Variants)) %>% left_join(non_missense_whitelist_simple, by="Hugo_Variant") %>% select(-Hugo_Variant)

#Complex cases, where the p.Position affects if the Frame Shift, Splice Site, or Non-Sense mutation is whitelisted
non_missense_mutect_complex <- mutect %>% filter(Variant_Classification != "missense") %>% filter(is_in(Hugo_Symbol, complex_nonmissense))
non_missense_whitelist_complex <- white_list %>% filter(Variant_types!= "missense") %>% filter(is_in(Gene_name, complex_nonmissense)) %>% select(-missense_type, -Initial_Protein, -Final_Protein)

non_missense_mutect_complex$Hugo_Variant_Position <- paste(non_missense_mutect_complex$Hugo_Symbol, non_missense_mutect_complex$Variant_Classification, non_missense_mutect_complex$Protein_Position, sep = "_")
non_missense_mutect_complex <- non_missense_mutect_complex %>% select(-Protein_Position)
non_missense_whitelist_complex$Hugo_Variant_Position <- paste(non_missense_whitelist_complex$Gene_name, non_missense_whitelist_complex$Variant_types, non_missense_whitelist_complex$Protein_Position, sep = "_")

WL_NonMissense_Complex_Gene_Variants <- non_missense_whitelist_complex$Hugo_Variant_Position


mutect_nonmissense_variants_complex <- non_missense_mutect_complex %>% filter(is_in(Hugo_Variant_Position, WL_NonMissense_Complex_Gene_Variants)) %>% left_join(non_missense_whitelist_complex, by="Hugo_Variant_Position") %>% select(-Hugo_Variant_Position)

####################MISSENSE####################################
#################################################################
#Simple Cases, where there is both a Initial_Protein and Final_Protein value
missense_mutect <- mutect %>% filter(Variant_Classification == "missense")
missense_whitelist_simple <- white_list %>% filter(Variant_types == "missense") %>% filter(!is.na(Initial_Protein)) %>% filter(!is.na(Final_Protein)) %>% select(-Initial_Protein, -Protein_Position, -Final_Protein)
missense_mutect$Hugo_ProteinChange <- paste(missense_mutect$Hugo_Symbol, missense_mutect$V2, sep = "_")
missense_whitelist_simple$Hugo_ProteinChange <- paste(missense_whitelist_simple$Gene_name, missense_whitelist_simple$missense_type, sep = "_")
missense_whitelist_simple <- missense_whitelist_simple %>% select(-missense_type)

WL_Missense_Simple_Gene_Variants <- missense_whitelist_simple$Hugo_ProteinChange

mutect_missense_variants_simple <- missense_mutect %>% filter(is_in(Hugo_ProteinChange, WL_Missense_Simple_Gene_Variants)) %>% left_join(missense_whitelist_simple, by="Hugo_ProteinChange") %>% select(-Hugo_ProteinChange)

#Account for complex cases where their is an insertion or deletion
missense_whitelist_complex_first <- white_list %>% filter(Variant_types == "missense") %>% filter(is.na(Initial_Protein)) 
missense_whitelist_complex_last <- white_list %>% filter(Variant_types == "missense") %>% filter(is.na(Final_Protein)) 
missense_whitelist_complex <- rbind(missense_whitelist_complex_first, missense_whitelist_complex_last) %>% select(-Initial_Protein, -Protein_Position, -Final_Protein)

missense_whitelist_complex$Hugo_ProteinPosition <- paste(missense_whitelist_complex$Gene_name, missense_whitelist_complex$missense_type, sep = "_")
missense_whitelist_complex <- missense_whitelist_complex %>% select(-missense_type)

missense_mutect$Hugo_ProteinPosition <- paste(missense_mutect$Hugo_Symbol, missense_mutect$Protein_Position, sep = "_")

WL_Missense_Complex_Gene_Variants <- missense_whitelist_complex$Hugo_ProteinPosition

mutect_missense_variants_complex <- missense_mutect %>% filter(is_in(Hugo_ProteinPosition, WL_Missense_Complex_Gene_Variants)) %>% left_join(missense_whitelist_complex, by="Hugo_ProteinPosition") %>% select(-Hugo_ProteinPosition, Hugo_ProteinChange)

###############JOIN MISSENSE AND NON-MISSENSE####################
#################################################################
whitelisted_variants <- rbind(mutect_nonmissense_variants_simple, mutect_nonmissense_variants_complex, mutect_missense_variants_simple, mutect_nonmissense_variants_complex) %>% unique() %>% select(-"Gene_name", -"Variant_types", -"start", -"end")
whitelisted_index <- as.numeric(unlist(whitelisted_variants$index))
mutect <- mutect %>% filter(!is_in(index, whitelisted_index))
mutect$whitelist <- 0
mutect$manual_review <- 0
mutect$notes <- NA

mutect_whitelist_annotated <- rbind(whitelisted_variants, mutect) %>% arrange(Hugo_Symbol)

annotated_file_location <- str_c(output, "/mutect_annotated.tsv")

if (file.exists(annotated_file_location) == TRUE) {
    file.remove(annotated_file_location)
}

write_tsv(mutect_whitelist_annotated, annotated_file_location)
