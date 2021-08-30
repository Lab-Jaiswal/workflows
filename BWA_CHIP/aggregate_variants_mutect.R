# Call the necessary libraries
library(maftools)
library(GenomicRanges)
library(readxl)
library(magrittr)
library(tidyverse)

# Function input: col_names and col_values as arguments
# Function goal: assigns the names of the col_values to be the col_names
names_set <- function(col_names, col_values) {
  names(col_values) <- col_names
  col_values
}

# Function input: .maf file
# Function goal: Extract the maf data 
make_data = function(maf){
  maf@data
}

# Function input: genomic content in the ref_content column
# Function goal: Get the longest repeat for each variant/ row
# Steps: (1) To get repeat_a: Replace all C, G, and T's with ";", split on every ";", extract and count the elements (a's) w/in the list
    # (2) Perform the above to get repeat_c, repeat_g, and repeat_t. 
    # (3) Find which repeat (a, t, g, or c) is the longest for every given variant.
detect_repeat <- function(ref_context) {
  repeat_a <- str_replace_all(ref_context, "C|G|T", ";") %>% str_split(";") %>% extract2(1) %>% nchar %>% max
  repeat_c <- str_replace_all(ref_context, "A|G|T", ";") %>% str_split(";") %>% extract2(1) %>% nchar %>% max
  repeat_g <- str_replace_all(ref_context, "A|C|T", ";") %>% str_split(";") %>% extract2(1) %>% nchar %>% max
  repeat_t <- str_replace_all(ref_context, "A|C|G", ";") %>% str_split(";") %>% extract2(1) %>% nchar %>% max
  max(repeat_a, repeat_c, repeat_g, repeat_t)
}

# Function input: one column containing values split up by ","'s, names of the new columns
# Function goals: Split the column on each "," to create as many new columns as mutect_vcf_ncol is equal to. Then, name the new columns.
split_columns <- function(column, nocolumns, name1, name2) {
  object <- str_replace_all(column, "\\[|\\]", "") %>%
    str_split_fixed(",", nocolumns)  %>%
    as_tibble() %>%
    mutate(across(everything(), as.numeric))  
  colnames(object) <- c(name1, str_c(name2, seq(1:(nocolumns- 1)))) 
  object
}

# Define the location of the FastQs and the twist panel
# Use the format of the end of the files to select the maf, vcf files, and sample names
command_args <- commandArgs(trailingOnly = TRUE)
panel_coordinates <- command_args[1]
mutect_directory <- command_args[2]

sample_names <- list.files(mutect_directory, pattern = "*_funcotator.vcf$") %>% str_remove_all("_.*$")

maf_files <- list.files(mutect_directory, pattern = "*_funcotator.maf$", full.names = TRUE)
vcf_files <- list.files(mutect_directory, pattern = "*_funcotator.vcf$", full.names = TRUE)

# Extract the vcf files, as well as the information in the INFO, FORMAT, and DATA columns
# Bind all of the newly formatted vcf files together horizontally
mutect_vcf_header <- read_lines(vcf_files[1]) 

vcf_colnames <- str_subset(mutect_vcf_header, "^#") %>% 
        str_subset("^##", negate = T) %>% 
        str_split("\\t") %>% 
        extract2(1)
vcf_colnames[length(vcf_colnames)] <- "DATA"
mutect_vcf_list <- map(vcf_files, read_lines) %>% 
        map(str_subset, "^# ", negate = TRUE) %>% 
        map(str_split_fixed, "\\t", length(vcf_colnames)) %>% 
        map(set_colnames, vcf_colnames) %>% 
        map(as_tibble) %>% 
        set_names(sample_names)
mutect_vcf_all <- bind_rows(mutect_vcf_list, .id = "Sample")

mutect_info_names <- str_split(mutect_vcf_all$INFO, ";") %>% map(str_remove_all, "=.*$") 
mutect_info <- str_split(mutect_vcf_all$INFO, ";") %>% map(str_remove_all, "^.*\\=") 
mutect_info_df <- map2(mutect_info_names, mutect_info, names_set) %>% 
        bind_rows %>% 
        select(-c( "DP", "ECNT", "RPA", "RU", "STR", "TLOD"))

mutect_data_names <- str_split(mutect_vcf_all$FORMAT, ":")  
mutect_data <- str_split(mutect_vcf_all$DATA, ":") 
mutect_data_df <- map2(mutect_data_names, mutect_data, names_set) %>% 
        bind_rows %>% 
        select(-DP)

mutect_vcf_bind <- select(mutect_vcf_all, Sample:FILTER) %>% 
         bind_cols(mutect_data_df) %>% 
         bind_cols(mutect_info_df) 

# Extract the maf data using the make_data function and bind the resulting data into one tibble
maf_columns <- c("Hugo_Symbol", "NCBI_Build", "Chromosome", "Start_Position", "End_Position", "Variant_Classification", "Variant_Type", "Protein_Change",
                 "tumor_f", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "Transcript_Exon", "Transcript_Position", "cDNA_Change", "Codon_Change", "gc_content", "DP", "MBQ", "MFRL", "MMQ","AS_SB_TABLE", 
                 "AS_FilterStatus", "ECNT", "GERMQ", "MPOS", "POPAF", "TLOD", "RPA", "RU", "STR", "STRQ", "OREGANNO_ID", "OREGANNO_Values", "Other_Transcripts","ref_context")
maf_column_types <- c("character", "character", "character", "numeric", "numeric", "factor", "factor", "character", "numeric", "character", "character", "character", "integer", 
                      "character", "character", "character", "numeric", "integer", "character", "character", "character", "character", "character", "integer", "integer", 
                      "character", "character", "character", "character", "character", "logical", "integer", "character", "character", "character", "character") 

list_of_mafs <- maf_files %>%
                map(read.maf) %>%
                map(make_data) %>%
                map(as_tibble) %>%
                set_names(sample_names) %>%
                map(select, maf_columns) %>%
                map(map2, paste0("as.", maf_column_types), ~ get(.y)(.x))

mutect_maf_all <- bind_rows(list_of_mafs, .id = "Sample")

# Remove all mitochondrial genes
maf_noMT <- subset(mutect_maf_all, mutect_maf_all$Chromosome != "MT")

# IF a gene's Reference_Allele column contains an A, T, G, or C AND their Tumor_Seq_Allele2 column is blank ("-")
    # THEN, subtract one from the Start_Position column
TSA_deletion <- filter(maf_noMT, Reference_Allele != "-" & Tumor_Seq_Allele2 == "-")
Start_Position_Mutated <- mutate(TSA_deletion, Start_Position = as.integer(Start_Position)) %>%
            mutate(VCF_Start_Position = Start_Position - 1)
# IF a gene's Reference_Allele, Tumor_Seq_Allele1, and Tumor_Seq_Allele2 columns contain an A, T, G, or C
    # THEN, do not modify the start position column
TSA_snp <- filter(maf_noMT, Reference_Allele != "-" & Tumor_Seq_Allele1 != "-" & Tumor_Seq_Allele2 != "-") %>%
            mutate(VCF_Start_Position = Start_Position)
# IF a gene's Reference_Allele column contains is blank ("-") 
    # THEN, do not modify the start position column
TSA_insertion <- filter(maf_noMT, Reference_Allele == "-") %>%
            mutate(VCF_Start_Position = Start_Position)
# Combine all partitioned dataframes
Modified_maf <- rbind(TSA_snp, TSA_insertion, Start_Position_Mutated)

# Combine a pruned version of the vcf tibble to the maf tibble 
vcf <- unite(mutect_vcf_bind, "Chrom_Pos", c("#CHROM", "POS", "Sample"), remove = FALSE)
maf <- unite(Modified_maf, "Chrom_Pos", c("Chromosome", "VCF_Start_Position", "Sample"), remove = FALSE)
vcf_pruned <- select(vcf, c(Chrom_Pos, FILTER, GT, AD, AF, F1R2, F2R1, PGT, PS, SB, PID))
combined <- vcf_pruned %>%
        merge(maf, by = "Chrom_Pos") %>%
        select(-Chrom_Pos)

# Apply the Detect_Repeat function to every row to find the longest repeat for every variant
combined$longest_repeat <- map_int(combined$ref_context, detect_repeat)

# Filter out all variant's based on Variant_Classification; remove all rows blank Protein_Change values
variant_classification <- filter(combined, is_in(Variant_Classification, c("Frame_Shift_Del", "Frame_Shift_Ins", "Missense_Mutation", "Nonsense_Mutation" ,"Splice_Site"))) %>% 
  filter(nchar(Protein_Change) > 0)

# Apply the Split_Columns function to AD, tumor_f, F2R1, F1R2, MBQ, MMQ, and MFRL columns
mutect_vcf_ncol <- str_split(variant_classification$AD, ",") %>% 
            map_int(length) %>% 
            max

mutect_vcf_ad <- split_columns(variant_classification$AD, mutect_vcf_ncol, "t_ref_count", "t_alt_count_")
mutect_vcf_af <- split_columns(variant_classification$tumor_f, mutect_vcf_ncol, "tumor_f", "")
mutect_vcf_f2r1 <- split_columns(variant_classification$F2R1, mutect_vcf_ncol, "f2r1_reference", "f2r1_alternate")
mutect_vcf_f1r2 <- split_columns(variant_classification$F1R2, mutect_vcf_ncol, "f1r2_reference", "f1r2_alternate")
mutect_vcf_mbq <- split_columns(variant_classification$MBQ, mutect_vcf_ncol, "mbq_reference", "mbq_alternate")
mutect_vcf_mmq <- split_columns(variant_classification$MMQ, mutect_vcf_ncol, "mmq_reference", "mmq_alternate")
mutect_vcf_mfrl <- split_columns(variant_classification$MFRL, mutect_vcf_ncol, "mfrl_reference", "mfrl_alternate")
mutect_vcf_sb <- str_replace_all(variant_classification$AS_SB_TABLE, "\\[|\\]", "") %>% str_split_fixed( "\\|", mutect_vcf_ncol) %>% as_tibble 
colnames(mutect_vcf_sb) <- c("sb_reference", str_c("sb_alt", seq(1:(mutect_vcf_ncol - 1))))

# Remove all rows with a t_ref_count == "0"
mutect_vcf_filter <- select(variant_classification, -AD, -tumor_f, -F2R1, -F1R2, -MBQ, -MMQ, -MFRL, -AS_SB_TABLE) %>% 
  bind_cols(mutect_vcf_ad, mutect_vcf_af, mutect_vcf_f2r1, mutect_vcf_f1r2, mutect_vcf_mbq, mutect_vcf_mfrl, mutect_vcf_mmq, mutect_vcf_sb) %>%
  filter(t_ref_count != "0")

# Select relevant columns for the final output
mutect_vcf_select <- select(mutect_vcf_filter, c(Sample, Hugo_Symbol, NCBI_Build, Chromosome, Start_Position,
                                                 End_Position, Variant_Classification, Variant_Type, Protein_Change,
                                                 FILTER, tumor_f, t_ref_count, t_alt_count_1, t_alt_count_2,
                                                 Reference_Allele, Tumor_Seq_Allele1, Transcript_Exon, Transcript_Position,
                                                 cDNA_Change, Codon_Change, gc_content, longest_repeat, DP, f1r2_reference,
                                                 gc_content, longest_repeat, f1r2_reference,
                                                 f1r2_alternate1, f1r2_alternate2,f2r1_reference,f2r1_alternate1,      
                                                 f2r1_alternate2,mbq_reference,mbq_alternate1,mbq_alternate2,  
                                                 mfrl_reference, mfrl_alternate1, mfrl_alternate2, mmq_reference,        
                                                 mmq_alternate1, mmq_alternate2, sb_reference,sb_alt1,               
                                                 sb_alt2, AS_FilterStatus,ECNT,GERMQ, MPOS,POPAF,TLOD, RPA, RU,STR,STRQ,GT,                  
                                                 PGT,PID,PS,OREGANNO_ID, OREGANNO_Values,Other_Transcripts,ref_context           
))

# Get twist panel (in same folder as aggregate variants mutect script)
twist_panel <- read_excel(panel_coordinates, col_names = F)

colnames(twist_panel) <- c("chr", "start", "end", "Transcript", "X5", "Strand", "Gene", "X8")
twist_panel_granges <- select(twist_panel, chr:end) %>% makeGRangesFromDataFrame

mutect_vcf_granges <- select(mutect_vcf_select, Chromosome:End_Position) %>% 
            set_colnames(c("chr", "start", "end")) %>% 
            makeGRangesFromDataFrame
mutect_overlaps <- findOverlaps(twist_panel_granges, mutect_vcf_granges) %>% as_tibble
mutect_overlaps_sorted <- sort(mutect_overlaps$subjectHits)

# filter for all variants within the genomic ranges specified in the twist_panel
mutect_vcf_filter <- slice(mutect_vcf_select, mutect_overlaps_sorted) %>% 
            filter(is_in(Hugo_Symbol, unique(twist_panel$Gene))) %>% 
            arrange(Sample, Chromosome, Start_Position)

mutect_basic_file <- str_c(mutect_directory, "/mutect_aggregated_simple.tsv")

# write output into tsv file
write_tsv(mutect_vcf_filter, mutect_basic_file)
