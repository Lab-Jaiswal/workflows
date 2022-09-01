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
panel <- command_args[3]
filtered <- command_args[4]

#if we filtered the vcf and maf files before starting the R script, then the files have the ending "funcotator_coding.maf or funcotator_coding.vcf"
#I have been using filtered <- 0 for this data
if (filtered == 1){
    vcf_pattern="*_funcotator_coding.vcf$"
    maf_pattern="*_funcotator_coding.maf$"
} else {
    vcf_pattern="*_funcotator.vcf$"
    maf_pattern="*_funcotator.maf$"
    }

#get a list of the vcf and maf funcotator file names
split_names_vcf <- list.files(mutect_directory, pattern = vcf_pattern) %>% str_remove_all("_.*$")
split_names_maf <- list.files(mutect_directory, pattern = maf_pattern) %>% str_remove_all("_.*$")
if (length(unique(split_names_vcf)) >= 1){
    sample_names_vcf <- list.files(mutect_directory, pattern = vcf_pattern) %>% str_remove_all("_G.*$")
    sample_names_maf <- list.files(mutect_directory, pattern = maf_pattern) %>% str_remove_all("_G.*$")

    } else {
        sample_names_vcf <- split_names_vcf
        sample_names_maf <- split_names_maf
    }

maf_files <- list.files(mutect_directory, pattern = maf_pattern, full.names = TRUE)
vcf_files <- list.files(mutect_directory, pattern = vcf_pattern, full.names = TRUE)

# Extract the vcf files, as well as the information in the INFO, FORMAT, and DATA columns
# Bind all of the newly formatted vcf files together horizontally
# Get the column names from the header of one of the files
mutect_vcf_header <- read_lines(vcf_files[1]) 

vcf_colnames <- str_subset(mutect_vcf_header, "^#") %>% 
  str_subset("^##", negate = T) %>% 
  str_split("\\t") %>% 
  extract2(1)
vcf_colnames[length(vcf_colnames)] <- "DATA"
mutect_vcf_list <- map(vcf_files, read_lines) %>%
  map(str_subset, pattern = "^#", negate = TRUE) %>%
  map(str_split_fixed, "\\t", length(vcf_colnames)) %>% 
  map(set_colnames, vcf_colnames) %>% 
  map(as_tibble) %>% 
  set_names(sample_names_vcf)
mutect_vcf_all <- bind_rows(mutect_vcf_list, .id = "Sample")

#mutect_info_names <- str_split(mutect_vcf_all$INFO, ";") %>% map(str_remove_all, "=.*$") 
#mutect_info <- str_split(mutect_vcf_all$INFO, ";") %>% map(str_remove_all, "^.*\\=") 
#mutect_info_df <- map2(mutect_info_names, mutect_info, names_set) %>% 
#  bind_rows %>% 
#  select(-c( "DP", "ECNT", "RPA", "RU", "STR", "TLOD"))

mutect_data_names <- str_split(mutect_vcf_all$FORMAT, ":")  
mutect_data <- str_split(mutect_vcf_all$DATA, ":") 
mutect_data_df <- map2(mutect_data_names, mutect_data, names_set) %>% 
  bind_rows %>% 
  select(-DP)

mutect_vcf_bind <- select(mutect_vcf_all, Sample:FILTER) %>% 
  bind_cols(mutect_data_df) #%>% 
  #bind_cols(mutect_info_df) 

# Extract the maf data using the make_data function and bind the resulting data into one tibble
maf_columns <- c("Hugo_Symbol", "NCBI_Build", "Chromosome", "Start_Position", "End_Position", "Variant_Classification", "Variant_Type", "Protein_Change", "tumor_f", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "Transcript_Exon", "Transcript_Position", "cDNA_Change", "Codon_Change", "gc_content", "DP", "MBQ", "MFRL", "MMQ","AS_SB_TABLE", "AS_FilterStatus", "ECNT", "GERMQ", "MPOS", "POPAF", "TLOD",  "STRQ", "OREGANNO_ID", "OREGANNO_Values", "Other_Transcripts","ref_context", "t_ref_count", "t_alt_count", "AS_UNIQ_ALT_READ_COUNT", "Strand", "Genome_Change", "Annotation_Transcript", "Transcript_Strand", "Refseq_mRNA_Id", "ref_context", "t_alt_count", "t_ref_count", "n_alt_count", "n_ref_count", "AS_SB_TABLE")

list_of_mafs <- maf_files %>%
    map(read.maf) %>%
    map(make_data) %>%
    map(as_tibble) %>%
    set_names(sample_names_maf) %>%
    map(select, all_of(maf_columns)) 

#no non-synomous mutations found
list_of_mafs_edited <- lapply(list_of_mafs, function(df) mutate_at(df, .vars = c("Transcript_Position", "Transcript_Exon", "gc_content", "MPOS", "POPAF", "TLOD"), as.character))

mutect_maf_all <- bind_rows(list_of_mafs_edited, .id = "Sample")

# Remove all mitochondrial genes

# IF a gene's Reference_Allele column contains an A, T, G, or C AND their Tumor_Seq_Allele2 column is blank ("-")
# THEN, subtract one from the Start_Position column
TSA_deletion <- filter(mutect_maf_all, Reference_Allele != "-" & Tumor_Seq_Allele2 == "-")
Start_Position_Mutated <- mutate(TSA_deletion, Start_Position = as.integer(Start_Position)) %>%
  mutate(VCF_Start_Position = Start_Position - 1)
# IF a gene's Reference_Allele, Tumor_Seq_Allele1, and Tumor_Seq_Allele2 columns contain an A, T, G, or C
# THEN, do not modify the start position column
TSA_snp <- filter(mutect_maf_all, Reference_Allele != "-" & Tumor_Seq_Allele1 != "-" & Tumor_Seq_Allele2 != "-") %>%
  mutate(VCF_Start_Position = Start_Position)
# IF a gene's Reference_Allele column contains is blank ("-") 
# THEN, do not modify the start position column
TSA_insertion <- filter(mutect_maf_all, Reference_Allele == "-") %>%
  mutate(VCF_Start_Position = Start_Position)
# Combine all partitioned dataframes
Modified_maf <- rbind(TSA_snp, TSA_insertion, Start_Position_Mutated)

# Combine a pruned version of the vcf tibble to the maf tibble 
vcf <- unite(mutect_vcf_bind, "Chrom_Pos", c("#CHROM", "POS", "Sample"), remove = FALSE)
maf <- unite(Modified_maf, "Chrom_Pos", c("Chromosome", "VCF_Start_Position", "Sample"), remove = FALSE)

columns <- c("Chrom_Pos", "FILTER", "AF","F1R2", "F2R1")
pruned_columns <- intersect(colnames(vcf), columns) 
vcf_pruned <- select(vcf, all_of(pruned_columns))
colnames(vcf_pruned) <- c("Chrom_Pos", "vcf_FILTER", "vcf_AF", "vcf_F1R2", "vcf_F2R1")
combined <- vcf_pruned %>%
  merge(maf, by = "Chrom_Pos") %>%
  select(-Chrom_Pos)

# Apply the Detect_Repeat function to every row to find the longest repeat for every variant
combined$longest_repeat <- map_int(combined$ref_context, detect_repeat)

# Filter out all variant's based on Variant_Classification; remove all rows blank Protein_Change values
if (filtered == 1){
    variant_classification <- filter(combined, is_in(Variant_Classification, c("Frame_Shift_Del", "Frame_Shift_Ins", "Missense_Mutation", "Nonsense_Mutation" ,"Splice_Site"))) %>% filter(nchar(Protein_Change) > 0)
    maf_noMT <- subset(mutect_maf_all, mutect_maf_all$Chromosome != "MT")

} else {
    variant_classification <- combined
    }

# Apply the Split_Columns function to AD, tumor_f, F2R1, F1R2, MBQ, MMQ, and MFRL columns
max_ncol <- function(argument) {
    number_args <- str_split(argument, ",") %>%
        map_int(length) %>%
        max
           }
  
#AD_vcf_ncol <- max(map_int(variant_classification$AD, max_ncol)) 
#mutect_vcf_ad <- split_columns(variant_classification$AD, AD_vcf_ncol, "t_ref_count", "t_alt_count_")

#tumor_f_vcf_ncol <- max(map_int(variant_classification$tumor_f, max_ncol)) 
#mutect_vcf_af <- split_columns(variant_classification$tumor_f, tumor_f_vcf_ncol, "maf_tumor_f", "")

F2R1_vcf_ncol <- max(map_int(variant_classification$vcf_F2R1, max_ncol)) 
mutect_vcf_f2r1 <- split_columns(variant_classification$vcf_F2R1, F2R1_vcf_ncol, "vcf_f2r1_reference", "vcf_f2r1_alternate")

F1R2_vcf_ncol <- max(map_int(variant_classification$vcf_F1R2, max_ncol)) 
mutect_vcf_f1r2 <- split_columns(variant_classification$vcf_F1R2, F1R2_vcf_ncol, "vcf_f1r2_reference", "vcf_f1r2_alternate")

MBQ_vcf_ncol <- max(map_int(variant_classification$MBQ, max_ncol)) 
mutect_vcf_mbq <- split_columns(variant_classification$MBQ, MBQ_vcf_ncol, "maf_mbq_reference", "maf_mbq_alternate")

MMQ_vcf_ncol <- max(map_int(variant_classification$MMQ, max_ncol)) 
mutect_vcf_mmq <- split_columns(variant_classification$MMQ, MMQ_vcf_ncol, "maf_mmq_reference", "maf_mmq_alternate")

MFRL_vcf_ncol <- max(map_int(variant_classification$MFRL, max_ncol)) 
mutect_vcf_mfrl <- split_columns(variant_classification$MFRL, MFRL_vcf_ncol, "maf_mfrl_reference", "maf_mfrl_alternate")

mutect_vcf_sb <- str_remove_all(variant_classification$AS_SB_TABLE, "\\[|\\]") %>% str_split_fixed( "\\|", 3) %>% as_tibble 
colnames(mutect_vcf_sb) <- c("maf_sb_reference", str_c("maf_sb_alt", seq(1:(3 - 1))))

# Remove all rows with a t_ref_count == "0"
mutect_vcf_select <- select(variant_classification,  -vcf_F2R1, -vcf_F1R2, -MBQ, -MMQ, -MFRL) %>% 
  bind_cols(mutect_vcf_f2r1, mutect_vcf_f1r2, mutect_vcf_mbq, mutect_vcf_mfrl, mutect_vcf_mmq, mutect_vcf_sb)


if (filtered == 1){
    mutect_vcf_filter <- mutect_vcf_select %>% filter(t_ref_count != "0")
} else {
    mutect_vcf_filter <- mutect_vcf_select
    }
  
# Select relevant columns for the final output
final_columns <- c("Sample", "Hugo_Symbol", "NCBI_Build", "Chromosome", "Start_Position", "End_Position", "Strand", "Variant_Classification", "Variant_Type", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "Genome_Change", "Annotation_Transcript", "Transcript_Strand", "Transcript_Exon", "Transcript_Position", "cDNA_Change", "Codon_Change","Protein_Change","Refseq_mRNA_Id", "ref_context", "tumor_f", "t_alt_count", "t_ref_count", "n_alt_count", "n_ref_count", "AS_FilterStatus", "vcf_FILTER", "AS_SB_TABLE", "TLOD",  "gc_content", "longest_repeat", "DP", "f1r2_reference",  "f1r2_reference", "f2r1_reference", "mbq_reference", "mfrl_reference",  "mmq_reference", "sb_reference", "AS_FilterStatus", "ECNT", "GERMQ", "MPOS", "POPAF", "TLOD", "OREGANNO_ID", "OREGANNO_Values", "Other_Transcripts")

pruned_final_columns <- intersect(colnames(mutect_vcf_filter), final_columns) 

mutect_vcf_select <- select(mutect_vcf_filter, c(all_of(pruned_final_columns), contains("t_alt_count"), contains("f1r2_alternate"), contains("sb_alt"), contains("mfrl_alternate"), contains("mmq_alternate"), contains("f2r1_alternate"), contains("mbq_alternate")))

directory <- dirname(mutect_directory)
mutect_simple_file <- str_c("/home/maurertm/smontgom/maurertm", "/mutect_aggregated_noWL_Apr25.tsv")

# write output into tsv file
write_tsv(mutect_vcf_select, mutect_simple_file)

# Get twist panel (in same folder as aggregate variants mutect script)
if (panel != "false") {

twist_panel <- read_excel(panel_coordinates, col_names = F)

colnames(twist_panel) <- c("chr", "start", "end", "Transcript", "X5", "Strand", "Gene", "X8")
twist_panel_granges <- select(twist_panel, chr:end) %>% makeGRangesFromDataFrame

mutect_vcf_granges <- select(mutect_vcf_select, Chromosome:End_Position) %>% 
  set_colnames(c("chr", "start", "end")) %>% 
  makeGRangesFromDataFrame
mutect_overlaps <- findOverlaps(twist_panel_granges, mutect_vcf_granges) %>% as_tibble
mutect_overlaps_sorted <- sort(mutect_overlaps$subjectHits)

## filter for all variants within the genomic ranges specified in the twist_panel
mutect_vcf_filter <- slice(mutect_vcf_select, mutect_overlaps_sorted) %>% 
  filter(is_in(Hugo_Symbol, unique(twist_panel$Gene))) %>% 
  arrange(Sample, Chromosome, Start_Position)

mutect_twist_file <- str_c(mutect_directory, "/mutect_aggregated_twist.tsv")

## write output into tsv file
write_tsv(mutect_vcf_filter, mutect_twist_file)

} else {
    print("Twist analysis not requested. Please add --twist to your ./submit_BWA_CHIP command to filter your results by the Twist panel.")
}
