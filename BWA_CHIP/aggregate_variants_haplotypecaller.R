library(GenomicRanges)
library(readxl)
library(magrittr)
library(tidyverse)
library(dplyr)

# Function input: col_names and col_values as arguments
# Function goal: assigns the names of the col_values to be the col_names
names_set <- function(col_names, col_values) {
  names(col_values) <- col_names
  col_values
}

# This function is not very elegant - each row of the INFO, FORMAT columns can have a different number of fields 
# We parse each row separately as a vector which we assign names to.  Then we use bind_rows to automatically figure out the right column names
# In subject with homozygous reference genotypes, most of the columns will be empty because they have no meaning when no variant is present
ParseHaplotypeCallerVCF <- function(haplotypecaller_vcf) {
  haplotypecaller_info_names <- str_split(haplotypecaller_vcf$INFO, ";") %>% map(str_remove_all, "=.*$") 
  haplotypecaller_info <- str_split(haplotypecaller_vcf$INFO, ";") %>% map(str_remove_all, "^.*\\=") 
  haplotypecaller_info_df <- map2(haplotypecaller_info_names, haplotypecaller_info, names_set) %>% bind_rows
  haplotypecaller_info_nums <- mutate(haplotypecaller_info_df, across(everything(), as.numeric)) %>% select(-DP)
  haplotypecaller_data_names <- str_split(haplotypecaller_vcf$FORMAT, ":")  
  haplotypecaller_data <- str_split(haplotypecaller_vcf$DATA, ":") 
  haplotypecaller_data_df <- map2(haplotypecaller_data_names, haplotypecaller_data, names_set) %>% bind_rows
  haplotypecaller_vcf_bind <- select(haplotypecaller_vcf, `#CHROM`:FILTER) %>% bind_cols(haplotypecaller_data_df) %>% bind_cols(haplotypecaller_info_nums)
  haplotypecaller_vcf_bind
}

# Define the location of the FastQs and the twist panel
command_args <- commandArgs(trailingOnly = TRUE)
panel_coordinates <- command_args[1]
mutect_directory <- command_args[2]

sample_names <- list.files(mutect_directory, pattern = "*_funcotator.vcf$") %>% str_remove_all("_.*$")

# Parse germline genotype VCFs from HaplotypeCaller
haplotypecaller_vcf_files <- list.files(mutect_directory, pattern = "*_haplotypecaller_genotypes.vcf$", full.names = TRUE)
haplotypecaller_vcf_header <- read_lines(haplotypecaller_vcf_files[1]) # Read first VCF to get header 
haplotypecaller_vcf_colnames <- str_subset(haplotypecaller_vcf_header, "^#") %>% # Only keep rows that begin with # because they are from the header
  str_subset("^##", negate = T) %>% str_split("\\t") %>% extract2(1) # Only keep line with a single # at the beginning because that is the column name
haplotypecaller_vcf_colnames[length(haplotypecaller_vcf_colnames)] <- "DATA" # Rename column that is sample specific - we will get the same name a different 
haplotypecaller_vcf_list <- map(haplotypecaller_vcf_files, read_lines) %>% 
  map(str_subset, "^#", negate = TRUE) %>% 
  map(str_split_fixed, "\\t", length(haplotypecaller_vcf_colnames)) %>% 
  map(set_colnames, haplotypecaller_vcf_colnames) %>% map(as_tibble) %>% set_names(sample_names)

haplotypecaller_vcf_parsed <- map(haplotypecaller_vcf_list, ParseHaplotypeCallerVCF) %>% set_names(sample_names) %>% bind_rows(.id = "Sample")
haplotypecaller_vcf_parsed$AC %<>% replace_na(0)
haplotypecaller_vcf_parsed$AF %<>% replace_na(0)
haplotypecaller_vcf_parsed$AN %<>% replace_na(1)

haplotypecaller_vcf_granges <- select(haplotypecaller_vcf_parsed, `#CHROM`:POS) %>% set_colnames(c("chr", "start")) %>% mutate(end = start) %>% makeGRangesFromDataFrame
mcols(haplotypecaller_vcf_granges) <- select(haplotypecaller_vcf_parsed, -`#CHROM`, -POS)

twist_panel <- read_excel(panel_coordinates, col_names = F)
colnames(twist_panel) <- c("chr", "start", "end", "Transcript", "X5", "Strand", "Gene", "X8")

twist_snps <- filter(twist_panel, str_detect(Gene, "^rs"))  
colnames(twist_snps)[7] <- "rsID"
twist_snps$rsID %<>% str_trim
colnames(twist_snps)[8] <- "Gene"
twist_granges <- select(twist_snps, chr:end) %>% makeGRangesFromDataFrame
mcols(twist_granges) <- select(twist_snps, rsID, Gene)

haplotypecaller_vcf_annot <- mergeByOverlaps(haplotypecaller_vcf_granges, twist_granges) %>% as_tibble %>% select(Sample, Gene, rsID, haplotypecaller_vcf_granges.seqnames, haplotypecaller_vcf_granges.start, REF, ALT, GT, AC, AF, AD:PL, QUAL, AN:RGQ) %>% arrange(Gene, Sample)
colnames(haplotypecaller_vcf_annot) %<>% str_replace_all("haplotypecaller_vcf_granges.seqnames", "Chromosome") %>%
  str_replace_all("haplotypecaller_vcf_granges.start", "Position") %>%
  str_replace_all("REF", "Reference_Allele") %>%
  str_replace_all("ALT", "Alternate_Allele") %>%
  str_replace_all("GT", "Genotype") %>%
  str_replace_all("AC", "Allele_Count") %>%
  str_replace_all("AF", "Allele_Frequency") %>%
  str_replace_all("AD", "Allele_Depths") %>%
  str_replace_all("DP", "Sequencing_Depth") 

haplotypecaller_vcf_annot$Genotype %<>% str_replace_all("\\/", "|")

# You have to explicitly call dplyr::slice because S4Vectors (part of GenomicRanges) has a different slice command
haplotypecaller_vcf_annot_filter <- arrange(haplotypecaller_vcf_annot, Sample, rsID) %>% group_by(rsID, Sample) %>% dplyr::slice(1)

haplotypecaller_germline_file <- str_c(mutect_directory, "/haplotypecaller_germline_genotypes.tsv")
write_tsv(haplotypecaller_vcf_annot_filter, haplotypecaller_germline_file)
