library(GenomicRanges)
library(magrittr)
library(tidyverse)

command_args <- commandArgs(trailingOnly = TRUE)
panel_coordinates <- command_args[1]
varscan_directory <- command_args[2]

varscan_files <- list.files(varscan_directory, pattern = "*multianno.txt", full.names = TRUE)
varscan_data_list <- map(varscan_files, read_tsv, skip = 1, col_names = F)

sample_names <- list.files(varscan_directory, pattern = "*multianno.txt") %>% str_remove_all("_.*$")
names(varscan_data_list) <- sample_names

varscan_data <- bind_rows(varscan_data_list, .id = "Sample")
varscan_data_basic <- select(varscan_data, Sample, X1:X10)
colnames(varscan_data_basic) <- c("Sample", "Chr", "Start", "End", "Ref", "Alt", "Func.ensGene", "Gene.ensGene", "GeneDetail.ensGene", "ExonicFunc.ensGene", "AAChange.ensGene")
varscan_data_basic$ExonicFunc.ensGene %<>% str_replace_all(" ", "_")

varscan_data_other <- select(varscan_data, X11:X21)
colnames(varscan_data_other) <- c("Genotype", "NC", "ADP", "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")

varscan_info_names <- unique(varscan_data$X22) %>% str_split(":") %>% extract2(1)
varscan_info <- select(varscan_data, X23) %>% separate(X23, sep = ":", into = varscan_info_names)
varscan_data_basic$VAF <- str_remove_all(varscan_info$FREQ, "%") %>% as.numeric %>% divide_by(100) %>% signif(digits = 5)

varscan_final <- bind_cols(varscan_data_basic, varscan_data_other, varscan_info)

# Get twist panel (in same folder as aggregate variants mutect script)
twist_panel <- read_excel(panel_coordinates, col_names = F)

colnames(twist_panel) <- c("chr", "start", "end", "Transcript", "X5", "Strand", "Gene", "X8")
twist_panel_granges <- select(twist_panel, chr:end) %>% makeGRangesFromDataFrame

varscan_granges <- select(varscan_final, Chromosome:End_Position) %>% 
  set_colnames(c("chr", "start", "end")) %>% 
  makeGRangesFromDataFrame
varscan_overlaps <- findOverlaps(twist_panel_granges, varscan_granges) %>% as_tibble
varscan_overlaps_sorted <- sort(varscan_overlaps$subjectHits)

# filter for all variants within the genomic ranges specified in the twist_panel
varscan_vcf_filter <- dplyr::slice(varscan_final, varscan_overlaps_sorted) %>% 
  dplyr::filter(is_in(Hugo_Symbol, unique(twist_panel$Gene))) %>% 
  arrange(Sample, Chromosome, Start_Position)

# write output into tsv file
varscan_final_file <- str_c(varscan_directory, "/varscan_aggregated.txt")
write_tsv(varscan_vcf_filter, varscan_final_file)
