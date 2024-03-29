library(GenomicRanges)
library(readxl)
library(magrittr)
library(tidyverse)

# Function input: row of semicolon separated variant effects
# Function goal: split the variant effects into separate columns,
# and each effect on its own line
MakeTibble <- function(rows) {
    # If there is only one variant effect, only make one row. Otherwise, return
    # a tibble of all variant effects.
    if (is.null(dim(rows))) {
        row_df <- data.frame("Gene" = NA, "Transcript" = NA, "Exon" = NA, "DNA_Change" = NA, "Protein_change" = NA)
        as_tibble_row(row_df, .name_repair = "minimal")
        
    } else {
        row_df <- set_colnames(rows, c("Gene", "Transcript", "Exon", "DNA_change", "Protein_change")) %>% 
            as_tibble
    }
    row_df
}

command_args <- commandArgs(trailingOnly = TRUE)
panel_coordinates <- command_args[1]
varscan_directory <- command_args[2]

varscan_files <- list.files(varscan_directory, pattern = "*multianno.txt", full.names = TRUE, recursive = TRUE)
varscan_data_list <- map(varscan_files, read_tsv, skip = 1, col_names = F)

sample_names <- list.files(varscan_directory, pattern = "*multianno.txt", recursive = TRUE) %>% str_remove_all("_.*$")
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
varscan_nofilter_file <- str_c(varscan_directory, "/varscan_aggregated_nofilter.tsv")
varscan_nofilter <- varscan_final
varscan_aachange_nofilter <- str_split(varscan_nofilter$AAChange.ensGene, ",") %>%
    map(str_split, ":") %>% map(reduce, rbind) %>% map(MakeTibble) 

varscan_nofilter$AAChange_split <- varscan_aachange_nofilter
# Use the unnest command to create separate rows
varscan_nofilter %<>% select(-AAChange.ensGene) %>% unnest(AAChange_split)
write_tsv(varscan_nofilter, varscan_nofilter_file)

# Get twist panel (in same folder as aggregate variants mutect script)
twist_panel <- read_excel(panel_coordinates, col_names = F)

colnames(twist_panel) <- c("chr", "start", "end", "Transcript", "X5", "Strand", "Gene", "X8")
twist_panel_granges <- select(twist_panel, chr:end) %>% makeGRangesFromDataFrame

varscan_granges <- select(varscan_final, Chr:End) %>% 
  set_colnames(c("chr", "start", "end")) %>% 
  makeGRangesFromDataFrame
varscan_overlaps <- findOverlaps(twist_panel_granges, varscan_granges) %>% as_tibble
varscan_overlaps_sorted <- sort(varscan_overlaps$subjectHits)

# filter for all variants within the genomic ranges specified in the twist_panel
varscan_vcf_filter <- dplyr::slice(varscan_final, varscan_overlaps_sorted) %>% 
  dplyr::filter(is_in(Gene.ensGene, unique(twist_panel$Gene))) %>% 
  arrange(Sample, Chr, Start)

varscan_aachange <- str_split(varscan_vcf_filter$AAChange.ensGene, ",") %>%
    map(str_split, ":") %>% map(reduce, rbind) %>% map(MakeTibble) 

varscan_vcf_filter$AAChange_split <- varscan_aachange
# Use the unnest command to create separate rows
varscan_vcf_filter %<>% select(-AAChange.ensGene) %>% unnest(AAChange_split)

# write output into tsv file
varscan_final_file <- str_c(varscan_directory, "/varscan_aggregated.tsv")
write_tsv(varscan_vcf_filter, varscan_final_file)
