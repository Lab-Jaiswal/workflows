library(argparser)
library(pbmcapply)
library(vroom)
library(readxl)
library(magrittr)
library(tidyverse)

subset_ranged <- function(grouped_df, grouping_df, mutect_df) {
    mutect_subset_ranged <- filter(mutect_df, is_in(gencode_28_hugosymbol, grouping_df$hgnc_symbol) & start_position >= grouping_df$start & end_position <= grouping_df$end)
    mutect_subset_ranged
}

subset_protein_change <- function(grouped_df, grouping_df, mutect_df) {
    grouped_df_noins <- filter(grouped_df, str_detect(protein_change, "ins", negate = TRUE))
    mutect_subset_protein_change <- filter(mutect_df, is_in(gencode_28_hugosymbol, grouping_df$hgnc_symbol) & is_in(gencode_28_proteinchange, grouped_df_noins$protein_change))

    grouped_df_ins <- filter(grouped_df, str_detect(protein_change, "ins"))
    if (nrow(grouped_df_ins) > 0) {
        grouped_df_ins$protein_change <- str_c(grouped_df_ins$protein_change, "*")
        mutect_subset_protein_change_ins <- filter(mutect_df, is_in(gencode_28_hugosymbol, grouping_df$hgnc_symbol) & str_detect(gencode_28_proteinchange, str_c(grouped_df_ins$protein_change, collapse = "|")))
        if (nrow(mutect_subset_protein_change_ins) > 0) {
            mutect_subset_protein_change <- bind_rows(mutect_subset_protein_change_ins, mutect_subset_protein_change)
        }
    }
    mutect_subset_protein_change
}

whitelist_annotate <- function(grouped_df, grouping_df, mutect_df) {
    mutect_subset <- filter(mutect_df, consequence_reformat == grouping_df$variant_type)

    whitelist_no_pos_no_pchange <- filter(grouped_df, is.na(start) & is.na(protein_change))
    mutect_subset_whitelist <- tibble()
    if (nrow(whitelist_no_pos_no_pchange) > 0) {
        mutect_subset_no_pos_no_pchange <- filter(mutect_subset, is_in(gencode_28_hugosymbol, whitelist_no_pos_no_pchange$hgnc_symbol))
        if (nrow(mutect_subset_no_pos_no_pchange) > 0) {
            mutect_subset_whitelist <- bind_rows(mutect_subset_whitelist, mutect_subset_no_pos_no_pchange)
        }
        mutect_subset <- filter(mutect_subset, !is_in(gencode_28_hugosymbol, whitelist_no_pos_no_pchange$hgnc_symbol))
    } 

    whitelist_no_pchange <- filter(grouped_df, !is.na(start) & is.na(protein_change))
    if (nrow(whitelist_no_pchange) > 0) {
        mutect_subset_no_pchange <- group_by(whitelist_no_pchange, hgnc_symbol, start, end) %>% 
            group_modify(subset_ranged, mutect_subset) %>%
            ungroup %>%
            select(-hgnc_symbol, -start, -end)
        if (nrow(mutect_subset_no_pchange) > 0) {
            mutect_subset_whitelist <- bind_rows(mutect_subset_whitelist, mutect_subset_no_pchange)
        }
        mutect_subset <- setdiff(mutect_subset, mutect_subset_no_pchange)
    }

    whitelist_no_pos <- filter(grouped_df, is.na(start) & !is.na(protein_change))
    if (nrow(whitelist_no_pos) > 0) {
        mutect_subset_no_pos <- group_by(whitelist_no_pos, hgnc_symbol) %>% group_modify(subset_protein_change, mutect_subset) %>%
            ungroup %>%
            select(-hgnc_symbol)
        if (nrow(mutect_subset_no_pos) > 0) {
            mutect_subset_whitelist <- bind_rows(mutect_subset_whitelist, mutect_subset_no_pos)
        }
        mutect_subset <- setdiff(mutect_subset, mutect_subset_no_pos)
    }
    mutect_subset_whitelist$whitelist <- 1
    mutect_subset$whitelist <- 0
    mutect_combined <- bind_rows(mutect_subset_whitelist, mutect_subset)
    mutect_combined
}

parser <- arg_parser("Script to annotate and filter CHIP calls") %>%
    add_argument("--output_directory", help = "Output directory") %>%
    add_argument("--chip_calls", help = "TSV of aggregated CHIP calls across saples") %>%
    add_argument("--pileup_regions", help = "TSV of aggregated calls for U2AF1 pileup_regions") %>%
    add_argument("--whitelist", help = "Excel file of variant whitelist") %>%
    add_argument("--remove_artifacts", help = "Remove common artifacts")

parsed_args <- parse_args(parser)
output_directory <- parsed_args$output_directory
aggregated_chip_calls <- parsed_args$chip_calls
aggregated_pileup_regions <- parsed_args$pileup_regions
white_list_file <- parsed_args$whitelist
remove_artifacts <- parsed_args$remove_artifacts

chip_variants <- vroom(aggregated_chip_calls, num_threads = 1)
colnames(chip_variants)[1] <- "eid"
chip_variants$eid %<>% str_remove_all("_.*$")

chip_variants %<>% separate_wider_delim(ad, delim = ",", names = c("ad_ref", "ad_alt"))
chip_variants$ad_ref %<>% as.integer
chip_variants$ad_alt %<>% as.integer
chip_variants %<>% separate_wider_delim(sb, delim = ",", names = c("sb_ref_forward", "sb_ref_reverse", "sb_alt_forward", "sb_alt_reverse"))
chip_variants$sb_ref_forward %<>% as.integer
chip_variants$sb_ref_reverse %<>% as.integer
chip_variants$sb_alt_forward %<>% as.integer
chip_variants$sb_alt_reverse %<>% as.integer
chip_variants %<>% separate_wider_delim(f1r2, delim = ",", names = c("f1r2_forward", "f1r2_reverse"))
chip_variants$f1r2_forward %<>% as.integer
chip_variants$f1r2_reverse %<>% as.integer
chip_variants %<>% separate_wider_delim(f2r1, delim = ",", names = c("f2r1_forward", "f2r1_reverse"))
chip_variants$f2r1_forward %<>% as.integer
chip_variants$f2r1_reverse %<>% as.integer

context_extract <- str_extract_all(chip_variants$gencode_28_referencecontext, "(.)\\1+") %>% map(nchar) %>% map(max) %>% unlist
context_extract[is.infinite(context_extract)] <- 1
chip_variants$longest_repeat <- context_extract

binom_p <- map2(chip_variants$ad_alt, chip_variants$dp, binom.test) %>% map_dbl(extract2, "p.value")
#write_rds(binom_p, "binom_p.rda")
chip_variants$binom_p <- binom_p

chip_variants %<>% select(eid, sample, filter, ad_ref:fad, sb_ref_forward:sb_alt_reverse, longest_repeat, binom_p, gencode_28_hugosymbol:gencode_xrefseq_prot_acc, gq:ps, as_filterstatus:tlod)

white_list <- read_excel(white_list_file, col_names = TRUE)
white_list$start %<>% as.integer
white_list$end %<>% as.integer

chip_variants$consequence_reformat <- chip_variants$consequence
chip_variants$consequence_reformat %<>% str_replace_all("MISSENSE", "missense") %>% str_replace_all("FRAME_SHIFT_DEL", "frameshift") %>% 
                                    str_replace_all("FRAME_SHIFT_INS", "frameshift") %>% 
                                    str_replace_all("IN_FRAME_DEL", "missense") %>%
                                    str_replace_all("IN_FRAME_INS", "missense") %>%
                                    str_replace_all("SPLICE_SITE", "splice-site") %>% 
                                    str_replace_all("NONSENSE", "nonsense")

chip_variants_singlep <- filter(chip_variants, !str_detect(gencode_28_proteinchange, ">") & !str_detect(gencode_28_proteinchange, "ins") & gencode_28_proteinchange != ".") 
chip_variants_singlep$start_position <- str_extract(chip_variants_singlep$gencode_28_proteinchange, "[0-9]+") %>% as.integer
chip_variants_singlep$end_position <- chip_variants_singlep$start_position

chip_variants_multiplep <- filter(chip_variants, str_detect(gencode_28_proteinchange, ">"))
chip_variants_multiplep_position_df <- str_extract(chip_variants_multiplep$gencode_28_proteinchange, "[0-9]+_[0-9]+") %>% str_split_fixed("_", 2)
chip_variants_multiplep$start_position <- chip_variants_multiplep_position_df[,1] %>% as.integer
chip_variants_multiplep$end_position <- chip_variants_multiplep_position_df[,2] %>% as.integer

chip_variants_insp <- filter(chip_variants, str_detect(gencode_28_proteinchange, "ins"))
chip_variants_insp_position_df <- str_extract(chip_variants_insp$gencode_28_proteinchange, "[0-9]+_[0-9]+") %>% str_split_fixed("_", 2)
chip_variants_insp$start_position <- chip_variants_insp_position_df[,1] %>% as.integer
chip_variants_insp$end_position <- chip_variants_insp_position_df[,2] %>% as.integer

chip_variants_nop <- filter(chip_variants, str_detect(gencode_28_proteinchange, "^\\.$"))
chip_variants_nop$start_position <- NA
chip_variants_nop$end_position <- NA

# Change this so that it checks that all dfs have more than 0 rows
# Plan is to write a recursive function that takes list of objects and only calls bind_rows if the next object has more than 0 rows
chip_variants_splitp <- bind_rows(
    chip_variants_singlep,
    chip_variants_multiplep,
    chip_variants_insp,
    chip_variants_nop
)

rm(chip_variants)
rm(chip_variants_singlep)
rm(chip_variants_multiplep)
rm(chip_variants_insp)
rm(chip_variants_nop)
gc()

chip_variants_whitelist <- group_by(white_list, variant_type) %>% group_modify(whitelist_annotate, chip_variants_splitp)
chip_variants_unfiltered_file <- str_c(output_directory, "chip_variants_unfiltered.rda", sep = "/")
write_rds(chip_variants_whitelist, chip_variants_unfiltered_file)

# Preliminary filtering
chip_variants_filtered1 <- filter(chip_variants_whitelist,
                                  dp >= 20 & 
                                  ad_ref >= 5 & 
                                  ad_alt >= 5 & 
                                  af >= 0.02)

chip_variants_filtered2 <- filter(chip_variants_filtered1, !str_detect(filter, "orientation|haplotype|base_qual|contamination|strand_bias|map_qual|position|germline|clustered_events"))

skipped_genes <- c("SF3A1", "GATA1", "GATA3", "PTEN", "SF1", "STAG1", "IKZF2", "IKZF3", "PDSS2", "LUC7L2", "JAK1", "JAK3", "GNA13", "KMT2A", "KMT2D", "CSF1R")

chip_variants_filtered3 <- filter(chip_variants_filtered2, !(is_in(gencode_28_hugosymbol, skipped_genes)))
chip_variants_filtered4 <- filter(chip_variants_filtered3, !(filter == "slippage" & af < 0.1))

if (remove_artifacts == "true") {
    chip_variants_count <- group_by(chip_variants_filtered4, gencode_28_hugosymbol, gencode_28_cdnachange, gencode_28_proteinchange) %>% summarise(n = n())
    chip_variants_r882h_count <- filter(chip_variants_count, gencode_28_hugosymbol == "DNMT3A" & gencode_28_proteinchange == "p.R882H")
    chip_variants_artifact <- filter(chip_variants_count, n > chip_variants_r882h_count$n & !(gencode_28_hugosymbol == "ASXL1" & gencode_28_cdnachange == "c.1926_1927insG")) 
    chip_variants_artifact_file <- str_c(output_directory, "artifacts.tsv", sep = "/")
    write_tsv(chip_variants_artifact, chip_variants_artifact_file)
    chip_variants_filtered_final <- filter(chip_variants_filtered4, !is_in(gencode_28_cdnachange, chip_variants_artifact$gencode_28_cdnachange))
} else {
    chip_variants_filtered_final <- chip_variants_filtered4
}
chip_variants_filtered_file <- str_c(output_directory, "chip_variants_filtered.rds", sep = "/")
write_rds(chip_variants_filtered_final, chip_variants_filtered_file)

#u2af1_pileups <- vroom(aggregated_pileup_regions)
#u2af1_pileups_agg <- group_by(u2af1_pileups, eid, sample, protein_change) %>% summarise(ref_count = sum(ref_count), alt_count = sum(alt_count))
#u2af1_pileups_filter <- filter(u2af1_pileups_agg, alt_count > 4)
#u2af1_pileups_filter_file <- str_c(output_directory, "u2af1_pileups_filtered.tsv", sep = "/")
#write_tsv(u2af1_pileups_filter, u2af1_pileups_filter_file)
