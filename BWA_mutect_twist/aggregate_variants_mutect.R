library(readxl)
library(GenomicRanges)
library(magrittr)
library(tidyverse)

SetNames <- function(col_names, col_values) {
    names(col_values) <- col_names
    col_values
}

command_args <- commandArgs(trailingOnly = TRUE)
panel_coordinates <- command_args[1]
mutect_directory <- command_args[2]

sample_names <- list.files(mutect_directory, pattern = "*_funcotator.vcf$") %>% str_remove_all("_.*$")

mutect_vcf_files <- list.files(mutect_directory, pattern = "*_funcotator.vcf$", full.names = TRUE)
mutect_vcf_header <- read_lines(mutect_vcf_files[1]) 
funcotator_columns <- str_subset(mutect_vcf_header, "ID=FUNCOTATION") %>% str_remove_all("^.*: ") %>% str_remove_all('\\".*$') %>% str_split("\\|") %>% extract2(1)

vcf_colnames <- str_subset(mutect_vcf_header, "^#") %>% str_subset("^##", negate = T) %>% str_split("\\t") %>% extract2(1)
vcf_colnames[length(vcf_colnames)] <- "DATA"
mutect_vcf_list <- map(mutect_vcf_files, read_lines) %>% map(str_subset, "^#", negate = TRUE) %>% map(str_split_fixed, "\\t", length(vcf_colnames)) %>% map(set_colnames, vcf_colnames) %>% map(as_tibble) %>% set_names(sample_names)
mutect_vcf_all <- bind_rows(mutect_vcf_list, .id = "Sample")

mutect_info_names <- str_split(mutect_vcf_all$INFO, ";") %>% map(str_remove_all, "=.*$") 
mutect_info <- str_split(mutect_vcf_all$INFO, ";") %>% map(str_remove_all, "^.*\\=") 
mutect_info_df <- map2(mutect_info_names, mutect_info, SetNames) %>% bind_rows

funcotator_info <- str_remove_all(mutect_info_df$FUNCOTATION, "^\\[") %>% str_remove_all("\\]$") %>% str_split_fixed("\\|", length(funcotator_columns))  %>% as_tibble %>% set_colnames(funcotator_columns)
mutect_info_df_final <- select(mutect_info_df, -FUNCOTATION) %>%  select(-DP)

mutect_data_names <- str_split(mutect_vcf_all$FORMAT, ":")  
mutect_data <- str_split(mutect_vcf_all$DATA, ":") 
mutect_data_df <- map2(mutect_data_names, mutect_data, SetNames) %>% bind_rows %>% select(-SB)

mutect_vcf_bind <- select(mutect_vcf_all, Sample:FILTER) %>% bind_cols(mutect_data_df) %>% bind_cols(funcotator_info) %>% bind_cols(mutect_info_df_final) 
colnames(mutect_vcf_bind) %<>% str_remove_all("Gencode_28_") %>% 
    str_replace_all("chromosome", "Chromosome") %>% 
    str_replace_all("ncbiBuild", "NCBI_Build") %>%
    str_replace_all("hugoSymbol", "Hugo_Symbol") %>%
    str_replace_all("variantClassification", "Variant_Classification") %>% 
    str_replace_all("variantType", "Variant_Type") %>% 
    str_replace_all("proteinChange", "Protein_Change") %>% 
    str_replace_all("refAllele", "Reference_Allele") %>%
    str_replace_all("tumorSeqAllele1", "Tumor_Seq_Allele1") %>% 
    str_replace_all("tumorSeqAllele2", "Tumor_Seq_Allele2") %>% 
    str_replace_all("genomeChange", "Genome_Change") %>% 
    str_replace_all("annotationTranscript", "Annotation_Transcript") %>%
    str_replace_all("transcriptStrand", "Transcript_Strand") %>% 
    str_replace_all("transcriptExon", "Transcript_Exon") %>%
    str_replace_all("transcriptPos", "Transcript_Position") %>% 
    str_replace_all("cDnaChange", "cDNA_Change") %>% 
    str_replace_all("codonChange", "Codon_Change") %>%
    str_replace_all("otherTranscripts", "Other_Transcripts") %>% 
    str_replace_all("gcContent", "gc_content") %>% 
    str_replace_all("referenceContext", "ref_context") %>%
    str_replace_all("^start$", "Start_Position") %>% 
    str_replace_all("^end$", "End_Position") %>% 
    str_replace_all("^AF$", "tumor_f") %>% 
    str_replace_all("Oreganno", "OREGANNO")

mutect_vcf_filter <- filter(mutect_vcf_bind, is_in(Variant_Classification, c("FRAME_SHIFT_DEL", "FRAME_SHIFT_INS", "MISSENSE", "NONSENSE", "SPLICE_SITE"))) %>% 
    filter(nchar(Protein_Change) > 0)

mutect_vcf_ncol <- str_split(mutect_vcf_filter$AD, ",") %>% map_int(length) %>% max

Split_Columns <- function(column, nocolumns, name1, name2) {
  object<-str_split_fixed(column, ",", nocolumns)  %>%
    as_tibble %>%
    mutate(across(everything(), as.numeric))
  colnames(object) <- c(name1, str_c(name2, seq(1:(nocolumns- 1))))
  object
}
mutect_vcf_ad <-Split_Columns(mutect_vcf_filter$AD, mutect_vcf_ncol, "t_ref_count", "t_alt_count_")
mutect_vcf_af <- Split_Columns(mutect_vcf_filter$tumor_f, mutect_vcf_ncol, "tumor_f", "")
mutect_vcf_f2r1 <- Split_Columns(mutect_vcf_filter$F2R1, mutect_vcf_ncol, "f2r1_reference", "f2r1_alternate")
mutect_vcf_f1r2 <- Split_Columns(mutect_vcf_filter$F1R2, mutect_vcf_ncol, "f1r2_reference", "f1r2_alternate")
mutect_vcf_mbq <- Split_Columns(mutect_vcf_filter$MBQ, mutect_vcf_ncol, "mbq_reference", "mbq_alternate")
mutect_vcf_mmq <- Split_Columns(mutect_vcf_filter$MMQ, mutect_vcf_ncol, "mmq_reference", "mmq_alternate")
mutect_vcf_mfrl <- Split_Columns(mutect_vcf_filter$MFRL, mutect_vcf_ncol, "mfrl_reference", "mfrl_alternate")
mutect_vcf_sb <- str_split_fixed(mutect_vcf_filter$AS_SB_TABLE, "\\|", mutect_vcf_ncol) %>% as_tibble 
colnames(mutect_vcf_sb) <- c("sb_reference", str_c("sb_alt", seq(1:(mutect_vcf_ncol - 1))))

#mutect_vcf_sa_map <- str_split_fixed(mutect_vcf_filter$SA_MAP_AF, ",", 3) %>% set_colnames(c("sa_map_forward", "sa_map_reverse", "sa_map_none")) %>% as_tibble %>% mutate(across(everything(), as.numeric))
#mutect_vcf_sa_post <- str_split_fixed(mutect_vcf_filter$SA_POST_PROB, ",", 3) %>% set_colnames(c("sa_post_forward", "sa_post_reverse", "sa_post_none")) %>% as_tibble %>% mutate(across(everything(), as.numeric))

#mutect_vcf_filter2 <- select(mutect_vcf_filter, -AD, -F2R1, -F1R2, -MBQ, -MFRL, -SA_MAP_AF, -SA_POST_PROB) %>% 
    #bind_cols(mutect_vcf_ad, mutect_vcf_f1r2, mutect_vcf_f2r1, mutect_vcf_mbq, mutect_vcf_mfrl, mutect_vcf_sa_map, mutect_vcf_sa_post) %>%
    #filter(t_ref_count > 0)

mutect_vcf_filter2 <- select(mutect_vcf_filter, -AD, -tumor_f, -F2R1, -F1R2, -MBQ, -MMQ, -MFRL, -AS_SB_TABLE) %>% 
    bind_cols(mutect_vcf_ad, mutect_vcf_af, mutect_vcf_f2r1, mutect_vcf_f1r2, mutect_vcf_mbq, mutect_vcf_mfrl, mutect_vcf_mmq, mutect_vcf_sb) %>%
    filter(t_ref_count > 0)

DetectRepeat <- function(ref_context) {
    repeat_a <- str_replace_all(ref_context, "C|G|T", ";") %>% str_split(";") %>% extract2(1) %>% nchar %>% max
    repeat_c <- str_replace_all(ref_context, "A|G|T", ";") %>% str_split(";") %>% extract2(1) %>% nchar %>% max
    repeat_g <- str_replace_all(ref_context, "A|C|T", ";") %>% str_split(";") %>% extract2(1) %>% nchar %>% max
    repeat_t <- str_replace_all(ref_context, "A|C|G", ";") %>% str_split(";") %>% extract2(1) %>% nchar %>% max

    max(repeat_a, repeat_c, repeat_g, repeat_t)
}

mutect_vcf_filter2$longest_repeat <- map_int(mutect_vcf_filter2$ref_context, DetectRepeat)

#mutect_vcf_select <- select(mutect_vcf_filter2, Sample, Hugo_Symbol, NCBI_Build:Variant_Classification, Variant_Type, Protein_Change, FILTER, tumor_f, t_ref_count, contains("t_alt_count"), Reference_Allele:Tumor_Seq_Allele2, Genome_Change:Codon_Change, gc_content, longest_repeat, contains("f1r2"), contains("f2r1"), contains("mbq"), contains("mfrl"), contains("sa_map"), contains("sa_post"), DP, ECNT, POP_AF:TLOD, RPA:STR, OREGANNO_ID, OREGANNO_Values, Other_Transcripts, ref_context)
mutect_vcf_select <- select(mutect_vcf_filter2, Sample, Hugo_Symbol:Variant_Classification, Variant_Type, Protein_Change, FILTER, contains("tumor_f"), t_ref_count, contains("t_alt_count"), 
    Reference_Allele:Tumor_Seq_Allele2, Genome_Change:Codon_Change, gc_content, longest_repeat, DP, 
    contains("f1r2"), contains("f2r1"), contains("mbq"), contains("mfrl"), contains("mmq"), matches("^sb_"), 
    AS_FilterStatus:STRQ, GT, PGT:PS, OREGANNO_ID, OREGANNO_Values, Other_Transcripts, ref_context)

twist_panel <- read_excel(panel_coordinates, col_names = F) 
colnames(twist_panel) <- c("chr", "start", "end", "Transcript", "X5", "Strand", "Gene", "X8"
twist_panel_granges <- select(twist_panel, chr:end) %>% makeGRangesFromDataFrame

mutect_vcf_granges <- select(mutect_vcf_select, Chromosome:End_Position) %>% set_colnames(c("chr", "start", "end")) %>% makeGRangesFromDataFrame
mutect_overlaps <- findOverlaps(twist_panel_granges, mutect_vcf_granges) %>% as_tibble
mutect_overlaps_sorted <- sort(mutect_overlaps$subjectHits)

mutect_vcf_filter3 <- dplyr::slice(mutect_vcf_select, mutect_overlaps_sorted) %>% dplyr::filter(is_in(Hugo_Symbol, unique(twist_panel$Gene))) %>% arrange(Sample, Chromosome, Start_Position)

mutect_basic_file <- str_c(mutect_directory, "/mutect_aggregated_simple.tsv")
write_tsv(mutect_vcf_filter3, mutect_basic_file)

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

ParseHaplotypeCallerVCF <- function(haplotypecaller_vcf) {
    # This function is not very elegant - each row of the INFO, FORMAT columns can have a different number of fields 
    # We parse each row separately as a vector which we assign names to.  Then we use bind_rows to automatically figure out the right column names
    # In subject with homozygous reference genotypes, most of the columns will be empty because they have no meaning when no variant is present

    haplotypecaller_info_names <- str_split(haplotypecaller_vcf$INFO, ";") %>% map(str_remove_all, "=.*$") 
    haplotypecaller_info <- str_split(haplotypecaller_vcf$INFO, ";") %>% map(str_remove_all, "^.*\\=") 
    haplotypecaller_info_df <- map2(haplotypecaller_info_names, haplotypecaller_info, SetNames) %>% bind_rows

    haplotypecaller_info_nums <- mutate(haplotypecaller_info_df, across(everything(), as.numeric)) %>% select(-DP)

    haplotypecaller_data_names <- str_split(haplotypecaller_vcf$FORMAT, ":")  
    haplotypecaller_data <- str_split(haplotypecaller_vcf$DATA, ":") 
    haplotypecaller_data_df <- map2(haplotypecaller_data_names, haplotypecaller_data, SetNames) %>% bind_rows

    haplotypecaller_vcf_bind <- select(haplotypecaller_vcf, `#CHROM`:FILTER) %>% bind_cols(haplotypecaller_data_df) %>% bind_cols(haplotypecaller_info_nums)
    haplotypecaller_vcf_bind
}

haplotypecaller_vcf_parsed <- map(haplotypecaller_vcf_list, ParseHaplotypeCallerVCF) %>% set_names(sample_names) %>% bind_rows(.id = "Sample")
haplotypecaller_vcf_parsed$AC %<>% replace_na(0)
haplotypecaller_vcf_parsed$AF %<>% replace_na(0)
haplotypecaller_vcf_parsed$AN %<>% replace_na(1)

haplotypecaller_vcf_granges <- select(haplotypecaller_vcf_parsed, `#CHROM`:POS) %>% set_colnames(c("chr", "start")) %>% mutate(end = start) %>% makeGRangesFromDataFrame
mcols(haplotypecaller_vcf_granges) <- select(haplotypecaller_vcf_parsed, -`#CHROM`, -POS)

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

haplotypecaller_vcf_annot_filter <- arrange(haplotypecaller_vcf_annot, Sample, rsID) %>% group_by(rsID, Sample) %>% slice(1)

haplotypecaller_germline_file <- str_c(mutect_directory, "/haplotypecaller_germline_genotypes.tsv")
write_tsv(haplotypecaller_vcf_annot_filter, haplotypecaller_germline_file)

#DiffBasePair <- function(reference_allele) {
    #reference_allele_split <- str_split(reference_allele, "") %>% extract2(1)
    #unique_base_pairs <- unique(reference_allele_split) %>% length
    #if(unique_base_pairs > 1) {
        #complex_deletion = TRUE
    #} else {
        #complex_deletion = FALSE
    #}
    #complex_deletion
#}

#DeletionRepeat <- function(ref_allele, ref_context) {
    #ref_allele_unique <- str_split(ref_allele, "") %>% extract2(1) %>% unique
    #extract_regex <- str_c("^[", ref_allele_unique, "]*")
    #ref_context_trimmed <- str_sub(ref_context, 12)
    #ref_context_repeat <- str_extract(ref_context_trimmed, extract_regex)
    #nchar(ref_context_repeat)
#}

#ComplexDeletionRepeat <- function(ref_allele, ref_context) {
    #extract_regex <- str_c("^(", ref_allele, ")*")
    #ref_context_trimmed <- str_sub(ref_context, 12)
    #ref_context_repeat <- str_extract(ref_context_trimmed, extract_regex)
    #nchar(ref_context_repeat) / nchar(ref_allele)
#}

#not_deletions <- filter(mutect_vcf_filter2, Variant_Classification != "FRAME_SHIFT_DEL")
#not_deletions$complex_deletion <- NA
#not_deletions$deletion_repeat <- NA

#snv_deletions <- filter(mutect_vcf_filter2, Variant_Classification == "FRAME_SHIFT_DEL" & nchar(Tumor_Seq_Allele1) == 1)
#snv_deletions$complex_deletion <- FALSE

#long_deletions <- filter(mutect_vcf_filter2, Variant_Classification == "FRAME_SHIFT_DEL" & nchar(Tumor_Seq_Allele1) > 1)
#long_deletions$complex_deletion <- map_lgl(long_deletions$Tumor_Seq_Allele1, DiffBasePair)

#deletions <- bind_rows(snv_deletions, long_deletions)
#simple_deletions <- filter(deletions, complex_deletion == FALSE)

#simple_deletions$deletion_repeat <- map2_int(simple_deletions$Tumor_Seq_Allele1, simple_deletions$ref_context, DeletionRepeat)

#complex_deletions <- filter(deletions, complex_deletion == TRUE)
#complex_deletions$deletion_repeat <- map2_dbl(complex_deletions$Tumor_Seq_Allele1, complex_deletions$ref_context, ComplexDeletionRepeat)

#deletions_final <- bind_rows(simple_deletions, complex_deletions)
#mutect_vcf_filter3 <- bind_rows(not_deletions, deletions_final)

