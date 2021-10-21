library(magrittr)
library(tidyverse)
library(dplyr)

command_args <- commandArgs(trailingOnly = TRUE) #Read in arguments from command line
project_directory <- command_args[1] #Extract first argument into project_directory
animal <- command_args[2]

if (animal == "human") {
    gtf_names <- read_tsv("/labs/sjaiswal/genomes/mm10/gencode_vM20/gencode.vM20.primary_assembly.annotation.genekey.txt")
} 

if (animal == "mouse"){
    gtf_names <- read_tsv("/labs/sjaiswal/genomes/mm10/gencode.vM25.annotation.genekey.txt") 
}

count_files <- list.files(str_c(project_directory, "/alignment_output/"), pattern = "*ReadsPerGene*", full.names = TRUE) #Get names of all of files containing counts for individual sampes
sample_names <- map_chr(count_files, basename) %>% str_replace_all("_L[0-9].*$", "") %>% str_replace_all("_L[0-9]{3}.*$", "") %>% str_replace_all("\\-", "_") %>% str_replace_all("\\-", " ") %>% make.names #Get samples names from file names
read_files <- map(count_files, read_tsv, col_names = FALSE) %>% #Read in each tab separated file
    map(select, X1, X2) %>% #Only keep first two columns (gene id and unstranded read count)
    map(set_colnames, c("Ensembl_ID", "Counts")) %>% #Set column names of each tibble to the same names
    reduce(left_join, by = "Ensembl_ID") %>% #Merge all tibbles with left join using only Ensembl_ID
    set_colnames(c("Ensembl_ID", sample_names)) %>% #Rename all columns after Ensembl_ID with sample names
    left_join(gtf_names) %>% #Left join with gene key to get gene symbols for gene IDs
    filter(str_detect(Ensembl_ID, 'ENSG|ENSMUSG')) %>% #Remove rows that do not contain gene IDs in Ensembl_ID column 
    select(Ensembl_ID, Symbol, !!first(sample_names):!!last(sample_names)) #Reorder columns so that Ensembl_ID and Symbol are first, followed by data

write_tsv(read_files, str_c(project_directory, "/counts.tsv")) #Write final tibble to disk as tsv file

#expr_matrix <- select(read_files, -Ensembl_ID, -Symbol) %>% as.matrix #Convert tibble into matrix for DESeq2
#rownames(expr_matrix) <- read_files$Symbol #Set rownames of expr_matrix to be gene symbols

#vst_normalized <- varianceStabilizingTransformation(expr_matrix) #Use varianceStabilizingTransformation (VST) to normalize counts
#vst_normalized <- vst_normalized - min(vst_normalized) #Subtract minimum of normalized values to make 0 lowest value

#normalized_expr <- data.frame(Ensembl_ID = read_files$Ensembl_ID, Symbol = read_files$Symbol, vst_normalized, stringsAsFactors = FALSE) %>% #Put normalized expression matrix back into data frame with gene IDs and symbols
    #mutate_if(is.numeric, round, digits = 3) #Round normalized values to 3 digits to educe file size on disk
#write_tsv(normalized_expr, str_c(project_directory, "/vst_normalized_expression.txt")) #Save tibble to disk as tsv file

#Limma voom not currently used
#library(limma)
#library(edgeR)

#counts_dge <- DGEList(expr_matrix) %>%
    #calcNormFactors(method = "TMM")
#expr_voom <- voom(counts_dge)
#normalized_expr <- data.frame(Ensembl_ID = read_files$Ensembl_ID, Symbol = read_files$Symbol, expr_voom$E, stringsAsFactors = FALSE) %>%
    #mutate_if(is.numeric, round, digits = 3)

#write_tsv(normalized_expr, str_c(project_directory, "/voom_normalized_expression.txt"))

#coldata <- data.frame(Sample = sample_names)
#expr_DESeqDataset <- DESeqDataSetFromMatrix(countData = expr_matrix, design = ~1, colData = coldata, tidy = FALSE) r
