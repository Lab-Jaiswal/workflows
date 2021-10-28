#BiocManager::install("amplican")
library(amplican)
library(stringr)
library(magrittr)

command_args <- commandArgs(trailingOnly = TRUE)
FASTQ_directory <- command_args[1]
output_directory <- command_args[2]
specify_config <- command_args[3]

if (specify_config == "NA"){   
    if (!str_detect(FASTQ_directory, "/$")){
        config <- paste0(FASTQ_directory, "/", "/config.csv")
    } else {
        config <- paste0(FASTQ_directory, "config.csv", sep="")
    } 
} 

if (specify_config != "NA"){ 
    if (!str_detect(FASTQ_directory, "/$") && (!str_detect(specify_config, "^/"))){
        config <- paste0(FASTQ_directory, "/", specify_config)
    } else if (str_detect(FASTQ_directory, "/$") && (str_detect(specify_config, "^/"))){
        specify_config <- str_remove(specify_config, "/")
        config <- paste0(FASTQ_directory, specify_config)
    } else {
        config <- paste0(FASTQ_directory, specify_config)
    }
}




# path to example fastq files
fastq_folder <- FASTQ_directory
# output folder, a full path
results_folder <- output_directory

#  run amplican
amplicanPipeline(config, fastq_folder, results_folder)

# results of the analysis can be found at
message(results_folder)
