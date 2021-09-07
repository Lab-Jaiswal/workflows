#' A tools4ukbb function
#' Function output: a list of all individuals (identified by eid) that have been diagnosed with the icd_codes of interest
#' 
#' @param icd_list a list of the icd10 codes you wish to investigate
#' @param dataframe the originial phenotype dataframe containing all individuals in the ukbiobank (~500,000 cols x 18,000 rows as of 09/07/2021)
#' @keywords with
#' @export
#' @examples
#' individuals_with_disease()

individuals_with_disease <- function(icd_list, dataframe) {
  
  diagnoses_1<-select(dataframe, c(eid, contains("diagnoses_icd10")))
  diagnoses_2<-select(dataframe,c(contains("diagnoses_secondary_icd10")))
  diagnoses_df<-cbind(diagnoses_1, diagnoses_2)
  
  icd_dataframe <- diagnoses_df %>% filter_all(any_vars(. %in% icd_list)) %>% as_tibble()
  indiv_with_disease <- indiv_with_disease$eid
  
  }