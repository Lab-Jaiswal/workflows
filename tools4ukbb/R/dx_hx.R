#' A tools4ukbb function
#' 
#' Function output: a column of name "Hx_of_{disease_name}", column of eid values
#' "Hx_of_{disease_name}" contains 0's (indicating an absence of the diagnosis) or 1's (indicating the presence of a diagnosis)
#' 
#' @param icd_list a list of the icd10 codes you wish to investigate
#' @param dataframe the originial phenotype dataframe containing all individuals in the ukbiobank (~500,000 cols x 18,000 rows as of 09/07/2021)
#' @param disease_name a string containing the name of the disease(s) of interest
#' @keywords hx
#' @export
#' @examples
#' dx_hx()

dx_hx<- function(icd_list, dataframe, disease_name){
  
  indiv_with_disease <- individuals_with_disease(icd_list, dataframe)
  indiv_without_disease <- individuals_without_disease(icd_list, dataframe)
  
  dx_positive<- select(dataframe, eid) %>% filter(is_in(eid, indiv_with_disease)) 
  dx_negative<- select(dataframe, eid) %>% filter(is_in(eid, indiv_without_disease))
  
  hx_diagnosis<- str_c("Hx_of_", disease_name)
  dx_positive[[hx_diagnosis]] <- 1
  dx_negative[[hx_diagnosis]] <- 0
  
  joined_df<-rbind(dx_positive, dx_negative)
}