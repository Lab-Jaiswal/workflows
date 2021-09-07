#' A tools4ukbb function
#'
#' Function output: a column of name "Date_of_first_{disease_name}_diagnosis", column of eid values
#' "Date_of_first_{disease_name}_diagnosis" contains the date an individual received the diagnosis (as indicated by the icd_code(s) of interest)
#' 
#' @param icd_list a list of the icd10 codes you wish to investigate
#' @param dataframe the originial phenotype dataframe containing all individuals in the ukbiobank (502462 x 18158 as of 09/07/2021)
#' @param disease_name a string containing the name of the disease(s) of interest
#' @keywords date
#' @export
#' @examples
#' dx_date()

dx_date<- function(icd_list, dataframe, disease_name){
  diagnoses_1<-select(dataframe, c(eid, contains("diagnoses_icd10")))
  diagnoses_2<-select(dataframe,c(contains("diagnoses_secondary_icd10")))
  diagnoses_df<-cbind(diagnoses_1, diagnoses_2)
  
  indiv_with_disease <- individuals_with_disease(icd_list, dataframe)
  dx_positive<- filter(dataframe, is_in(eid, indiv_with_disease))
  icd_table <- filter(diagnoses_df, is_in(eid, indiv_with_disease))
  
  diagnosis_long <- pivot_longer(icd_table, -eid, names_to = "Diagnosis_Column", values_to = "diagnosis")
  diagnosis_long$visit_num <- str_remove_all(diagnosis_long$Diagnosis_Column, "diagnoses_icd10_f41270_")
  diagnosis_long_filter <- filter(diagnosis_long, is_in(diagnosis, icd_list))
  
  date_dataframe<-select(dx_positive, c(eid, contains("date_of_first_inpatient_diagnosis_icd10"))) %>% lapply(as.character ) %>% as_tibble()
  date_long <- pivot_longer(date_dataframe, -eid, names_to = "Date_Column", values_to = "date_of_diagnosis")  
  date_long$visit_num <- str_remove_all(date_long$Date_Column, "date_of_first_inpatient_diagnosis_icd10_f41280_")
  
  date_long[, 1]<- sapply(date_long[, 1], as.numeric)
  diagnosis_date <- left_join(diagnosis_long_filter, date_long)
  
  get_first__date <- function(diagnosis_date_rows, diagnosis_date_group) {
    first_diag <- arrange(diagnosis_date_rows, date_of_diagnosis) %>% slice(1)
    select(first_diag, date_of_diagnosis)
  }
  diagnosis_first_date <- group_by(diagnosis_date, eid) %>% group_modify(get_first__date)
  names(diagnosis_first_date)[names(diagnosis_first_date) == "date_of_diagnosis"] <- str_c("Date_of_first_", disease_name, "_diagnosis")
  diagnosis_first_date <- diagnosis_first_date
}