#' A tools4ukbb function
#'
#' Function output: a column of name "Age_at_first_{disease_name}_diagnosis", column of eid values
#' "Date_of_first_{disease_name}_diagnosis" contains the date an individual received the diagnosis (as indicated by the icd_code(s) of interest)
#' 
#' @param icd_list a list of the icd10 codes you wish to investigate
#' @param dataframe the originial phenotype dataframe containing all individuals in the ukbiobank (502462 x 18158 as of 09/07/2021)
#' @param disease_name a string containing the name of the disease(s) of interest
#' @keywords age
#' @export
#' @examples
#' dx_age()

dx_age<- function(icd_list, dataframe, disease_name){
  indiv_with_disease <- individuals_with_disease(icd_list, dataframe)
  dx_positive<- filter(dataframe, is_in(eid, indiv_with_disease))
  
  first_dx_date<-dx_date(icd_list, dataframe, "X")
  first_dx_date$diagnosis_date<- ymd(first_dx_date$Date_of_first_X_diagnosis) 
  DOB_col<- date_of_birth(dx_positive)
  
  
  first_dx_date_with_DOB<- left_join(DOB_col, first_dx_date, by="eid")
  first_dx_date_with_DOB[, str_c("Age_at_first_", disease_name, "_dx")] <- interval(start= first_dx_date_with_DOB$DOB, end=first_dx_date_with_DOB$Date_of_first_X_diagnosis)/                  
    duration(n=1, unit="years")
  first_dx_date_with_age <- select(first_dx_date_with_DOB, -DOB, -Date_of_first_X_diagnosis)
}

