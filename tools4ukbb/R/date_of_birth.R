#' A tools4ukbb function
#' Function output: a column of name "DOB" (in yyyy-mm-dd date format), column of eid values
#' 
#' @param dataframe the originial phenotype dataframe containing all individuals in the ukbiobank (~500,000 cols x 18,000 rows as of 09/07/2021)
#' @keywords date of birth
#' @export
#' @examples
#' date_of_birth()

date_of_birth <- function(dataframe) {
  
  birth<-select(dataframe, eid, year_of_birth_f34_0_0, month_of_birth_f52_0_0)
  birth[,"birth_month"]<- as.integer(factor(birth$month_of_birth_f52_0_0, levels = month.name))
  birth[, "DOB_numeric"]<- str_c(birth$year_of_birth_f34_0_0, "-",birth$birth_month, "-01")
  birth[,"DOB"] <- ymd(birth$DOB_numeric)
  
  Date_of_birth<-select(birth, eid, DOB)

}