#' @title Reading a CSV file and remove spaces from the column names
#' @description Reading a CSV file and remove spaces from the column names
#' @param fileName a CSV file name
#' @importFrom dplyr %>%
#' @importFrom stringr str_replace_all
#' @importFrom magrittr set_colnames
#' @importFrom readr read_csv
read_csv_checkColumnNames <- function(fileName){
  readr::read_csv(fileName) %>%
    (function(d){magrittr::set_colnames(d, stringr::str_replace_all(colnames(d), " ", "."))})
}
