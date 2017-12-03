#' A compiled dataset of MHC class I epitopes and ligands.
#' 
#' \code{EpitopeDataset} is a dataframe containing 6159 epitopes and 13374 MHC ligands. Entries were retrieved from various sources. The columns Cluster60 and Cluster80 contain the clustering results by the IEDB peptide clustering tool, with threshold of 60% and 80%, respectively.
#' 
#' @format A tibble with 19533 rows and 7 variables
#' @docType data
#' @keywords datasets
#' @name Dataset_EpitopeDataset
#' @rdname Dataset_EpitopeDataset
#' @usage EpitopeDataset
#' @examples
#' dplyr::filter(EpitopeDataset, Dataset=="Calis1")
#' dplyr::filter(EpitopeDataset, stringr::str_detect(MHC, stringr::fixed("A*02"))|stringr::str_detect(MHC, "A02")|stringr::str_detect(MHC, "A2$")|stringr::str_detect(MHC, stringr::fixed("A2|")))
"EpitopeDataset"