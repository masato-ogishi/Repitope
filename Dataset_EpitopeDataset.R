#' A compiled dataset of MHC class I epitopes and ligands.
#'
#' \code{EpitopeDataset} is a dataframe containing 7106 epitopes and 14662 MHC ligands. Entries were retrieved from various sources. The columns Cluster60 and Cluster80 contain the clustering results by the IEDB peptide clustering tool, with threshold of 60% and 80%, respectively.
#'
#' @format A tibble with 21768 rows and 9 variables
#' @docType data
#' @keywords datasets
#' @name Dataset_EpitopeDataset
#' @rdname Dataset_EpitopeDataset
#' @examples
#' ## Data derived from the "Calis1" dataset.
#' dplyr::filter(EpitopeDataset, stringr::str_detect(Dataset, "Calis1"))
#'
#' ## Data restricted on the HLA-A2 molecule.
#' dplyr::filter(EpitopeDataset, stringr::str_detect(MHC, stringr::fixed("A*02"))|stringr::str_detect(MHC, "A02")|stringr::str_detect(MHC, "A2$")|stringr::str_detect(MHC, stringr::fixed("A2|")))
#'
"EpitopeDataset"
