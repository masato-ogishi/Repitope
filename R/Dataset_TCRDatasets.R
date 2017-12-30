#' Compiled datasets of human TCR beta CDR3 sequences
#'
#' \code{TCRSet} is a combined human TCR repertoire derived from unsorted T cells. TCR clonotypes observed in at least 22 out of 218 datasets were retained.
#' \code{TCRSet_CD8} is a combined human CD8+ TCR repertoire. TCR clonotypes observed in at least 2 out of 6 donors were retained.\cr
#' \code{TCRSet_PubClones10000} is a public TCR repertoire returned from \code{publicClonotypeAnalysis}.
#' \code{TCRSet_CD8_PubClones10000} is a public CD8+ TCR repertoire returned from \code{publicClonotypeAnalysis}. The size threshold was set 10000.
#' For details, see \code{help("TCRAnalysis_PublicClonotypes")} as well.
#'
#' @format A character vector.
#' @docType data
#' @keywords datasets
#' @name Dataset_TCRDatasets
#' @rdname Dataset_TCRDatasets
#' @examples
#' ## Calculate features.
#' FeatureDFList_CD8 <- Features(
#'   peptideSet=EpitopeDataset$"Peptide",
#'   TCRSet=TCRSet_CD8_PubClones10000,
#'   aaIndexIDSet="all",
#'   alignTypeSet="global-local",
#'   fragLenSet=3:8,
#'   TCRFragDepthSet=10000,
#'   seedSet=1:5
#' )
#'
"TCRSet"

#' @format A character vector.
#' @docType data
#' @keywords datasets
#' @name Dataset_TCRDatasets
#' @rdname Dataset_TCRDatasets
"TCRSet_CD8"

#' @format A character vector.
#' @docType data
#' @keywords datasets
#' @name Dataset_TCRDatasets
#' @rdname Dataset_TCRDatasets
"TCRSet_PubClones10000"

#' @format A character vector.
#' @docType data
#' @keywords datasets
#' @name Dataset_TCRDatasets
#' @rdname Dataset_TCRDatasets
"TCRSet_CD8_PubClones10000"
