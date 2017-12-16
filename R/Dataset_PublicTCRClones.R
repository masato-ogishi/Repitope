#' A compiled dataset of human public TCR beta CDR3 clonotypes.
#'
#' \code{TCRSet_CD8_PubClones10000} is a TCR repertoire derived from CD8+ human T cells. The size threshold was set 10000.
#' \code{TCRSet_PubClones10000} is a TCR repertoire derived from unsorted human T cells. The size threshold was set 10000.
#' For details, type \code{help("TCRAnalysis_PublicClonotypes")}.
#'
#' @format A character vector.
#' @docType data
#' @keywords datasets
#' @name Dataset_PublicTCRClones
#' @rdname Dataset_PublicTCRClones
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
"TCRSet_CD8_PubClones10000"
"TCRSet_PubClones10000"
