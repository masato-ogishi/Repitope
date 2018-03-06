#' Compiled datasets of human TCR beta CDR3 sequences
#'
#' \code{TCRSet} and \code{TCRSet_CD8} are the lists of 100,000 randomly selected human TCR clonotypes with defined random seeds derived from unsorted T cells or CD8+ T cells, respectively. The entire clonotypes are not included in this package.\cr
#' \code{TCRSet_Public} and \code{TCRSet_Public_CD8} are the sets of public TCR clonotypes. The criteria of publicity are being observed in at least 22 out of 218 datasets or being observed in at least 2 out of 6 datasets, respectively.\cr
#' \code{TCRSet_Central10000} and \code{TCRSet_Central10000_CD8} are the sets of central TCR clonotypes returned from \code{centralClonotypeAnalysis}. The size was set 10000.\cr
#' For details, see \code{help("TCRAnalysis_CentralClonotypes")} as well.
#'
#' @format A character vector.
#' @docType data
#' @keywords datasets
#' @name Dataset_TCRDatasets
#' @rdname Dataset_TCRDatasets
#' @examples
#' ## Calculate features using CD8+ T-cell public clonotypes.
#' FeatureDFList_CD8 <- Features(
#'   peptideSet=EpitopeDataset$"Peptide",
#'   TCRSet=TCRSet_Public_CD8,
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
"TCRSet_Public"

#' @format A character vector.
#' @docType data
#' @keywords datasets
#' @name Dataset_TCRDatasets
#' @rdname Dataset_TCRDatasets
"TCRSet_Public_CD8"

#' @format A character vector.
#' @docType data
#' @keywords datasets
#' @name Dataset_TCRDatasets
#' @rdname Dataset_TCRDatasets
"TCRSet_Central10000"

#' @format A character vector.
#' @docType data
#' @keywords datasets
#' @name Dataset_TCRDatasets
#' @rdname Dataset_TCRDatasets
"TCRSet_Central10000_CD8"
