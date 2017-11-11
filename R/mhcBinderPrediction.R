#' Prediction of MHC binders using the stabilized matrix method (SMM) approach
#' 
#' returns peptides that bind to at least one of the 50 HLA alleles implemented in the package \code{EpitopePrediction}.
#' 
#' @param peptideSet A set of peptide sequences.
#' @param ic50.threshold Peptides with predicted IC50 values lower than this will be considered binders. A threshold of 500 nM has been used frequently, according to Sette et al. https://www.ncbi.nlm.nih.gov/pubmed/7527444
#' @importFrom dplyr %>%
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr bind_rows
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom magrittr set_colnames
#' @importFrom parallel detectCores
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterEvalQ
#' @importFrom parallel clusterExport
#' @importFrom parallel stopCluster
#' @importFrom pbapply pbapply
#' @importFrom EpitopePrediction supportedMHCs
#' @importFrom EpitopePrediction smm
#' @export
#' @rdname mhcBinderPrediction
#' @name mhcBinderPrediction
mhcBinderPrediction <- function(peptideSet, ic50.threshold=500){
  cl <- parallel::makeCluster(parallel::detectCores(), type="SOCK")
  invisible(parallel::clusterEvalQ(cl, {
    library(EpitopePrediction)
  }))
  
  peptideSet <- sequenceFilter(peptideSet)
  peptideSet.8 <- peptideSet[nchar(peptideSet)==8]
  peptideSet.9 <- peptideSet[nchar(peptideSet)==9]
  peptideSet.10 <- peptideSet[nchar(peptideSet)==10]
  peptideSet.11 <- peptideSet[nchar(peptideSet)==11]
  sink(tempfile())
  mhcDF <- EpitopePrediction::supportedMHCs() %>% dplyr::filter(stringr::str_detect(mhc, "HLA-"))
  sink()
  pMHCparameterSet <- suppressWarnings(dplyr::bind_rows(
    expand.grid(peptideSet.8, dplyr::filter(mhcDF, l==8)[["mhc"]]),
    expand.grid(peptideSet.9, dplyr::filter(mhcDF, l==9)[["mhc"]]),
    expand.grid(peptideSet.10, dplyr::filter(mhcDF, l==10)[["mhc"]]),
    expand.grid(peptideSet.11, dplyr::filter(mhcDF, l==11)[["mhc"]])
  )) %>% magrittr::set_colnames(c("Peptide", "MHC"))
  pred.ic50 <- pbapply::pbapply(pMHCparameterSet, 1, function(v){EpitopePrediction::smm(v[[1]], v[[2]])}, cl=cl)
  pred.bind <- ifelse(pred.ic50<ic50.threshold, 1, 0)
  pMHCDF <- data.frame(pMHCparameterSet, "PredIC50"=pred.ic50, "PredBinding"=pred.bind)
  pMHCDF.summary <- pMHCDF %>% 
    dplyr::group_by(Peptide) %>% 
    dplyr::summarise(PredBinding=ifelse(mean(PredBinding)>0, 1, 0))
  predictedBinderSet <- dplyr::filter(pMHCDF.summary, PredBinding==1)[["Peptide"]]
  
  parallel::stopCluster(cl)
  return(predictedBinderSet)
}