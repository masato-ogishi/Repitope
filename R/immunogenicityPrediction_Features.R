#' Generate peptide features for immunogenicity prediction.
#' 
#' A wrapper function of \code{peptideDescriptorAnalysis} and \code{rTPCPAnalysis}.
#' 
#' @param peptideSet A set of peptide sequences.
#' @param TCRSet A set of TCR sequences.
#' @param aaIndexIDSet A set of AAIndexIDs indicating the pairwise matching score matrix to be used. A set of AAIndex-derived matrices can be retrieved by \code{AACPMatrix}. Set "all" to select all AAIndexIDs.
#' @param fragLenSet A set of the lengths of TCR sequence fragments to be matched against the peptide set.
#' @param alignTypeSet A set of alignment-type strings directly passed to the \code{type} argument of the \code{pairwiseAlignment} function in the \code{Biostrings} package.
#' @param seedSet A set of random seeds.
#' @param TCRFragDepthSet A set of the numbers of TCR fragments to be matched. This should be kept constant for comparison.
#' @export
#' @rdname immunogenicityPrediction_Features
#' @name immunogenicityPrediction_Features
immunogenicityPrediction_Features <- function(
  peptideSet, TCRSet, 
  aaIndexIDSet="all", alignTypeSet="global-local", 
  fragLenSet=5, TCRFragDepthSet=10000, seedSet=1:5
){
  message("Peptide descriptor analysis.")
  df_feature_pepDesc <- peptideDescriptorAnalysis(peptideSet, fragLenSet)
  gc();gc()
  message("rTPCP analysis.")
  df_feature_rTPCP <- rTPCPAnalysis(peptideSet, TCRSet, aaIndexIDSet, alignTypeSet, fragLenSet, TCRFragDepthSet, seedSet)
  gc();gc()
  list(df_feature_pepDesc, df_feature_rTPCP)
}