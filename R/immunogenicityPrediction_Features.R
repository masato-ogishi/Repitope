#' Generate peptide features for immunogenicity prediction.
#' 
#' A wrapper function of \code{peptideDescriptorAnalysis} and \code{rTPCPAnalysis}.
#' 
#' @param peptideSet A set of peptide sequences.
#' @param TCRSet A set of TCR sequences.
#' @param aaIndexIDSet A set of AAIndexIDs indicating the pairwise matching score matrix to be used. A set of AAIndex-derived matrices can be retrieved by \code{AACPMatrix}. Set "all" to select all AAIndexIDs.
#' @param fragLen The length of TCR sequence fragments to be matched against the peptide set.
#' @param alignTypeSet A set of alignment-type strings directly passed to the \code{type} argument of the \code{pairwiseAlignment} function in the \code{Biostrings} package.
#' @param seed A random seed.
#' @param TCRFragDepth The number of TCR fragments to be matched. This should be kept constant throughout an analysis to keep repertoire-dependent fluctuation of descriptive statistics minimized.
#' @export
#' @rdname immunogenicityPrediction_Features
#' @name immunogenicityPrediction_Features
immunogenicityPrediction_Features <- function(
  peptideSet, TCRSet, 
  aaIndexIDSet="all", fragLen=5, alignTypeSet="global-local", 
  seed=12345,
  TCRFragDepth=10000
){
  message("Peptide descriptor analysis.")
  df_feature_pepDesc <- peptideDescriptorAnalysis(peptideSet, fragLen)
  gc();gc()
  message("rTPCP analysis.")
  df_feature_rTPCP <- rTPCPAnalysis(peptideSet, TCRSet, aaIndexIDSet, fragLen, alignTypeSet, seed, TCRFragDepth)
  gc();gc()
  dplyr::left_join(df_feature_pepDesc, df_feature_rTPCP, by="Peptide")
}