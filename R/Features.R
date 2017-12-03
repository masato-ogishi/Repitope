#' Generate peptide features for immunogenicity prediction.
#' 
#' A wrapper of \code{Features_PeptideDescriptor} and \code{Features_rTPCP}, generating a list of dataframes of peptide features calculated with different parameter sets.
#' 
#' @param peptideSet A set of peptide sequences.
#' @param TCRSet Either a set of TCR sequences (as a character vector) or a list of sets of TCR sequences. If provided as a list, it must be the same length with the seedSet.
#' @param aaIndexIDSet A set of AAIndexIDs indicating the pairwise matching score matrix to be used. A set of AAIndex-derived matrices can be retrieved by \code{AACPMatrix}. Set "all" to select all AAIndexIDs.
#' @param fragLenSet A set of the lengths of TCR sequence fragments to be matched against the peptide set.
#' @param alignTypeSet A set of alignment-type strings directly passed to the \code{type} argument of the \code{pairwiseAlignment} function in the \code{Biostrings} package.
#' @param seedSet A set of random seeds.
#' @param TCRFragDepthSet A set of the numbers of TCR fragments to be matched.
#' @export
#' @rdname Features
#' @name Features
Features <- function(
  peptideSet, TCRSet, 
  aaIndexIDSet="all", alignTypeSet="global-local", 
  fragLenSet=4:6, TCRFragDepthSet=10000, seedSet=1:5
){
  # Feature calculation
  message("Peptide descriptor analysis.")
  df_feature_peptDesc <- Features_PeptideDescriptor(peptideSet, fragLenSet) ## A single dataframe
  gc();gc()
  message("rTPCP analysis.")
  dt_feature_rTPCP <- Features_rTPCP(peptideSet, TCRSet, aaIndexIDSet, alignTypeSet, fragLenSet, TCRFragDepthSet, seedSet) ## A list of datatables
  gc();gc()
  
  # Output
  df_feature_list <- lapply(dt_feature_rTPCP, function(dt){cbind(df_feature_peptDesc, dt[,"Peptide":=NULL])})
  return(df_feature_list)
}