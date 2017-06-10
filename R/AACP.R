# Amino Acid Pairwise Contact Potential (AACP)
#' Format AAIndex Amino Acid Pairwise Contact Potential (AACP) matrix
#' @param aacp.aaindex.matrix AAIndex AACP matrix
#' @importFrom dplyr select
#' @importFrom IRanges reverse
AACP.AAIndex.Matrix.Format <- function(aacp.aaindex.matrix){
  df <- aacp.aaindex.matrix
  aaindex.names <- df[["AAIndexID"]]
  df <- dplyr::select(df, -AAIndexID, -Description, -PMID, -Author, -Title, -Reference)
  df_rev <- df
  colnames(df_rev) <- IRanges::reverse(colnames(df_rev))
  df <- data.frame(df, df_rev, check.names=F)
  df <- df[unique(colnames(df))]
  rownames(df) <- aaindex.names
  return(df)
}
