#' Utility functions for dataframe/datatables.
#'
#' @param df Any dataframe.
#' @param groupColumnName The column name for the grouping variables.
#' @param valueColumnName The column name for values.
#' @param descending Logical. Whether the dataframe should be reordered in a descending order.
#' @importFrom dplyr %>%
#' @importFrom dplyr select
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr arrange
#' @importFrom dplyr ungroup
#' @importFrom dplyr left_join
#' @export
#' @rdname Utility_DF
#' @name Utility_DF
sortByMedian <- function(df, groupColumnName="FeatureID", valueColumnName="Importance", descending=F){
  df.ord <- dplyr::select(df, c(groupColumnName, valueColumnName))
  colnames(df.ord) <- c("Group", "Value")
  df.ord <- df.ord %>%
    dplyr::group_by(Group) %>%
    dplyr::summarise(Value=median(Value)) %>%
    dplyr::arrange(Value) %>%
    dplyr::ungroup() %>%
    dplyr::select(Group)
  colnames(df.ord) <- groupColumnName
  res <- suppressMessages(dplyr::left_join(df.ord, df))
  if(descending==T) res <- res[nrow(res):1,]
  return(res)
}
