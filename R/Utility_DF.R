#' Utility functions for dataframes.
#'
#' @param df Any dataframe.
#' @param groupColumnName The column name for the grouping variables.
#' @param valueColumnName The column name for values.
#' @param descending Logical. Whether the dataframe should be reordered in a descending order.
#' @param .predictate A predictate function to be applied to the rows or a logical vector. Note: columns to be mutated have to be created in advance.
#' @param ... Other parameters passed to \code{dplyr::mutate}.
#' @importFrom dplyr %>%
#' @importFrom dplyr select
#' @importFrom dplyr mutate
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr arrange
#' @importFrom dplyr ungroup
#' @importFrom dplyr left_join
#' @importFrom rlang enquo
#' @importFrom rlang eval_tidy
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

#' @export
#' @rdname Utility_DF
#' @name Utility_DF
mutate_if_byrow <- function(df, .predictate, ...) {
  .predictate <- rlang::enquo(.predictate)
  .predictate_lgl <- rlang::eval_tidy(.predictate, df)
  df[.predictate_lgl, ] <- df[.predictate_lgl, ] %>% dplyr::mutate(...)
  return(df)
}
