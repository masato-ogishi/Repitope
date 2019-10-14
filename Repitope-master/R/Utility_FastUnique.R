#' The fastest "unique" function.
#'
#' @param x Any vector.
#' @export
#' @rdname Utility_FastUnique
#' @name Utility_FastUnique
fastUnique <- function(x){
  dt <- data.table::data.table("X"=as.vector(x))
  data.table::setkeyv(dt, "X")
  return(dt[, logical(1), keyby=X]$X)
}


