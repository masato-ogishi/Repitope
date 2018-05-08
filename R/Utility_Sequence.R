#' Sequence utility functions.
#'
#' \code{sequenceFilter} filters amino acid sequences so that those containing non-standard letters are excluded.\cr
#' \code{sequenceSlidingWindow} does a sliding window with a fixed window size.\cr
#'
#' @param sequenceSet A set of amino acid sequences.
#' @param windowSize A size of the sliding window.
#' @export
#' @rdname Utility_Sequence
#' @name Utility_Sequence
sequenceFilter <- function(sequenceSet){
  s <- sequenceSet[!is.na(sequenceSet)]
  s <- toupper(s)
  letters <- unique(unlist(stringr::str_split(s, "")))
  letters.exclude <- setdiff(letters, Biostrings::AA_STANDARD)
  for(l in letters.exclude){
    s <- s[!stringr::str_detect(s, stringr::fixed(l))]
  }
  return(s)
}

#' @export
#' @rdname Utility_Sequence
#' @name Utility_Sequence
sequenceSlidingWindow <- function(sequenceSet, windowSize){
  f <- sapply(1:(max(nchar(sequenceSet), na.rm=T)-windowSize+1),
              function(i){stringr::str_sub(sequenceSet, i, i+windowSize-1)})
  f <- f[nchar(f)==windowSize]
  return(f)
}
