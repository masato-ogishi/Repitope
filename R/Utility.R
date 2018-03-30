#' Miscellaneous utility functions.
#'
#' \code{sequenceFilter} filters amino acid sequences so that those containing non-standard letters are excluded.\cr
#' \code{sequenceSlidingWindow} does a sliding window with a fixed window size.\cr
#' \code{compressedToLongFormat} converts a compresed dataframe into a long-format dataframe. The compressed strings should be separated by "|".\cr
#'
#' @param sequenceSet A set of amino acid sequences.
#' @param windowSize A size of the sliding window.
#' @param df A dataframe like \code{Repitope::EpitopeDataset} which has compressed columns.
#' @param compressedColumnName A string indicating the names of the compressed column to be converted into a long format.
#' @importFrom Biostrings AA_STANDARD
#' @importFrom zoo coredata
#' @importFrom stringr str_sub
#' @importFrom stringr str_split
#' @importFrom stringr str_detect
#' @importFrom stringr fixed
#' @importFrom stringr str_replace
#' @importFrom stringr str_replace_all
#' @importFrom dplyr %>%
#' @importFrom dplyr first
#' @importFrom dplyr filter
#' @importFrom data.table rbindlist
#' @importFrom data.table transpose
#' @importFrom purrr flatten
#' @importFrom readr read_csv
#' @importFrom pbapply pblapply
#' @export
#' @rdname Utility
#' @name Utility
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
#' @rdname Utility
#' @name Utility
sequenceSlidingWindow <- function(sequenceSet, windowSize){
  f <- sapply(1:(max(nchar(sequenceSet), na.rm=T)-windowSize+1),
              function(i){stringr::str_sub(sequenceSet, i, i+windowSize-1)})
  f <- f[nchar(f)==windowSize]
  return(f)
}

#' @export
#' @rdname Utility
#' @name Utility
compressedToLongFormat <- function(df, compressedColumnName){
  cmp <- stringr::str_split(df[[compressedColumnName]], stringr::fixed("|"))
  l <- sapply(cmp, length)
  df_long <- lapply(1:nrow(df), function(i){replicate(l[i], zoo::coredata(df[i,]), simplify=F)})
  df_long <- data.table::rbindlist(purrr::flatten(df_long))
  df_long[[compressedColumnName]] <- unlist(cmp)
  return(df_long)
}
