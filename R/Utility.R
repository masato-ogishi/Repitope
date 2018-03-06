#' Miscellaneous utility functions.
#'
#' \code{sequenceSlidingWindow} does a sliding window with a fixed window size.
#' \code{sequenceFilter} filters amino acid sequences so that those containing non-standard letters are excluded.
#' \code{compressedToLongFormat} converts a compresed dataframe into a long-format dataframe. The compressed strings should be separated by "|".
#' \code{standardizeHLAString} standardizes HLA class I strings.
#'
#' @param sequenceSet A set of amino acid sequences.
#' @param windowSize A size of the sliding window.
#' @param df A dataframe like \code{Repitope::EpitopeDataset} which has compressed columns.
#' @param compressedColumnName A string indicating the names of the compressed column to be converted into a long format.
#' @param HLAStrings A character vector of HLA class I strings to be standardized.
#' @importFrom Biostrings AA_STANDARD
#' @importFrom stringr str_split
#' @importFrom stringr str_detect
#' @importFrom stringr fixed
#' @importFrom stringr str_replace
#' @importFrom zoo coredata
#' @importFrom data.table rbindlist
#' @importFrom purrr flatten
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
compressedToLongFormat <- function(df, compressedColumnName){
  cmp <- stringr::str_split(df[[compressedColumnName]], stringr::fixed("|"))
  l <- sapply(cmp, length)
  df_long <- lapply(1:nrow(df), function(i){replicate(l[i], zoo::coredata(df[i,]), simplify=F)})
  df_long <- data.table::rbindlist(purrr::flatten(df_long))
  df_long[[compressedColumnName]] <- unlist(cmp)
  return(df_long)
}

#' @export
#' @rdname Utility
#' @name Utility
standardizeHLAString <- function(HLAStrings){
  HLAStrings <- as.character(HLAStrings)
  HLAStrings <- stringr::str_replace(paste0("HLA-", HLAStrings), "HLA-HLA", "HLA")
  HLAStrings[grep("^HLA-$", HLAStrings, value=F)] <- NA
  HLAStrings[grep("^HLA-NA$", HLAStrings, value=F)] <- NA
  HLAStrings[grep("HLA-?", HLAStrings, value=F, fixed=T)] <- NA
  HLAStrings[grep("^HLA-allele undetermined$", HLAStrings, value=F)] <- NA
  HLAStrings[grep("^HLA-Class I$", HLAStrings, value=F)] <- NA
  HLAStrings[grep("^HLA class I$", HLAStrings, value=F)] <- NA
  HLAStrings <- stringr::str_split(HLAStrings, " ", simplify=T)[,1]
  HLAStrings <- stringr::str_replace(HLAStrings, stringr::fixed("*"), "")
  HLAStrings <- stringr::str_replace(HLAStrings, stringr::fixed(":"), "")
  return(HLAStrings)
}
