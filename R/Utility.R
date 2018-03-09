#' Miscellaneous utility functions.
#'
#' \code{sequenceSlidingWindow} does a sliding window with a fixed window size.\cr
#' \code{sequenceFilter} filters amino acid sequences so that those containing non-standard letters are excluded.\cr
#' \code{compressedToLongFormat} converts a compresed dataframe into a long-format dataframe. The compressed strings should be separated by "|".\cr
#' \code{standardizeHLAString} standardizes HLA class I strings.\cr
#' \code{convertToHLASupertypes} classifies standardized HLA class I strings into 10 major HLA supertypes. Note: this function is intended for use with the internally stored \code{EpitopeDataset} and therefore does not cover all HLAs.\cr
#'
#' @param sequenceSet A set of amino acid sequences.
#' @param windowSize A size of the sliding window.
#' @param df A dataframe like \code{Repitope::EpitopeDataset} which has compressed columns.
#' @param compressedColumnName A string indicating the names of the compressed column to be converted into a long format.
#' @param HLAStrings A character vector of HLA class I strings.
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

#' @export
#' @rdname Utility
#' @name Utility
convertToHLASupertypes <- function(HLAStrings){
  HLAStrings[which(HLAStrings %in% c("HLA-A1","HLA-A01","HLA-A0101","HLA-A01011"))] <- "A1"
  HLAStrings[which(HLAStrings %in% c("HLA-A2","HLA-A02","HLA-A0201","HLA-A0202","HLA-A0203","HLA-A0204","HLA-A0205","HLA-A0206","HLA-A0207","HLA-A0209","HLA-A0210","HLA-A0211","HLA-A0213","HLA-A0214","HLA-A0217","HLA-A0220"))] <- "A2"
  HLAStrings[which(HLAStrings %in% c("HLA-A3","HLA-A03","HLA-A0301","HLA-A0302","HLA-A3/11"))] <- "A3"
  HLAStrings[which(HLAStrings %in% c("HLA-A24","HLA-A2402","HLA-A2403","HLA-A2407"))] <- "A24"
  HLAStrings[which(HLAStrings %in% c("HLA-A26","HLA-A2601","HLA-A2602","HLA-A2603"))] <- "A26"
  HLAStrings[which(HLAStrings %in% c("HLA-B7","HLA-B07","HLA-B0701","HLA-B0702"))] <- "B7"
  HLAStrings[which(HLAStrings %in% c("HLA-B8","HLA-B08","HLA-B0801","HLA-B0802","HLA-B0804"))] <- "B8"
  HLAStrings[which(HLAStrings %in% c("HLA-B27","HLA-B2701","HLA-B2702","HLA-B2703","HLA-B2704","HLA-B2705","HLA-B2706","HLA-B2707","HLA-B2709"))] <- "B27"
  HLAStrings[which(HLAStrings %in% c("HLA-B44","HLA-B4401","HLA-B4402","HLA-B4403","HLA-B4405","HLA-B4415","HLA-B4427"))] <- "B44"
  HLAStrings[which(HLAStrings %in% c("HLA-B62"))] <- "B62"
  HLAStrings[which(stringr::str_detect(HLAStrings, "HLA-"))] <- "Others"
  HLAStrings[which(is.na(HLAStrings))] <- "Others"
  return(HLAStrings)
}
