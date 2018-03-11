#' Miscellaneous utility functions.
#'
#' \code{sequenceFilter} filters amino acid sequences so that those containing non-standard letters are excluded.\cr
#' \code{sequenceSlidingWindow} does a sliding window with a fixed window size.\cr
#' \code{compressedToLongFormat} converts a compresed dataframe into a long-format dataframe. The compressed strings should be separated by "|".\cr
#' \code{hlaGenotypeToSupertype} translates HLA class I genotype strings to corresponding supertypes. Uses the genotype-supertype relation dataset from Sidney et al., 2008. \cr
#'
#' @param sequenceSet A set of amino acid sequences.
#' @param windowSize A size of the sliding window.
#' @param df A dataframe like \code{Repitope::EpitopeDataset} which has compressed columns.
#' @param compressedColumnName A string indicating the names of the compressed column to be converted into a long format.
#' @param HLAStrings A character vector of HLA class I strings.
#' @param as.data.frame Logical. Whether the result should be summarized into a dataframe of dummy variables.
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

#' @export
#' @rdname Utility
#' @name Utility
hlaGenotypeToSupertype <- function(HLAStrings, as.data.frame=T){
  # Stanadrdize HLA strings
  hla <- HLAStrings %>%
    stringr::str_replace_all("HLA-", "") %>%
    stringr::str_replace_all("HLA ", "") %>%
    stringr::str_replace_all(":", "") %>%
    stringr::str_split(stringr::fixed("|"))

  hlaGeneAB <- function(HLAStrings){
    hla.a <- stringr::str_detect(HLAStrings, "^A.+")
    hla.b <- stringr::str_detect(HLAStrings, "^B.+")
    hla <- hla.a | hla.b
    HLAStrings[which(hla==F)] <- NA
    return(HLAStrings)
  }
  hla <- pbapply::pblapply(hla, hlaGeneAB)

  removeSpace <- function(HLAStrings){
    unique(unlist(lapply(stringr::str_split(HLAStrings, " "), dplyr::first)))
  }
  hla <- pbapply::pblapply(hla, removeSpace)

  # Supertype matching
  HLA_ST <- suppressMessages(readr::read_csv(system.file("HLASupertypeTable.csv", package="Repitope")))  ## Sidney et al., 2008. Additional File 1.
  HLA_ST <- HLA_ST[1:2]
  colnames(HLA_ST) <- c("Allele", "Supertype")
  HLA_ST$Supertype[which(HLA_ST$Supertype=="Unclassified")] <- NA
  HLA_ST$Supertype <- stringr::str_replace_all(HLA_ST$Supertype, " ", "|")
  hlaSupertype <- function(HLAStrings){
    dplyr::filter(HLA_ST, Allele %in% HLAStrings)$Supertype
  }
  hla <- pbapply::pblapply(hla, hlaSupertype)
  hla <- pbapply::pblapply(hla, function(m){sort(unique(unlist(stringr::str_split(m, stringr::fixed("|")))))})
  hla[which(sapply(hla, length)==0)] <- NA
  if(as.data.frame==F) return(hla)

  # Convert to a dataframe of dummy variables
  HLA_Supertypes <- sort(unique(unlist(hla)))
  hla_df <- data.table::transpose(as.data.frame(pbapply::pblapply(hla, function(m){as.numeric(table(factor(m, levels=HLA_Supertypes)))})))
  colnames(hla_df) <- HLA_Supertypes
  return(hla_df)
}
