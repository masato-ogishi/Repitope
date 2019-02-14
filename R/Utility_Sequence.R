#' Utility functions for sequence analysis.
#'
#' \code{sequenceFilter} filters amino acid sequences so that those containing non-standard letters are excluded.\cr
#' \code{sequenceSlidingWindow} does a sliding window with a fixed window size.\cr
#' \code{InSilicoMutagenesis} generates single-aa-substituted sequences. Currently insertions and deletions are not supported.\cr
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

#' @export
#' @rdname Utility_Sequence
#' @name Utility_Sequence
InSilicoMutagenesis <- function(sequenceSet){
  coreN <- parallel::detectCores(logical=F)
  cl <- parallel::makeCluster(coreN, type="PSOCK")
  mut <- pbapply::pblapply(
    sequenceSet,
    function(s){
      apply(expand.grid(1:nchar(s), Biostrings::AA_STANDARD), 1, function(v){substr(s, v[1], v[1]) <- v[2]; return(s)})
    },
    cl=cl
  )
  parallel::stopCluster(cl)
  gc();gc()
  mut <- unique(unlist(mut))
  return(mut)
}
