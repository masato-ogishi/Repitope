#' Utility functions for sequence analysis.
#'
#' \code{sequenceFilter} filters amino acid sequences so that those containing non-standard letters are excluded.\cr
#' \code{sequenceSlidingWindow} splits the input sequences in a sliding window basis with a fixed window size.\cr
#' \code{InSilicoMutagenesis} generates single-aa-substituted sequences. Currently insertions and deletions are not supported.\cr
#'
#' @param sequenceSet A set of amino acid sequences.
#' @param windowSize A size of the sliding window.
#' @param coreN The number of cores to be used for parallelization.
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
  s <- s[nchar(s)>=8] ## minimum length
  return(s)
}

#' @export
#' @rdname Utility_Sequence
#' @name Utility_Sequence
sequenceSlidingWindow <- function(sequenceSet, windowSize=3){
  if(min(nchar(sequenceSet))<windowSize){
    warning("The windowSize parameter exceeds the minimum length of the input sequence! The parameter was adjusted.")
    windowSize <- min(nchar(sequenceSet))
  }
  f <- sapply(1:(max(nchar(sequenceSet))-windowSize+1),
              function(i){stringr::str_sub(sequenceSet, i, i+windowSize-1)})
  f <- f[nchar(f)==windowSize]
  return(f)
}

#' @export
#' @rdname Utility_Sequence
#' @name Utility_Sequence
InSilicoMutagenesis <- function(sequenceSet, coreN=parallel::detectCores(logical=F)){
  cl <- parallel::makeCluster(coreN, type="PSOCK")
  doSNOW::registerDoSNOW(cl)
  sequenceList <- split(sequenceSet, nchar(sequenceSet))
  cat("In silico substitutions...\n")
  sequenceSet_sub <- foreach::foreach(i=1:length(sequenceList), .inorder=F)%dopar%{
    l <- nchar(sequenceList[[i]][1])
    apply(data.table::CJ(1:l, Biostrings::AA_STANDARD), 1,
          function(v){
            s <- sequenceList[[i]]
            substr(s, v[1], v[1]) <- v[2]
            return(s)
          })
  }
  cat("In silico insertions...\n")
  sequenceSet_ins <- foreach::foreach(i=1:(length(sequenceList)-1), .inorder=F)%dopar%{
    l <- nchar(sequenceList[[i]][1])
    apply(data.table::CJ(1:l, Biostrings::AA_STANDARD), 1,
          function(v){
            s <- sequenceList[[i]]
            stringi::stri_sub(s, v[1], as.numeric(v[1])-1) <- v[2]
            return(s)
          })
  }
  cat("In silico deletions...\n")
  sequenceSet_del <- foreach::foreach(i=2:length(sequenceList), .inorder=F)%dopar%{
    l <- nchar(sequenceList[[i]][1])
    lapply(1:l, function(p){
      s <- sequenceList[[i]]
      stringi::stri_sub(s, p, p) <- ""
      return(s)
    })
  }
  parallel::stopCluster(cl)
  gc();gc()
  cat("Merging...\n")
  mut <- list(sequenceSet, sequenceSet_sub, sequenceSet_ins, sequenceSet_del)
  mut <- sort(fastUnique(unlist(mut)))
  return(mut)
}
