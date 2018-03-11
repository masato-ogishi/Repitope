#' Generate sequence fragment libraries.
#'
#' @param sequenceSet A set of amino acid sequences.
#' @param windowSizeSet A set of sliding window sizes.
#' @importFrom Biostrings reverse
#' @importFrom stringr str_sub
#' @importFrom dplyr %>%
#' @importFrom dplyr full_join
#' @importFrom dplyr mutate
#' @importFrom dplyr transmute
#' @importFrom data.table rbindlist
#' @importFrom parallel detectCores
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterExport
#' @importFrom parallel stopCluster
#' @importFrom pbapply pblapply
#' @export
#' @rdname ContactPotentialProfiling_FragmentLibrary
#' @name ContactPotentialProfiling_FragmentLibrary
CPP_FragmentLibrary <- function(sequenceSet, windowSizeSet){
  fragmentLibrary_default <- function(sequenceSet, fragLen){
    cat("Fragment length = ", fragLen, "\n", sep="")
    maxLen <- max(nchar(sequenceSet), na.rm=T)
    cl <- parallel::makeCluster(parallel::detectCores()-1)
    parallel::clusterExport(cl, varlist=c("sequenceSet", "fragLen", "maxLen"), envir=environment())
    f <- unlist(pbapply::pblapply(seq(1, maxLen-fragLen+1),
                                  function(i){stringr::str_sub(sequenceSet, i, i+fragLen-1)}, cl=cl))
    f <- f[nchar(f)==fragLen]
    parallel::stopCluster(cl)
    df.fr <- as.data.frame(table(f))
    colnames(df.fr) <- c("Fragment", "Count")
    df.fr$Fragment <- as.character(df.fr$Fragment)
    df.fr <- dplyr::full_join(df.fr, dplyr::mutate(df.fr, Fragment=Biostrings::reverse(Fragment)), by="Fragment")
    df.fr <- dplyr::transmute(df.fr, Fragment, Count=rowSums(df.fr[2:3], na.rm=T), FragmentLength=fragLen)
    frCt <- sum(df.fr$Count)
    df.fr <- dplyr::transmute(df.fr, Fragment, FragmentLength, Count, Freq=Count/frCt)
    return(df.fr)
  }
  fragLibDF <- data.table::rbindlist(lapply(windowSizeSet, function(l){fragmentLibrary_default(sequenceSet, l)}))
  return(fragLibDF)
}

