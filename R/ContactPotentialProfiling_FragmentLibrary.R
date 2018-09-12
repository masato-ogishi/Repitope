#' Generate sequence fragment libraries.
#'
#' \code{CPP_FragmentLibrary} generates a formatted fragment library datatable.
#'
#' @param sequenceSet A set of amino acid sequences.
#' @param fragLenSet A set of sliding window sizes. Must be between 3 and 8.
#' @param maxFragDepth The maximum size of the fragment library.
#' @param seedSet A set of random seeds.
#' @export
#' @rdname ContactPotentialProfiling_FragmentLibrary
#' @name ContactPotentialProfiling_FragmentLibrary
CPP_FragmentLibrary <- function(sequenceSet, fragLenSet=3:8, maxFragDepth=100000, seedSet=1:5){
  # Generate fragment libraries from the entire sequenece dataset provided
  fragmentLibrary_default <- function(sequenceSet, fragLen){
    cat("Fragment length = ", fragLen, "\n", sep="")
    maxLen <- max(nchar(sequenceSet), na.rm=T)
    cl <- parallel::makeCluster(parallel::detectCores(logical=F))
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
    df.fr <- dplyr::transmute(df.fr, Fragment, Count, Freq=Count/frCt)
    return(df.fr)
  }
  cat("Fragmenting...\n", sep="")
  fragLibList <- lapply(fragLenSet, function(l){fragmentLibrary_default(sequenceSet, l)})
  names(fragLibList) <- fragLenSet

  # Format the libraries for contact potential profiling analysis
  format_FragLib <- function(fragLen, fragDepth=NULL, seed=12345){
    FragLib <- fragLibList[[as.character(fragLen)]]
    if(is.null(fragDepth)) fragDepth <- 100000
    set.seed(seed); FragSet.Dedup <- sample(FragLib$Fragment, size=fragDepth, replace=T)
    set.seed(seed); FragSet.Weighted <- sample(FragLib$Fragment, size=fragDepth, replace=T, prob=FragLib$Freq)
    set.seed(seed); FragSet.Mock <- sapply(1:fragDepth, function(i){paste0(sample(Biostrings::AA_STANDARD, size=fragLen, replace=T), collapse="")})
    FragDF <- data.frame("ID"=1:fragDepth, "FragLen"=fragLen, "FragDepth"=fragDepth, "Seed"=seed, "Deduplicated"=FragSet.Dedup, "Weighted"=FragSet.Weighted, "Mock"=FragSet.Mock, stringsAsFactors=F)
    return(FragDF)
  }
  cat("Formatting...\n", sep="")
  paramGrid <- expand.grid(fragLenSet, seedSet)
  cl <- parallel::makeCluster(parallel::detectCores(logical=F))
  parallel::clusterExport(cl, varlist=c("format_FragLib", "fragLibList", "paramGrid"), envir=environment())
  fragLibDT <- pbapply::pblapply(1:nrow(paramGrid), function(i){
    format_FragLib(paramGrid[i, 1], maxFragDepth, paramGrid[i, 2])
  }, cl=cl) %>%
    data.table::rbindlist() %>%
    data.table::as.data.table() %>%
    data.table::melt.data.table(id.vars=c("ID", "FragLen", "FragDepth", "Seed"),
                                measure.vars=c("Deduplicated", "Weighted", "Mock"),
                                variable.name="Library", value.name="Fragment")
  fragLibDT[,"Library":=paste0(Library, "_", FragLen, "_", Seed)][,FragLen:=NULL][,FragDepth:=NULL][,Seed:=NULL]
  fragLibDT <- data.table::dcast.data.table(fragLibDT, ID~Library, value.var="Fragment")
  fragLibDT[,"ID":=NULL]
  parallel::stopCluster(cl)
  return(fragLibDT)
}

