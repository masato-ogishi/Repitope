#' Generate sequence fragment libraries.
#'
#' \code{CPP_FragmentLibrary} generates a formatted fragment library datatable.
#'
#' @param tcrSequenceSet A set of TCR CDR3B amino acid sequences.
#' @param fragLenSet A set of sliding window sizes.
#' @param fragDepth The depth of fragments used to compute repertoire-wide TCR-peptide contact potentials.
#' @param seedSet A set of random seeds.
#' @export
#' @rdname ContactPotentialProfiling_FragmentLibrary
#' @name ContactPotentialProfiling_FragmentLibrary
CPP_FragmentLibrary <- function(tcrSequenceSet, fragLenSet=3:8, fragDepth=100000, seedSet=1:5){
  # Generate fragment libraries from the entire sequenece dataset provided
  fragmentLibrary_default <- function(tcrSequenceSet, fragLen){
    cat("Fragment length = ", fragLen, "\n", sep="")
    maxLen <- max(nchar(tcrSequenceSet), na.rm=T)
    cl <- parallel::makeCluster(parallel::detectCores(logical=F))
    parallel::clusterExport(cl, varlist=c("tcrSequenceSet", "fragLen", "maxLen"), envir=environment())
    f <- unlist(pbapply::pblapply(seq(1, maxLen-fragLen+1),
                                  function(i){stringr::str_sub(tcrSequenceSet, i, i+fragLen-1)}, cl=cl))
    f <- f[nchar(f)==fragLen]
    parallel::stopCluster(cl)
    f <- c(f, stringi::stri_reverse(f))
    dt.fr <- data.table::as.data.table(table(f))
    colnames(dt.fr) <- c("Fragment", "Count")
    dt.fr[,Freq:=Count/length(f)]
    return(dt.fr)
  }
  cat("Fragmenting...\n", sep="")
  fragLibList <- lapply(fragLenSet, function(l){fragmentLibrary_default(tcrSequenceSet, l)})
  names(fragLibList) <- fragLenSet

  # Format the libraries for contact potential profiling analysis
  format_FragLib <- function(fragLen, seed=12345){
    FragLib <- fragLibList[[as.character(fragLen)]]
    set.seed(seed); FragSet.Dedup <- sample(FragLib$Fragment, size=fragDepth, replace=T)
    set.seed(seed); FragSet.Weighted <- sample(FragLib$Fragment, size=fragDepth, replace=T, prob=FragLib$Freq)
    set.seed(seed); FragSet.Mock <- sapply(1:fragDepth, function(i){paste0(sample(Biostrings::AA_STANDARD, size=fragLen, replace=T), collapse="")})
    FragDF <- data.table::data.table("ID"=1:fragDepth, "FragLen"=fragLen, "Seed"=seed, "Deduplicated"=FragSet.Dedup, "Weighted"=FragSet.Weighted, "Mock"=FragSet.Mock)
    return(FragDF)
  }
  cat("Formatting...\n", sep="")
  paramGrid <- data.table::CJ("FragLen"=fragLenSet, "Seed"=seedSet)
  cl <- parallel::makeCluster(parallel::detectCores(logical=F))
  parallel::clusterExport(cl, varlist=c("format_FragLib", "fragLibList", "paramGrid"), envir=environment())
  fragLibDT <- pbapply::pblapply(1:nrow(paramGrid), function(i){
    format_FragLib(paramGrid$"FragLen"[i], paramGrid$"Seed"[i])
  }, cl=cl) %>%
    data.table::rbindlist()
  parallel::stopCluster(cl)
  gc();gc()
  fragLibDT <- data.table::melt.data.table(
    fragLibDT,
    id.vars=c("ID", "FragLen", "Seed"),
    measure.vars=c("Deduplicated", "Weighted", "Mock"),
    variable.name="Library", value.name="Fragment"
  )
  gc();gc()
  fragLibDT[,"Library":=paste0(Library, "_", FragLen, "_", Seed)][,FragLen:=NULL][,Seed:=NULL]
  fragLibDT <- data.table::dcast.data.table(fragLibDT, ID~Library, value.var="Fragment")
  fragLibDT[,"ID":=NULL]
  gc();gc()
  return(fragLibDT)
}

