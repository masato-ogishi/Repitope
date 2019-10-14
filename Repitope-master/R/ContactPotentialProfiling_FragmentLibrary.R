#' Generate sequence fragment libraries.
#'
#' \code{CPP_FragmentLibrary} generates a formatted fragment library datatable.
#'
#' @param tcrSequenceSet A set of TCR CDR3B amino acid sequences.
#' @param fragLenSet A set of sliding window sizes.
#' @param maxFragDepth The maximum depth of fragments used to compute repertoire-wide TCR-peptide contact potentials.
#' @param seedSet A set of random seeds.
#' @param coreN The number of cores to be used for parallelization.
#' @param tmpDir Destination directory to save intermediate files.
#' @export
#' @rdname ContactPotentialProfiling_FragmentLibrary
#' @name ContactPotentialProfiling_FragmentLibrary
CPP_FragmentLibrary <- function(
  tcrSequenceSet,
  fragLenSet=3:8,
  maxFragDepth=100000,
  seedSet=1:5,
  coreN=parallel::detectCores(logical=F),
  tmpDir=file.path(tempdir(), "FragLibDT", format(Sys.time(), "%Y.%m.%d.%H.%M.%S"))
){
  # Temporary directory
  dir.create(tmpDir, showWarnings=F, recursive=T)

  # Generate fragment libraries from the entire sequenece dataset provided
  maxLen <- max(nchar(tcrSequenceSet), na.rm=T)
  cl <- parallel::makeCluster(coreN, type="PSOCK")
  parallel::clusterExport(cl, varlist=c("tcrSequenceSet", "fragLenSet", "maxLen"), envir=environment())
  cat("Fragmenting...\n", sep="")
  fragLibList <- pbapply::pblapply(fragLenSet, function(l){
    f <- unlist(lapply(seq(1, maxLen-l+1), function(i){stringr::str_sub(tcrSequenceSet, i, i+l-1)}))
    f <- f[nchar(f)==l]
    f <- c(f, stringi::stri_reverse(f))
    dt.fr <- data.table::as.data.table(table(f))
    colnames(dt.fr) <- c("Fragment", "Count")
    dt.fr[,Freq:=Count/length(f)]
    return(dt.fr)
  }, cl=cl)
  names(fragLibList) <- fragLenSet
  parallel::stopCluster(cl)
  gc();gc()

  # Format the libraries for contact potential profiling analysis
  cat("Formatting...\n", sep="")
  for(s in seedSet){
    cat("Random seed = ", s, "\n", sep="")
    set.seed(s)
    fragLibDT <- foreach::foreach(l=fragLenSet)%do%{
      FragLib <- fragLibList[[as.character(l)]]
      FragSet.Dedup <- sample(FragLib$Fragment, size=maxFragDepth, replace=T)
      FragSet.Weighted <- sample(FragLib$Fragment, size=maxFragDepth, replace=T, prob=FragLib$Freq)
      FragSet.Mock <- sapply(1:maxFragDepth, function(i){paste0(sample(Biostrings::AA_STANDARD, size=l, replace=T), collapse="")})
      FragDF <- data.table::data.table("ID"=1:maxFragDepth, "FragLen"=l, "Seed"=s, "Deduplicated"=FragSet.Dedup, "Weighted"=FragSet.Weighted, "Mock"=FragSet.Mock)
      return(FragDF)
    } %>%
      data.table::rbindlist()
    fst::write_fst(fragLibDT, file.path(tmpDir, paste0("FragLibDT_Seed", s, ".fst")))
    gc();gc()
  }
  fragLibDT_filenames <- list.files(path=tmpDir, pattern="FragLibDT.+fst$", full.names=T)
  fragLibDT <- data.table::rbindlist(lapply(fragLibDT_filenames, fst::read_fst, as.data.table=T))
  fragLibDT <- data.table::melt.data.table(
    fragLibDT,
    id.vars=c("ID", "FragLen", "Seed"),
    measure.vars=c("Deduplicated", "Weighted", "Mock"),
    variable.name="Library", value.name="Fragment"
  )
  fragLibDT[,Library:=paste0(Library, "_", FragLen, "_", Seed)][,FragLen:=NULL][,Seed:=NULL]
  fragLibDT <- data.table::dcast.data.table(fragLibDT, ID~Library, value.var="Fragment")
  fragLibDT[,ID:=NULL]
  gc();gc()
  return(fragLibDT)
}

