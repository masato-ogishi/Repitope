#' Repertoire-wide TCR-peptide contact profile (rTPCP) analysis.
#' 
#' \code{fragmentation} is a utility function which generates fragments (substrings of fixed length) from input sequences.\cr
#' \code{tcrFragmentDictionary} creates the TCR fragment dictionary to be matched.\cr
#' \code{rTPCPAnalysis} calculates descriptive statistics of repertoire-wide TCR-peptide pairwise contact potentials.
#' 
#' @param sequenceSet A set of amino acid sequences.
#' @param peptideSet A set of peptide sequences.
#' @param TCRSet A set of TCR sequences.
#' @param aaIndexIDSet A set of AAIndexIDs indicating the pairwise matching score matrix to be used. A set of AAIndex-derived matrices can be retrieved by \code{AACPMatrix}. Set "all" to select all AAIndexIDs.
#' @param fragLen The length of TCR sequence fragments to be matched against the peptide set.
#' @param alignTypeSet A set of alignment-type strings directly passed to the \code{type} argument of the \code{pairwiseAlignment} function in the \code{Biostrings} package.
#' @param seed A random seed.
#' @param TCRFragDepth The number of TCR fragments to be matched. This should be kept constant throughout an analysis to keep repertoire-dependent fluctuation of descriptive statistics minimized.
#' @importFrom stringr str_sub
#' @importFrom Kmisc str_rev
#' @importFrom purrr set_names
#' @importFrom dplyr %>%
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom dplyr bind_cols
#' @importFrom tidyr gather
#' @importFrom tidyr spread
#' @importFrom tidyr separate
#' @importFrom tidyr unite
#' @importFrom magrittr set_colnames
#' @importFrom Biostrings pairwiseAlignment
#' @importFrom Biostrings AA_STANDARD
#' @importFrom psych describe
#' @importFrom parallel splitIndices
#' @importFrom parallel detectCores
#' @importFrom snowfall sfInit
#' @importFrom snowfall sfLibrary
#' @importFrom snowfall sfExport
#' @importFrom snowfall sfApply
#' @importFrom snowfall sfStop
#' @importFrom ff ff 
#' @export
#' @rdname rTPCP
#' @name rTPCP
fragmentation <- function(sequenceSet, fragLen){
  s <- sequenceFilter(sequenceSet)
  f <- sapply(1:(max(nchar(s))-fragLen+1), function(i){stringr::str_sub(s, i, i+fragLen-1)})
  f[nchar(f)==fragLen]
}

#' @export
#' @rdname rTPCP
#' @name rTPCP
tcrFragmentDictionary <- function(TCRSet, fragLen=5, TCRFragDepth=10000, seed=12345){
  set.seed(seed)
  c(TCRSet, Kmisc::str_rev(TCRSet)) %>% 
    fragmentation(fragLen) %>% 
    (function(tcr){l <- length(tcr); pr <- table(tcr)/l; sample(names(pr), size=TCRFragDepth, replace=T, prob=pr)})
}

#' @export
#' @rdname rTPCP
#' @name rTPCP
rTPCPAnalysis <- function(peptideSet, TCRSet, 
                          aaIndexIDSet="all", fragLen=5, alignTypeSet="global-local", 
                          seed=12345,
                          TCRFragDepth=10000){
  # Preparation of TCR fragment dictionaries
  time.start <- proc.time()
  set.seed(seed)
  TCRSet.FR <- tcrFragmentDictionary(TCRSet, fragLen, TCRFragDepth, seed)
  TCRSet.mock <- tcrFragmentDictionary(sapply(TCRSet, function(tcr){paste0(sample(Biostrings::AA_STANDARD, size=nchar(tcr), replace=T), collapse="")}),
                                       fragLen, TCRFragDepth, seed)
  TCRFragDictSet <- list("FR"=TCRSet.FR, "Mock"=TCRSet.mock)
  
  # Preparation of AAIndex-derived pairwise matching matrices
  aaIndexMatrixFormat <- function(aaIndexID, pairMatInverse=T){
    pairMat <- dplyr::filter(AACPMatrix, AAIndexID==aaIndexID)$data[[1]]
    pairMat <- as.matrix(pairMat)
    rownames(pairMat) <- colnames(pairMat)
    if(pairMatInverse){ pairMat <- -pairMat }
    pairMat
  }
  suppressWarnings(if(aaIndexIDSet=="all"){ id <- AACPMatrix$AAIndexID }else{ id <- aaIndexIDSet })
  pairMatSet <- c(lapply(id, function(a){aaIndexMatrixFormat(a, F)}), lapply(id, function(a){aaIndexMatrixFormat(a, T)}))
  id <- c(id, paste0(id, "inv"))
  pairMatSet <- purrr::set_names(pairMatSet, nm=id)
  
  # Fragment matching
  snowfall::sfInit(parallel=T, cpus=parallel::detectCores(), type="SOCK")
  sink(tempfile())
  suppressMessages(snowfall::sfLibrary(ff))
  suppressMessages(snowfall::sfLibrary(Biostrings))
  suppressMessages(snowfall::sfLibrary(psych))
  sink()
  fragmentMatchingStats <- function(v){ ## v = vector(peptide, AAIndexID, TCRRepDictID, alignType)
    s <- as.numeric(psych::describe(
      Biostrings::pairwiseAlignment(pattern=TCRFragDictSet[[v[3]]], subject=v[1], 
                                    substitutionMatrix=pairMatSet[[v[2]]], type=v[4], scoreOnly=T), 
      trim=.1, interp=F, skew=T, type=3, ranges=T, IQR=T, quant=c(.10, .25, .75, .90)
    ))[c(6, 7, 14, 15, 16, 5, 17, 18, 11, 12)]
    ## trimmed, mad, IQR, Q0.1, Q0.25, median, Q0.75, Q0.9, skew, kurtosis
    ff::ff(s, vmode="double")
  }
  snowfall::sfExport(list=c("fragmentMatchingStats","TCRFragDictSet","pairMatSet"))
  parameterMatrix <- expand.grid(peptideSet, id, c("FR", "Mock"), alignTypeSet, stringsAsFactors=F)
  message(paste0("Number of parameter combinations = ", nrow(parameterMatrix)))
  message(paste0("Parallelized fragment matching was started. (Memory occupied = ", memory.size(), "[Mb])"))
  df_feature <- snowfall::sfApply(parameterMatrix, 1, fragmentMatchingStats)
  message(paste0("Parallelized fragment matching was finished. (Memory occupied = ", memory.size(), "[Mb])"))
  gc();gc();
  message(paste0("Integration of fragment matching results was started. (Memory occupied = ", memory.size(), "[Mb])"))
  df_feature <- as.numeric(unlist(lapply(df_feature, function(ffVector){ffVector[]}))) ## Convert a list of ff vectors into a combined normal numeric vector
  df_feature <- dplyr::bind_cols(parameterMatrix, as.data.frame(t(matrix(df_feature, nrow=10)))) ## The row number 10 corresponds to the number of descriptive statistics retained.
  message(paste0("Integration of fragment matching results was finished. (Memory occupied = ", memory.size(), "[Mb])"))
  gc();gc();
  message(paste0("Data formatting..."))
  df_feature <- df_feature %>%
    magrittr::set_colnames(c("Peptide","AAIndexID","TCRRep","Alignment","TrimmedMean10","MedAbsDev","IQR","Q10","Q25","Q50","Q75","Q90","Skew","Kurtosis")) %>%
    dplyr::mutate(FragLen=fragLen) %>%
    tidyr::gather(Stat, Value, -Peptide, -AAIndexID, -TCRRep, -FragLen, -Alignment) %>%
    tidyr::spread(TCRRep, Value) %>%
    dplyr::mutate(Diff=FR-Mock) %>%
    dplyr::select(-Mock) %>%
    tidyr::gather(TCRRep, Value, -Peptide, -AAIndexID, -FragLen, -Alignment, -Stat) %>%
    tidyr::unite(Feature, c(AAIndexID, TCRRep, FragLen, Alignment, Stat), sep="_") %>%
    tidyr::spread(Feature, Value) %>%
    magrittr::set_colnames(gsub("-","",colnames(.),fixed=T))
  snowfall::sfStop()
  
  # Output
  time.end <- proc.time()
  message(paste0("Overall time required = ", (time.end-time.start)[3], "[sec]"))
  return(df_feature)
}
