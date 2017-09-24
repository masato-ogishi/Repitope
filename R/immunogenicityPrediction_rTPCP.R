#' Repertoire-wide TCR-peptide contact profile (rTPCP) analysis.
#' 
#' Calculates descriptive statistics of repertoire-wide TCR-peptide pairwise contact potentials.
#' 
#' @param peptideSet A set of peptide sequences.
#' @param TCRSet A set of TCR sequences.
#' @param aaIndexIDSet A set of AAIndexIDs indicating the pairwise matching score matrix to be used. A set of AAIndex-derived matrices can be retrieved by \code{AACPMatrix}. Set "all" to select all AAIndexIDs.
#' @param alignTypeSet A set of alignment-type strings directly passed to the \code{type} argument of the \code{pairwiseAlignment} function in the \code{Biostrings} package.
#' @param fragLenSet A set of the lengths of TCR sequence fragments to be matched against the peptide set.
#' @param seedSet A set of random seeds.
#' @param TCRFragDepthSet A set of the numbers of TCR fragments to be matched. This should be kept constant for comparison.
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
#' @importFrom purrr set_names
#' @importFrom magrittr set_colnames
#' @importFrom Biostrings pairwiseAlignment
#' @importFrom Biostrings AA_STANDARD
#' @importFrom psych describe
#' @importFrom parallel splitIndices
#' @importFrom parallel detectCores
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterEvalQ
#' @importFrom pbapply pbapply
#' @importFrom parallel clusterExport
#' @importFrom parallel stopCluster
#' @importFrom ff ff 
#' @export
#' @rdname rTPCP
#' @name rTPCP
rTPCPAnalysis <- function(peptideSet, TCRSet, 
                          aaIndexIDSet="all", alignTypeSet="global-local", 
                          fragLenSet=5, TCRFragDepthSet=10000, seedSet=1:5){
  time.start <- proc.time()
  
  # TCR fragment dictionary to be matched. 
  # Note: parallelization is avoided because this process needs random seeds...
  fragmentDictionary <- function(sequenceSet, fragLen, depth, seed){
    set.seed(seed)
    s <- sequenceFilter(c(sequenceSet, Kmisc::str_rev(sequenceSet)))
    f <- sapply(1:(max(nchar(s), na.rm=T)-fragLen+1), function(i){stringr::str_sub(s, i, i+fragLen-1)})
    f <- f[nchar(f)==fragLen]
    l <- length(f)
    pr <- table(f)/l
    sample(names(pr), size=depth, replace=T, prob=pr)
  }
  TCRFragDictSet <- function(TCRSet, fragLen, TCRFragDepth, seed){
    TCRSet.FR <- fragmentDictionary(TCRSet, fragLen, TCRFragDepth, seed)
    TCRSet.mock <- fragmentDictionary(
      sapply(TCRSet, function(tcr){paste0(sample(Biostrings::AA_STANDARD, size=nchar(tcr), replace=T), collapse="")}),
      fragLen, TCRFragDepth, seed
    )
    list("FR"=TCRSet.FR, "Mock"=TCRSet.mock)
  }
  TCRFragDictSet.List <- apply(expand.grid(fragLenSet, TCRFragDepthSet, seedSet, stringsAsFactors=F), 1,
                               function(v){TCRFragDictSet(TCRSet, as.numeric(v[[1]]), as.numeric(v[[2]]), as.numeric(v[[3]]))})
  TCRParameterSet <- apply(expand.grid(fragLenSet, TCRFragDepthSet, seedSet, stringsAsFactors=F), 1,
                           function(v){paste0(v[[1]], "_", v[[2]], "_", v[[3]])})
  TCRFragDictSet.List <- purrr::set_names(TCRFragDictSet.List, nm=TCRParameterSet)
    
  # AAIndex-derived pairwise matching matrices
  aaIndexMatrixFormat <- function(aaIndexID, pairMatInverse=T){
    pairMat <- dplyr::filter(AACPMatrix, AAIndexID==aaIndexID)$data[[1]]
    pairMat <- as.matrix(pairMat)
    rownames(pairMat) <- colnames(pairMat)
    if(pairMatInverse){ pairMat <- -pairMat }
    pairMat
  }
  if(identical(aaIndexIDSet, "all")) aaIndexIDSet <- AACPMatrix$AAIndexID
  pairMatSet <- c(lapply(aaIndexIDSet, function(a){aaIndexMatrixFormat(a, F)}), lapply(aaIndexIDSet, function(a){aaIndexMatrixFormat(a, T)}))
  aaIndexIDSet <- c(aaIndexIDSet, paste0(aaIndexIDSet, "inv"))
  pairMatSet <- purrr::set_names(pairMatSet, nm=aaIndexIDSet)
  
  # Parallelization
  cl <- parallel::makeCluster(parallel::detectCores(), type="SOCK")
  invisible(parallel::clusterEvalQ(cl, {
    library(ff)
    library(Biostrings)
    library(psych)
  }))
  
  # Fragment matching
  fragmentMatchingStats <- function(peptide, AAIndexID, alignType, TCRParameterString){
    s1 <- as.numeric(psych::describe(
      Biostrings::pairwiseAlignment(pattern=TCRFragDictSet.List[[TCRParameterString]][["FR"]], subject=peptide, 
                                    substitutionMatrix=pairMatSet[[AAIndexID]], type=alignType, scoreOnly=T), 
      trim=.1, interp=F, skew=T, type=3, ranges=T, IQR=T, quant=c(.10, .25, .75, .90)
    ))[c(6, 7, 14, 15, 16, 5, 17, 18, 11, 12)]
    s2 <- as.numeric(psych::describe(
      Biostrings::pairwiseAlignment(pattern=TCRFragDictSet.List[[TCRParameterString]][["Mock"]], subject=peptide, 
                                    substitutionMatrix=pairMatSet[[AAIndexID]], type=alignType, scoreOnly=T), 
      trim=.1, interp=F, skew=T, type=3, ranges=T, IQR=T, quant=c(.10, .25, .75, .90)
    ))[c(6, 7, 14, 15, 16, 5, 17, 18, 11, 12)]
    ## trimmed, mad, IQR, Q0.1, Q0.25, median, Q0.75, Q0.9, skew, kurtosis
    ff::ff(c(s1, s1-s2), vmode="double") ## A vector of length 20
  }
  statNameSet <- c("TrimmedMean10","MedAbsDev","IQR","Q10","Q25","Q50","Q75","Q90","Skew","Kurtosis")
  statNameSet <- expand.grid(statNameSet, c("FR", "Diff"))
  statNameSet <- apply(statNameSet, 1, function(v){paste0(v[1], "_", v[2])})
  parameterGrid <- expand.grid(peptideSet, aaIndexIDSet, alignTypeSet, TCRParameterSet, stringsAsFactors=F)
  message(paste0("Number of parameter combinations = ", nrow(parameterGrid)))
  message(paste0("Parallelized fragment matching was started. (Memory occupied = ", memory.size(), "[Mb])"))
  df_feature <- pbapply::pbapply(parameterGrid, 1, 
                                 function(v){fragmentMatchingStats(v[1], v[2], v[3], v[4])},
                                 cl=cl)
  message(paste0("Parallelized fragment matching was finished. (Memory occupied = ", memory.size(), "[Mb])"))
  gc();gc();
  message(paste0("Integration of fragment matching results was started. (Memory occupied = ", memory.size(), "[Mb])"))
  df_feature <- as.numeric(unlist(lapply(df_feature, function(ffVector){ffVector[]}))) ## Convert a list of ff vectors into a combined normal numeric vector
  df_feature <- dplyr::bind_cols(parameterGrid, as.data.frame(t(matrix(df_feature, nrow=length(statNameSet))))) ## The row number 10 corresponds to the number of descriptive statistics retained.
  message(paste0("Integration of fragment matching results was finished. (Memory occupied = ", memory.size(), "[Mb])"))
  gc();gc();
  message(paste0("Data formatting..."))
  df_feature <- df_feature %>%
    magrittr::set_colnames(c("Peptide","AAIndexID","Alignment","TCRParameter",statNameSet)) %>%
    tidyr::gather(Stat, Value, -Peptide, -AAIndexID, -Alignment, -TCRParameter) %>%
    dplyr::transmute(Peptide, Alignment, TCRParameter, Feature=paste0(AAIndexID, "_", Stat), Value) %>%
    tidyr::separate(TCRParameter, into=c("FragLen", "Depth", "Seed"), sep="_") %>%
    tidyr::spread(Feature, Value)
  
  # Output
  parallell::stopCluster(cl)
  time.end <- proc.time()
  message(paste0("Overall time required = ", (time.end-time.start)[3], "[sec]"))
  return(df_feature)
}
