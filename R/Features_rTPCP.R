#' Repertoire-wide TCR-peptide contact profile (rTPCP) analysis.
#'
#' Calculates descriptive statistics of repertoire-wide TCR-peptide pairwise contact potentials.
#'
#' @param peptideSet A set of peptide sequences.
#' @param TCRSet Either a set of TCR sequences (as a character vector) or a list of sets of TCR sequences. If provided as a character vector, TCR sequences are randomly down-sampled to 10000, using \code{seedSet} to create a list. If provided as a list, it must be the same length with the seedSet.
#' @param fragLenSet A set of the lengths of TCR sequence fragments to be matched against the peptide set.
#' @param aaIndexIDSet A set of AAIndexIDs indicating the pairwise matching score matrix to be used. A set of AAIndex-derived matrices can be retrieved by \code{AACPMatrix}. Set "all" to select all AAIndexIDs.
#' @param alignTypeSet A set of alignment-type strings directly passed to the \code{type} argument of the \code{pairwiseAlignment} function in the \code{Biostrings} package.
#' @param TCRFragDepthSet A set of the numbers of TCR fragments to be matched. This should be kept constant for comparison.
#' @param seedSet A set of random seeds.
#' @importFrom dplyr %>%
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom magrittr set_colnames
#' @importFrom stringr str_sub
#' @importFrom stringr str_replace_all
#' @importFrom stringr str_split
#' @importFrom stringr fixed
#' @importFrom purrr flatten
#' @importFrom data.table as.data.table
#' @importFrom data.table transpose
#' @importFrom Kmisc str_rev
#' @importFrom Biostrings pairwiseAlignment
#' @importFrom Biostrings AA_STANDARD
#' @importFrom psych describe
#' @importFrom scales rescale
#' @importFrom parallel detectCores
#' @importFrom parallel splitIndices
#' @importFrom parallel detectCores
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterEvalQ
#' @importFrom parallel clusterExport
#' @importFrom parallel stopCluster
#' @importFrom pbapply pbapply
#' @importFrom pbapply pblapply
#' @export
#' @rdname Features
#' @name Features
Features_rTPCP <- function(
  peptideSet, TCRSet,
  fragLenSet=3:8, aaIndexIDSet="all",
  alignTypeSet="global-local", TCRFragDepthSet=10000,
  seedSet=1:5
){
  time.start <- proc.time()

  # AAIndex-derived pairwise matching matrices
  aaIndexMatrixFormat <- function(aaIndexID, pairMatInverse=T){
    pairMat <- dplyr::filter(AACPMatrix, AAIndexID==aaIndexID)$data[[1]]
    pairMat <- as.matrix(pairMat)
    rownames(pairMat) <- colnames(pairMat)
    if(pairMatInverse){ pairMat <- -pairMat }
    return(scales::rescale(pairMat))
  }
  if(identical(aaIndexIDSet, "all")) aaIndexIDSet <- AACPMatrix$AAIndexID
  pairMatSet <- c(lapply(aaIndexIDSet, function(a){aaIndexMatrixFormat(a, pairMatInverse=F)}),
                  lapply(aaIndexIDSet, function(a){aaIndexMatrixFormat(a, pairMatInverse=T)}))
  aaIndexIDSet <- c(aaIndexIDSet, paste0(aaIndexIDSet, "inv"))
  names(pairMatSet) <- aaIndexIDSet

  # TCR fragment dictionary.
  ## Note: parallelization is avoided because this process uses random seeds...
  fragmentDictionary <- function(sequenceSet, fragLen, depth, seed){
    set.seed(seed)
    s <- c(sequenceSet, Kmisc::str_rev(sequenceSet))
    f <- sapply(1:(max(nchar(s), na.rm=T)-fragLen+1), function(i){stringr::str_sub(s, i, i+fragLen-1)})
    f <- f[nchar(f)==fragLen]
    l <- length(f)
    pr <- table(f)/l
    sample(names(pr), size=depth, replace=T, prob=pr)
  }
  TCRFragDictSet <- function(TCRSet, fragLen, depth, seed){
    set.seed(seed)
    TCRSet.FR <- fragmentDictionary(TCRSet, fragLen, depth, seed)
    TCRSet.mock <- fragmentDictionary(
      sapply(TCRSet, function(tcr){paste0(sample(Biostrings::AA_STANDARD, size=nchar(tcr), replace=T), collapse="")}),
      fragLen, depth, seed
    )
    list("FR"=TCRSet.FR, "Mock"=TCRSet.mock)
  }
  message("Preparation of TCR fragment dictionary was started.")
  if(is.character(TCRSet)){
    TCRSet <- lapply(seedSet, function(s){set.seed(s); sample(TCRSet, length(TCRSet))})
  }else if(is.list(TCRSet)){
    if(length(TCRSet)!=length(seedSet)){
      print("The length of TCRSet (list) and that of seedSet are not matched!")
      return(NULL)
    }else{
      TCRSet <- mapply(function(tcr, s){set.seed(s); sample(tcr, length(tcr))}, TCRSet, seedSet, SIMPLIFY=F)
    }
  }
  names(TCRSet) <- seedSet
  TCRFragDictSet.List <- pbapply::pbapply(
    expand.grid(fragLenSet, TCRFragDepthSet, seedSet, stringsAsFactors=F), 1,
    function(v){TCRFragDictSet(TCRSet[[as.character(v[3])]], as.numeric(v[1]), as.numeric(v[2]), as.numeric(v[3]))}
  )
  TCRParameterSet <- apply(expand.grid(fragLenSet, TCRFragDepthSet, seedSet, stringsAsFactors=F), 1,
                           function(v){paste0(v[1], "_", v[2], "_", v[3])})
  names(TCRFragDictSet.List) <- TCRParameterSet
  message("Preparation of TCR fragment dictionary was finished.")

  # Fragment matching analysis

  ## Working function
  statSet <- c("mean","sd","median","trimmed","mad","skew","kurtosis","se","IQR","Q0.1","Q0.9")
  statNameSet <- c("Mean","SD","Median","TrimmedMean10","MedAbsDev","Skew","Kurtosis","SE","IQR","Q10","Q90")
  statNameSet <- expand.grid(statNameSet, c("FR", "Diff"))
  statNameSet <- apply(statNameSet, 1, function(v){paste0(v[1], "_", v[2])}) ## A character vector of length 22
  fragmentMatchingStats <- function(peptide, AAIndexID, alignType, TCRParameterString){
    s1 <- try(psych::describe(
      Biostrings::pairwiseAlignment(pattern=TCRFragDictSet.List[[TCRParameterString]][["FR"]], subject=peptide,
                                    substitutionMatrix=pairMatSet[[AAIndexID]], type=alignType, gapOpening=100, gapExtension=100,
                                    scoreOnly=T),
      trim=.1, interp=F, skew=T, type=3, ranges=T, IQR=T, quant=c(.10, .90)
    )[statSet], silent=T)
    s2 <- try(psych::describe(
      Biostrings::pairwiseAlignment(pattern=TCRFragDictSet.List[[TCRParameterString]][["Mock"]], subject=peptide,
                                    substitutionMatrix=pairMatSet[[AAIndexID]], type=alignType, gapOpening=100, gapExtension=100,
                                    scoreOnly=T),
      trim=.1, interp=F, skew=T, type=3, ranges=T, IQR=T, quant=c(.10, .90)
    )[statSet], silent=T)
    ## psych::describe returns the following values: vars,n,mean,sd,median,trimmed,mad,min,max,range,skew,kurtosis,se,IQR,Q0.1,Q0.9
    ## Of which the followings are kept: mean,sd,median,trimmed,mad,skew,kurtosis,se,IQR,Q0.1,Q0.9
    if(any(c(class(s1),class(s2))=="try-error")){
      message(paste0("Fragment matching on the peptide '", peptide, "' returned an unexpected error!"))
      return(rep(NA, length(statSet)*2))
    }
    return(as.numeric(c(s1, s1-s2))) ## A numeric vector of length 22
  }

  ## Combinations of parameters
  parameterGrid <- expand.grid(peptideSet, aaIndexIDSet, alignTypeSet, TCRParameterSet, stringsAsFactors=F) %>%
    magrittr::set_colnames(c("Peptide","AAIndexID","Alignment","TCRParameter"))
  paramCombN <- nrow(parameterGrid)
  message("Number of parameter combinations = ", paramCombN)
  gc();gc()

  ## Parallelized fragment matching
  message("Fragment matching was started. (Memory occupied = ", memory.size(), "[Mb])")
  cl <- parallel::makeCluster(parallel::detectCores(), type="SOCK")
  invisible(parallel::clusterEvalQ(cl, {
    library(psych)
    library(Biostrings)
  }))
  parallel::clusterExport(
    cl,
    list("parameterGrid", "fragmentMatchingStats", "TCRFragDictSet.List", "pairMatSet", "statSet"),
    envir=environment()
  )
  dt_feature <- pbapply::pblapply(
    1:paramCombN,
    function(i){fragmentMatchingStats(parameterGrid[i,1], parameterGrid[i,2], parameterGrid[i,3], parameterGrid[i,4])},
    cl=cl
  ) %>% data.table::as.data.table() %>% data.table::transpose() %>% magrittr::set_colnames(statNameSet)
  parallel::stopCluster(cl)
  message("Fragment matching was finished. (Memory occupied = ", memory.size(), "[Mb])")
  gc();gc()

  ## Final formatting
  message("Data formatting...")
  parameterGrid <- cbind(
    data.table::as.data.table(parameterGrid),
    data.table::as.data.table(stringr::str_split(parameterGrid$"TCRParameter", "_", simplify=T)) %>%
      magrittr::set_colnames(c("FragLen","Depth","Seed"))
  )
  parameterGrid <- parameterGrid[,"TCRParameter":=NULL]
  dt_feature <- cbind(parameterGrid, dt_feature)
  dt_feature_list <- dt_feature %>% split(by=c("AAIndexID", "FragLen"), sorted=T)
  dt_feature_list <- mapply(
    function(dt, nm){
      dt.meta <- dt[, c("Peptide", "AAIndexID", "Alignment", "FragLen", "Depth", "Seed")]
      dt.num <- dt[, c("Peptide", "AAIndexID", "Alignment", "FragLen", "Depth", "Seed"):=NULL]
      colnames(dt.num) <- stringr::str_replace_all(paste0(colnames(dt.num), "_", nm), stringr::fixed("."), "_")
      dt <- cbind(dt.meta, dt.num)[, c("AAIndexID", "FragLen"):=NULL]
      return(dt)
    }, dt_feature_list, names(dt_feature_list), SIMPLIFY=F)
  dt_feature_list <- purrr::flatten(lapply(dt_feature_list, function(dt){split(dt, by=c("Alignment", "Depth", "Seed"), keep.by=F, sorted=T)}))
  paramSet <- unique(names(dt_feature_list))
  dt_feature_list <- lapply(paramSet, function(paramString){
    dt.lst <- dt_feature_list[grep(paramString, names(dt_feature_list))]
    pept <- dt.lst[[1]][["Peptide"]]
    dt.lst <- lapply(dt.lst, function(dt){dt[,"Peptide":=NULL]})
    names(dt.lst) <- NULL
    dt <- do.call("cbind", dt.lst)
    colnames(dt) <- paste0("rTPCP_", colnames(dt))
    dt <- dt[, "Peptide":=pept]
    return(dt)
  })
  names(dt_feature_list) <- paramSet
  gc();gc()

  # Output
  time.end <- proc.time()
  message("Overall time required = ", format((time.end-time.start)[3], nsmall=2), "[sec]")
  return(dt_feature_list)
}
