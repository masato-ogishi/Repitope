#' Repertoire-wide TCR-peptide contact profile (rTPCP) analysis.
#' 
#' Calculates descriptive statistics of repertoire-wide TCR-peptide pairwise contact potentials.
#' 
#' @param peptideSet A set of peptide sequences.
#' @param TCRSet Either a set of TCR sequences (as a character vector) or a list of sets of TCR sequences. If provided as a character vector, TCR sequences are randomly down-sampled to 10000, using \code{seedSet} to create a list. If provided as a list, it must be the same length with the seedSet.
#' @param aaIndexIDSet A set of AAIndexIDs indicating the pairwise matching score matrix to be used. A set of AAIndex-derived matrices can be retrieved by \code{AACPMatrix}. Set "all" to select all AAIndexIDs.
#' @param alignTypeSet A set of alignment-type strings directly passed to the \code{type} argument of the \code{pairwiseAlignment} function in the \code{Biostrings} package.
#' @param fragLenSet A set of the lengths of TCR sequence fragments to be matched against the peptide set.
#' @param seedSet A set of random seeds.
#' @param TCRFragDepthSet A set of the numbers of TCR fragments to be matched. This should be kept constant for comparison.
#' @importFrom dplyr %>%
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom dplyr bind_cols
#' @importFrom tidyr separate
#' @importFrom magrittr set_colnames
#' @importFrom lubridate seconds_to_period
#' @importFrom lubridate now
#' @importFrom stringr str_sub
#' @importFrom stringr str_replace_all
#' @importFrom stringr fixed
#' @importFrom purrr flatten
#' @importFrom data.table as.data.table
#' @importFrom Kmisc str_rev
#' @importFrom Biostrings pairwiseAlignment
#' @importFrom Biostrings AA_STANDARD
#' @importFrom psych describe
#' @importFrom parallel splitIndices
#' @importFrom pbapply pbapply
#' @importFrom pbapply pblapply
#' @importFrom pbapply timerProgressBar
#' @importFrom pbapply setTimerProgressBar
#' @export
#' @rdname Features_rTPCP
#' @name Features_rTPCP
Features_rTPCP <- function(peptideSet, TCRSet, 
                           aaIndexIDSet="all", alignTypeSet="global-local", 
                           fragLenSet=5, TCRFragDepthSet=10000, seedSet=1:5){
  time.start <- proc.time()
  
  # TCR fragment dictionary to be matched. 
  # Note: parallelization is avoided because this process uses random seeds...
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
  names(pairMatSet) <- aaIndexIDSet

  # Fragment matching
  
  ## Working function
  statSet <- c("mean","sd","median","trimmed","mad","skew","kurtosis","se","IQR","Q0.1","Q0.9")
  statNameSet <- c("Mean","SD","Median","TrimmedMean10","MedAbsDev","Skew","Kurtosis","SE","IQR","Q10","Q90")
  statNameSet <- expand.grid(statNameSet, c("FR", "Diff"))
  statNameSet <- apply(statNameSet, 1, function(v){paste0(v[1], "_", v[2])}) ## A character vector of length 22
  fragmentMatchingStats <- function(peptide, AAIndexID, alignType, TCRParameterString){
    s1 <- try(psych::describe(
      Biostrings::pairwiseAlignment(pattern=TCRFragDictSet.List[[TCRParameterString]][["FR"]], subject=peptide, 
                                    substitutionMatrix=pairMatSet[[AAIndexID]], type=alignType, scoreOnly=T), 
      trim=.1, interp=F, skew=T, type=3, ranges=T, IQR=T, quant=c(.10, .90)
    )[statSet], silent=T)
    s2 <- try(psych::describe(
      Biostrings::pairwiseAlignment(pattern=TCRFragDictSet.List[[TCRParameterString]][["Mock"]], subject=peptide, 
                                    substitutionMatrix=pairMatSet[[AAIndexID]], type=alignType, scoreOnly=T), 
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
  parameterGrid <- expand.grid(peptideSet, aaIndexIDSet, alignTypeSet, TCRParameterSet, stringsAsFactors=F)
  paramCombN <- nrow(parameterGrid)
  message("Number of parameter combinations = ", paramCombN)
  gc();gc()
  
  ## Non-parallelized fragment matching calculation
  ### Note: parallelization using foreach is slower than substitution-only loop...
  ### Note: pbapply::timerProgressBar is a bit slow...
  message("Fragment matching was started. (Memory occupied = ", memory.size(), "[Mb])")
  mat <- matrix(nrow=paramCombN, ncol=length(statNameSet))
  id.lst <- parallel::splitIndices(paramCombN, 100)
  tm.elapsed <- 0
  for(j in 1:100){
    tm.start <- proc.time()
    for(i in id.lst[[j]]){
      mat[i,] <- fragmentMatchingStats(parameterGrid[i,1], parameterGrid[i,2], parameterGrid[i,3], parameterGrid[i,4])
    }
    tm.end <- proc.time()
    tm.elapsed <- round(tm.elapsed + (tm.end-tm.start)[3])
    tm.remain <- round(tm.elapsed*100/j - tm.elapsed)
    cat("\r ", j, "% elapsed = ", as.character(lubridate::seconds_to_period(tm.elapsed)), ", remaining ~ ", as.character(lubridate::seconds_to_period(tm.remain)), sep="")
  }
  mat <- as.data.frame(mat) %>% magrittr::set_colnames(statNameSet)
  parameterGrid <- parameterGrid %>%
    magrittr::set_colnames(c("Peptide","AAIndexID","Alignment","TCRParameter"))
  ymdt <- stringr::str_replace_all(stringr::str_replace_all(stringr::str_replace_all(lubridate::now(), stringr::fixed(":"), "."), " ", "."), "-", ".")
  readr::write_csv(dplyr::bind_cols(parameterGrid, mat), paste0("FragmentMatchMatrix_", ymdt, ".csv.gz"))
  message("Fragment matching was finished. (Memory occupied = ", memory.size(), "[Mb])")
  gc();gc()
  
  ## Final formatting
  message("Data formatting...")
  parameterGrid <- parameterGrid %>%
    tidyr::separate(TCRParameter, into=c("FragLen", "Depth", "Seed"), sep="_") %>%
    dplyr::mutate(FragLen=as.numeric(FragLen), Depth=as.numeric(Depth), Seed=as.numeric(Seed))
  gc();gc()
  dt_feature_list <- dplyr::bind_cols(parameterGrid, mat) %>%
    data.table::as.data.table() %>%
    split(by=c("AAIndexID", "FragLen"), sorted=T)
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
