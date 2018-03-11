#' Generate peptide features.
#'
#' \code{Features_PeptDesc} calculates descriptive statistics of peptide descriptors.\cr
#' \code{Features_CPP} calculates descriptive statistics of interresidue pairwise contact potential profiles.\cr
#' \code{Features} is a wrapper function.
#'
#' @param peptideSet A set of peptide sequences.
#' @param fragLenSet A set of sliding window sizes. Must be between 3 and 8.
#' @param aaIndexIDSet A set of AAIndex IDs indicating the AACP scales to be used. Set "all" to shortcut the selection of all available AACP scales.
#' @param seedSet A set of random seeds.
#' @param coreN The number of cores to be used for parallelization. Set \code{NULL} to skip parallelization.
#' @param tmpDir Destination directory to save intermediate files.
#' @importFrom dplyr %>%
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom dplyr left_join
#' @importFrom tidyr gather
#' @importFrom tidyr spread
#' @importFrom tidyr unite
#' @importFrom data.table :=
#' @importFrom data.table as.data.table
#' @importFrom data.table rbindlist
#' @importFrom data.table setorder
#' @importFrom data.table setcolorder
#' @importFrom data.table melt.data.table
#' @importFrom data.table dcast.data.table
#' @importFrom fst read_fst
#' @importFrom fst write_fst
#' @importFrom magrittr set_names
#' @importFrom magrittr set_colnames
#' @importFrom Biostrings reverse
#' @importFrom Biostrings pairwiseAlignment
#' @importFrom Biostrings AA_STANDARD
#' @importFrom stringr str_sub
#' @importFrom stringr str_replace_all
#' @importFrom stringr str_split
#' @importFrom stringr str_detect
#' @importFrom stringr fixed
#' @importFrom psych describe
#' @importFrom matrixStats rowMins
#' @importFrom matrixStats rowMaxs
#' @importFrom matrixStats rowMeans2
#' @importFrom matrixStats rowMedians
#' @importFrom parallel detectCores
#' @importFrom parallel detectCores
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterEvalQ
#' @importFrom parallel clusterExport
#' @importFrom parallel stopCluster
#' @importFrom pbapply pbapply
#' @importFrom pbapply pblapply
#' @importFrom snow parLapply
#' @importFrom bigmemory as.big.matrix
#' @importFrom bigmemory attach.big.matrix
#' @import Peptides
#' @export
#' @rdname Features
#' @name Features_PeptDesc
Features_PeptDesc <- function(
  peptideSet,
  fragLenSet=3:8
){
  # Start calculation
  set.seed(12345)
  time.start <- proc.time()

  # Working functions
  peptideDescriptor.Batch <- function(peptide){
    unlist(list(
      Peptides::aIndex(peptide),
      Peptides::blosumIndices(peptide),
      Peptides::boman(peptide),
      Peptides::charge(peptide),
      Peptides::crucianiProperties(peptide),
      Peptides::fasgaiVectors(peptide),
      Peptides::hmoment(peptide),
      Peptides::hydrophobicity(peptide),
      Peptides::instaIndex(peptide),
      Peptides::kideraFactors(peptide),
      Peptides::mswhimScores(peptide),
      Peptides::pI(peptide),
      Peptides::protFP(peptide),
      Peptides::vhseScales(peptide),
      Peptides::zScales(peptide)
    ))
  }
  peptideDescriptor.NameSet <- c("AliphaticIndex","BLOSUM1","BLOSUM2","BLOSUM3","BLOSUM4","BLOSUM5","BLOSUM6","BLOSUM7","BLOSUM8","BLOSUM9","BLOSUM10","Boman","Charge","PP1","PP2","PP3","F1","F2","F3","F4","F5","F6","HydrophobicMoment","Hydrophobicity","Instability","KF1","KF2","KF3","KF4","KF5","KF6","KF7","KF8","KF9","KF10","MSWHIM1","MSWHIM2","MSWHIM3","pI","ProtFP1","ProtFP2","ProtFP3","ProtFP4","ProtFP5","ProtFP6","ProtFP7","ProtFP8","VHSE1","VHSE2","VHSE3","VHSE4","VHSE5","VHSE6","VHSE7","VHSE8","Z1","Z2","Z3","Z4","Z5")
  peptideDescriptor.FragStat.Single <- function(peptide, fragLen){
    f <- sapply(1:(max(nchar(peptide))-fragLen+1), function(i){stringr::str_sub(peptide, i, i+fragLen-1)})
    d <- sapply(f[nchar(f)==fragLen], function(s){peptideDescriptor.Batch(s)})
    data.frame("Peptide"=peptide,
               "FragLen"=fragLen,
               "AADescriptor"=peptideDescriptor.NameSet,
               "Min"=matrixStats::rowMins(d),
               "Max"=matrixStats::rowMaxs(d),
               "Mean"=matrixStats::rowMeans2(d),
               "Median"=matrixStats::rowMedians(d))
  }

  # Parallelized calculation of descriptive statistics
  parameterGrid <- expand.grid(peptideSet, fragLenSet, stringsAsFactors=F)
  cl <- parallel::makeCluster(parallel::detectCores(), type="SOCK")
  parallel::clusterExport(
    cl=cl,
    list("parameterGrid", "peptideDescriptor.FragStat.Single", "peptideDescriptor.Batch", "peptideDescriptor.NameSet"),
    envir=environment()
  )
  df_feature <- pbapply::pbapply(
    parameterGrid, 1,
    function(v){peptideDescriptor.FragStat.Single(v[1], as.numeric(v[2]))},
    cl=cl
  )
  parallel::stopCluster(cl)
  gc();gc();

  # Formatting
  df_feature <- data.table::rbindlist(df_feature) %>%
    tidyr::gather(Stat, Value, -Peptide, -FragLen, -AADescriptor) %>%
    tidyr::unite(Feature, AADescriptor, Stat, FragLen, sep="_") %>%
    tidyr::spread(Feature, Value)
  colnames(df_feature) <- paste0("PeptDesc_", colnames(df_feature))
  colnames(df_feature)[1] <- "Peptide"

  # Basic peptide information
  df_basic <- data.frame("Peptide"=peptideSet, "Peptide_Length"=nchar(peptideSet))
  for(aa in Biostrings::AA_STANDARD){
    df_basic[[paste0("Peptide_Contain", aa)]] <- as.numeric(stringr::str_detect(peptideSet, aa))
  }
  df_feature <- suppressWarnings(dplyr::left_join(df_basic, df_feature, by="Peptide"))
  rownames(df_feature) <- 1:nrow(df_feature)

  # Output
  time.end <- proc.time()
  message("Overall time required = ", format((time.end-time.start)[3], nsmall=2), "[sec]")
  return(data.table::as.data.table(df_feature))
}

#' @export
#' @rdname Features
#' @name Features_CPP
Features_CPP <- function(
  peptideSet,
  aaIndexIDSet="all",
  fragLenSet=3:8,
  fragDepthSet=100000,
  seedSet=1:5,
  coreN=parallel::detectCores(),
  tmpDir=tempdir()
){
  # Temp directory check
  dir.create(tmpDir, showWarnings=F, recursive=T)

  # Start calculation
  set.seed(12345)
  time.start <- proc.time()

  # Pairwise contact potential matrix
  AACPMatrix <- CPP_AACPMatrix()
  if(identical(aaIndexIDSet, "all")){
    aaIndexIDSet <- unique(AACPMatrix$AAIndexID)
  }else{
    AACPMatrix <- AACPMatrix[AAIndexID %in% aaIndexIDSet,]
  }
  AACPMatrix[["AAIndexID"]] <- as.numeric(factor(AACPMatrix[["AAIndexID"]], levels=aaIndexIDSet))
  AACPMatrix_binfile <- "AACPMatrix.bin"
  AACPMatrix_descfile <- "AACPMatrix.desc"
  if(!file.exists(file.path(tmpDir, AACPMatrix_descfile))){
    AACPMatrix <- bigmemory::as.big.matrix(
      as.matrix(AACPMatrix),
      type="double",
      separated=F,
      backingpath=tmpDir,
      backingfile=AACPMatrix_binfile,
      descriptorfile=AACPMatrix_descfile
    )
  }

  # Fragment library for matching
  FragLib_TCR <- fst::read_fst(system.file("TCRFragmentLibrary.fst", package="Repitope"), as.data.table=T) ## fragmentLibrary(TCRSet_Public, 3:8)
  FragLib_TCR <- FragLib_TCR[FragmentLength %in% fragLenSet,]
  Frag_Levels <- FragLib_TCR$Fragment
  FragLib_TCR$Fragment <- as.numeric(factor(FragLib_TCR$Fragment, levels=Frag_Levels))
  FragLib_binfile <- "FragLib.bin"
  FragLib_descfile <- "FragLib.desc"
  if(!file.exists(file.path(tmpDir, FragLib_descfile))){
    FragLib_TCR <- bigmemory::as.big.matrix(
      as.matrix(FragLib_TCR), type="double", separated=F,
      backingpath=tmpDir,
      backingfile=FragLib_binfile,
      descriptorfile=FragLib_descfile
    )
  }

  # Contact potential profiling
  ## Working functions
  load_AACPMatrix <- function(aaIndexID){
    AACPMatrix <- bigmemory::attach.big.matrix(obj=AACPMatrix_descfile, path=tmpDir)
    AACPMatrix <- AACPMatrix[AACPMatrix[, "AAIndexID"]==grep(paste0(aaIndexID, "$"), aaIndexIDSet),]
    AACPMatrix <- AACPMatrix[, Biostrings::AA_STANDARD]
    rownames(AACPMatrix) <- Biostrings::AA_STANDARD
    return(AACPMatrix)
  }
  load_FragLib <- function(fragLen, fragDepth=NULL, seed=12345){
    FragLib <- bigmemory::attach.big.matrix(obj=FragLib_descfile, path=tmpDir)
    FragLib <- FragLib[FragLib[,"FragmentLength"]==fragLen,]
    FragSet <- Frag_Levels[FragLib[,"Fragment"]]
    if(is.null(fragDepth)) fragDepth <- sum(FragLib[,"Count"])
    set.seed(seed); FragSet.Weighted <- sample(FragSet, size=fragDepth, replace=T, prob=FragLib[,"Freq"])
    set.seed(seed); FragSet <- sample(FragSet, size=fragDepth, replace=T)
    set.seed(seed); FragSet.Mock <- sapply(1:fragDepth, function(i){paste0(sample(Biostrings::AA_STANDARD, size=fragLen, replace=T), collapse="")})
    return(list(FragSet, FragSet.Weighted, FragSet.Mock))
  }
  calculate_stats_single <- function(peptide, fragSet, aacpMat, seed=12345){
    set.seed(seed)
    al <- Biostrings::pairwiseAlignment(
      pattern=fragSet, subject=peptide,
      substitutionMatrix=aacpMat, type="global-local",
      gapOpening=100, gapExtension=100,
      scoreOnly=T
    )
    statSet <- c("mean","sd","median","trimmed","mad","skew","kurtosis","se","IQR","Q0.1","Q0.9")
    statDF <- psych::describe(al, trim=.1, interp=F, skew=T, type=3, ranges=T, IQR=T, quant=c(.10, .90))[statSet]
    colnames(statDF) <- c("Mean","SD","Med","TrM","MAD","Skew","Kurt","SE","IQR","Q10","Q90")
    return(statDF)
  }
  contact_potential_profiling <- function(peptideSet, aaIndexID, fragLen, fragDepth, seed, cl=cl){
    parallel::clusterExport(
      cl=cl,
      list("aaIndexID", "fragLen", "fragDepth", "seed"),
      envir=environment()
    )
    dt <- snow::parLapply(
      cl=cl,
      peptideSet,
      function(pept){
        set.seed(seed)
        AACPMatrix_Tmp <- load_AACPMatrix(aaIndexID)
        FragLib_Tmp <- load_FragLib(fragLen, fragDepth, seed)
        df_meta <- data.frame("Peptide"=pept, "AAIndexID"=aaIndexID, "FragLen"=fragLen, "FragDepth"=fragDepth, "Seed"=seed, stringsAsFactors=F)
        df <- data.table::rbindlist(lapply(FragLib_Tmp, function(d){cbind(df_meta, calculate_stats_single(pept, d, AACPMatrix_Tmp, seed))}))
        df$"Library" <- c("Deduplicated", "Weighted", "Mock")
        return(df)
      }
    )
    dt <- data.table::as.data.table(data.table::rbindlist(dt))
    data.table::setcolorder(dt, c("Peptide", "AAIndexID", "FragLen", "FragDepth", "Library", "Seed", "Mean","SD","Med","TrM","MAD","Skew","Kurt","SE","IQR","Q10","Q90"))
    return(dt)
  }

  ## Main
  parameterGrid <- expand.grid(aaIndexIDSet, fragLenSet, fragDepthSet, seedSet, stringsAsFactors=F) %>%
    magrittr::set_colnames(c("AAIndexID","FragLen","FragDepth","Seed"))
  paramCombN <- nrow(parameterGrid)
  message("Number of parameter combinations = ", paramCombN)
  message("Launching parallel clusters...")
  if(is.null(coreN)) coreN <- 1
  cl <- parallel::makeCluster(coreN, type="SOCK")
  parallel::clusterExport(
    cl=cl,
    list("load_AACPMatrix", "AACPMatrix_descfile", "aaIndexIDSet",
         "load_FragLib", "FragLib_descfile", "Frag_Levels",
         "calculate_stats_single", "tmpDir"),
    envir=environment()
  )
  message("Contact potential profiling...")
  dt_cpp <- pbapply::pblapply(
    1:paramCombN,
    function(i){
      out <- file.path(tmpDir, paste0("dt_feature_cpp_", i, ".fst"))
      if(!file.exists(out)){
        dt <- contact_potential_profiling(peptideSet, parameterGrid[i,1], parameterGrid[i,2], parameterGrid[i,3], parameterGrid[i,4], cl=cl)
        rownames(dt) <- 1:nrow(dt)
        fst::write_fst(dt, out)
      }
    }
  )
  parallel::stopCluster(cl)

  ## Final formatting
  message("Data formatting...")
  dt_cpp <- pbapply::pblapply(
    list.files(pattern="^dt_feature_cpp.+fst$", path=tmpDir, full.names=T), fst::read_fst, as.data.table=T
  ) %>% data.table::rbindlist() %>% data.table::as.data.table()
  data.table::setorder(dt_cpp, Peptide, AAIndexID, FragLen, FragDepth, Library, Seed)
  col_id <- c("Peptide", "AAIndexID", "FragLen", "FragDepth", "Library", "Seed")
  col_val <- setdiff(colnames(dt_cpp), col_id)
  dt_cpp <- data.table::melt.data.table(dt_cpp, id=col_id, measure=col_val, variable.name="Stat", value.name="Value")
  dt_cpp[,"Feature":=paste0("CPP_", AAIndexID, "_", Stat, "_", FragLen)][,"AAIndexID":=NULL][,"FragLen":=NULL][,"Stat":=NULL]
  dt_cpp <- data.table::dcast.data.table(dt_cpp, Peptide+FragDepth+Library~Feature, value.var="Value", fun.aggregate=mean)

  # Finish the timer
  time.end <- proc.time()
  message("Overall time required = ", format((time.end-time.start)[3], nsmall=2), "[sec]")

  # Clear logs
  message("Erase the temporary folder...")
  rm(list=setdiff(ls(), c("tmpDir", "AACPMatrix_binfile", "FragLib_binfile", "dt_cpp")))
  gc();gc()
  file.remove(file.path(tmp.dir, AACPMatrix_binfile))
  file.remove(file.path(tmp.dir, FragLib_binfile))
  file.remove(list.files(pattern="^dt_feature_cpp.+fst$", path=tmpDir, full.names=T))

  # Output
  return(dt_cpp)
}

#' @export
#' @rdname Features
#' @name Features_CPP
Features <- function(
  peptideSet,
  aaIndexIDSet="all",
  fragLenSet=3:8,
  fragDepthSet=100000,
  seedSet=1:5,
  coreN=parallel::detectCores(),
  tmpDir=tempdir()
){
  dt_peptdesc <- Features_PeptDesc(peptideSet, fragLenSet)
  dt_cpp <- Features_CPP(peptideSet, aaIndexIDSet, fragLenSet, fragDepthSet, seedSet, coreN, tmpDir)
  dt <- merge(dt_peptdesc, dt_cpp, by="Peptide")
  dt[,"FragDepth":=format(FragDepth, scientific=F)]
  return(split(dt, by=c("Library", "FragDepth")))
}
