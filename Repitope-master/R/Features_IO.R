#' Importing and exporting utilities for feature dataframes.
#'
#' @param featureDFFileNames A character vector of binary file names with either ".fst" or ".rds" extensions.
#' @param featureDFList A named list of feature dataframes.
#' @param fileNameHeader A file name header.
#' @export
#' @rdname Features_IO
#' @name Features_IO
readFeatureDFList <- function(featureDFFileNames){
  paramSet <- gsub(".rds$", "", gsub(".fst$", "", basename(featureDFFileNames)))
  paramSet <- data.table::transpose(as.data.frame(strsplit(paramSet, "_", fixed = T), stringsAsFactors = F))
  paramSet <- as.character(dplyr::last(paramSet))
  pos.fst <- which(stringr::str_detect(basename(featureDFFileNames), ".fst"))
  pos.rds <- which(stringr::str_detect(basename(featureDFFileNames), ".rds"))
  files.fst <- featureDFFileNames[pos.fst]
  files.rds <- featureDFFileNames[pos.rds]
  featureDFList.fst <- NULL
  featureDFList.rds <- NULL
  if(length(files.fst)>=1){
    message("Reading fst files...")
    featureDFList.fst <- pbapply::pblapply(files.fst, function(f){fst::read.fst(f, as.data.table=T)})
  }
  if(length(files.rds)>=1){
    message("Reading rds files...")
    featureDFList.rds <- pbapply::pblapply(files.rds, function(f){readRDS(f)})
  }
  featureDFList <- c(featureDFList.fst, featureDFList.rds)
  names(featureDFList) <- c(paramSet[pos.fst], paramSet[pos.rds])
  return(featureDFList)
}
#' @export
#' @rdname Features_IO
#' @name Features_IO
saveFeatureDFList <- function(featureDFList, fileNameHeader){
  dir.create(dirname(fileNameHeader), showWarnings=F, recursive=T)
  dataFrameQ <- sapply(featureDFList, is.data.frame)
  if(any(dataFrameQ==F)){
    message("Saving rds files...")
    pbapply::pblapply(
      1:length(featureDFList),
      function(i){saveRDS(featureDFList[[i]], paste0(fileNameHeader, names(featureDFList)[i], ".rds"))}
    )
  }else{
    message("Saving fst files...")
    pbapply::pblapply(
      1:length(featureDFList),
      function(i){fst::write.fst(featureDFList[[i]], paste0(fileNameHeader, names(featureDFList)[i], ".fst"))}
    )
  }
  invisible(NULL)
}
