#' Feature dataframes In/Out.
#'
#' @param featureDFFileNames A character vector of binary file names with either ".fst" or ".rds" extensions.
#' @param featureDFList A named list of feature dataframes.
#' @param fileNameHeader A file name header.
#' @importFrom fst read.fst
#' @importFrom fst write.fst
#' @importFrom pbapply pblapply
#' @export
#' @rdname Features_IO
#' @name Features_IO
readFeatureDFList <- function(featureDFFileNames){
  paramSet <- gsub(".rds$", "", gsub(".fst$", "", basename(featureDFFileNames)))
  paramSet <- as.character(t(as.data.frame(strsplit(paramSet, "_", fixed=T), stringsAsFactors=F))[,2])
  featureDFList <- pbapply::pblapply(
    featureDFFileNames,
    function(f){
      if(grep(".fst$", basename(f))==1){
        fst::read.fst(f, as.data.table=T)
      }else if(grep(".rds$", basename(f))==1){
        readRDS(f)
      }else{
        NULL
      }
    }
  )
  names(featureDFList) <- paramSet
  return(featureDFList)
}
#' @export
#' @rdname Features_IO
#' @name Features_IO
saveFeatureDFList <- function(featureDFList, fileNameHeader){
  pbapply::pblapply(
    1:length(featureDFList),
    function(i){
      if(is.data.frame(featureDFList[[i]])){
        fst::write.fst(featureDFList[[i]], paste0(fileNameHeader, names(featureDFList)[i], ".fst"))
      }else{
        saveRDS(featureDFList[[i]], paste0(fileNameHeader, names(featureDFList)[i], ".rds"))
      }
    }
  )
  invisible(return(NULL))
}
