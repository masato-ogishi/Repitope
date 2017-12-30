#' Feature dataframes In/Out.
#'
#' @param featureDFFileNameList A set of binary file names with ".fst" extensions.
#' @param featureDFList A named list of feature dataframes.
#' @param fileNameHeader A file name header.
#' @importFrom fst read.fst
#' @importFrom fst write.fst
#' @importFrom pbapply pblapply
#' @export
#' @rdname Features_IO
#' @name Features_IO
readFeatureDFList <- function(featureDFFileNameList){
  paramSet <- gsub(".fst$", "", basename(featureDFFileNameList))
  paramSet <- as.character(t(as.data.frame(strsplit(paramSet, "_", fixed=T), stringsAsFactors=F))[,2])
  featureDFList <- pbapply::pblapply(featureDFFileNameList, fst::read.fst, as.data.table=T)
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
      fst::write.fst(featureDFList[[i]], paste0(fileNameHeader, names(featureDFList)[i], ".fst"))
    }
  )
}
