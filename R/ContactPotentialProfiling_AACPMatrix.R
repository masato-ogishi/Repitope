#' Amino acid pairwise contact potential (AACP) matrix.
#'
#' \code{CPP_AACPMatrix} generates a formatted AACP datatable.
#'
#' @param sourceFile A source file. By default, the AACP dataset downloaded from the official AAIndex database will be used. Alternatively, users can provide modified AACP scales.
#' @export
#' @rdname ContactPotentialProfiling_AACPMatrix
#' @name ContactPotentialProfiling_AACPMatrix
CPP_AACPMatrix <- function(sourceFile=system.file("AACPMatrix.csv", package="Repitope")){
  AACP_DF <- suppressMessages(readr::read_csv(sourceFile, na=c(""))) %>%
    dplyr::select(-Description, -PMID, -Author, -Title, -Reference) %>%
    tidyr::gather(AAPair, Value, -AAIndexID) %>%
    tidyr::spread(AAIndexID, Value)
  AACP_DF_Rev <- AACP_DF %>%
    dplyr::mutate("AAPair"=Biostrings::reverse(AACP_DF$"AAPair"))
  AACP_DF <- dplyr::bind_rows(AACP_DF, AACP_DF_Rev) %>% dplyr::distinct(AAPair, .keep_all=T)
  AACP_DF <- AACP_DF %>%
    tidyr::gather(AAIndexID, Value, -AAPair) %>%
    dplyr::mutate(AA.1=stringr::str_sub(AAPair, 1, 1), AA.2=stringr::str_sub(AAPair, 2, 2)) %>%
    dplyr::select(-AAPair) %>%
    tidyr::spread(AA.2, Value) %>%
    dplyr::select(-AA.1)
  formatting <- function(aaIndexID, pairMatInverse=T){
    pairMat <- AACP_DF %>%
      dplyr::filter(AAIndexID==aaIndexID) %>%
      dplyr::select(-AAIndexID)
    pairMat <- as.matrix(pairMat)
    rownames(pairMat) <- colnames(pairMat)
    if(pairMatInverse){ pairMat <- -pairMat }
    pairMat <- scales::rescale(pairMat)
    return(pairMat)
  }
  aaIndexIDSet <- unique(c(AACP_DF$AAIndexID, paste0(AACP_DF$AAIndexID, "inv")))
  AACP_DT <- c(lapply(unique(AACP_DF$AAIndexID), function(a){formatting(a, pairMatInverse=F)}),
               lapply(unique(AACP_DF$AAIndexID), function(a){formatting(a, pairMatInverse=T)}))
  names(AACP_DT) <- aaIndexIDSet
  AACP_DT <- data.table::rbindlist(lapply(AACP_DT, as.data.frame))
  AACP_DT[["AAIndexID"]] <- unlist(lapply(aaIndexIDSet, function(id){rep(id, length(Biostrings::AA_STANDARD))}))
  AACPMatrixList <- split(AACP_DT, by="AAIndexID", keep.by=F)
  AACPMatrixList <- lapply(AACPMatrixList, as.matrix)
  AACPMatrixList <- lapply(AACPMatrixList, magrittr::set_rownames, sort(Biostrings::AA_STANDARD))
  return(AACPMatrixList)
}
