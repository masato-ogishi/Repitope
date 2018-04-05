#' Utilities for NetMHC.
#'
#' @param peptideSet A set of peptide sequences.
#' @param outDir A path to the output directory.
#' @param netMHCFileNames A set of tab-delimited NetMHC result files.
#' @param hlaSuperTypeSet A set of HLA supertypes for which binding was predicted.
#' @importFrom purrr flatten
#' @importFrom seqinr write.fasta
#' @importFrom BBmisc chunk
#' @export
#' @rdname DataPreparation_NetMHC
#' @name DataPreparation_NetMHC
NetMHC_PeptideExport <- function(peptideSet, outDir="./NetMHC/"){
  peptList <- split(peptideSet, nchar(peptideSet))
  peptList <- lapply(peptList, BBmisc::chunk, chunk.size=5000)
  lengList <- names(peptList)
  nchunkList <- sapply(peptList, length)
  peptList <- purrr::flatten(peptList)
  names(peptList) <- unlist(mapply(function(x, y){paste0(x, "mer_", 1:y)}, lengList, nchunkList, SIMPLIFY=F))
  invisible(mapply(
    function(peptSeq, name){seqinr::write.fasta(as.list(peptSeq), peptSeq, paste0(outDir, name, ".fasta"))},
    peptList, names(peptList), SIMPLIFY=F
  ))
}

#' @export
#' @rdname DataPreparation_NetMHC
#' @name DataPreparation_NetMHC
NetMHC_Import <- function(netMHCFileNames, hlaSuperTypeSet=c("A01","A02","A03","A24","A26","B07","B08","B27","B39","B44","B58","B62")){
  binderCategory <- function(rankPercent){
    if(rankPercent<=0.5) return("StrongBinder")
    if(rankPercent<=2) return("WeakBinder")
    return("NonBinder")
  }
  df.netmhc.all <- data.table::rbindlist(lapply(
    netMHCFileNames,
    function(netMHCFileName){suppressWarnings(suppressMessages(readr::read_delim(netMHCFileName, "\t", escape_double=F, trim_ws=T, skip=1)))}
  )) %>%
    dplyr::select(Peptide, dplyr::matches("Rank"), -H_Avg_Ranks) %>%
    magrittr::set_colnames(c("Peptide", hlaSuperTypeSet)) %>%
    tidyr::gather(HLARestriction, Rank, -Peptide) %>%
    dplyr::mutate(BinderCategory=factor(sapply(Rank, binderCategory), levels=c("NonBinder", "WeakBinder", "StrongBinder")))
  return(df.netmhc.all)
}

