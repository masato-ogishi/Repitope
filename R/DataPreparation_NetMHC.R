#' Utility functions for NetMHC.
#'
#' @param peptideSet A set of peptide sequences.
#' @param outDir A path to the output directory.
#' @param outFileName A file name for an output.
#' @param netMHCFastaFileNames A set of FASTA files for NetMHC analysis.
#' @param netMHCOutputFileNames A set of tab-delimited NetMHC result files.
#' @importFrom purrr flatten
#' @importFrom seqinr write.fasta
#' @importFrom BBmisc chunk
#' @importFrom readr write_delim
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
NetMHC_Import <- function(netMHCOutputFileNames){
  hlaSuperTypeSet <- c("A01","A02","A03","A24","A26","B07","B08","B27","B39","B44","B58","B62")
  binderCategory <- function(rankPercent){
    if(rankPercent<=0.5) return("StrongBinder")
    if(rankPercent<=2) return("WeakBinder")
    return("NonBinder")
  }
  df.netmhc.all <- data.table::rbindlist(lapply(
    netMHCOutputFileNames,
    function(netMHCFileName){suppressWarnings(suppressMessages(readr::read_delim(netMHCFileName, "\t", escape_double=F, trim_ws=T, skip=1)))}
  )) %>%
    dplyr::select(Peptide, dplyr::matches("Rank"), -H_Avg_Ranks) %>%
    magrittr::set_colnames(c("Peptide", hlaSuperTypeSet)) %>%
    tidyr::gather(HLARestriction, Rank, -Peptide) %>%
    dplyr::mutate(BinderCategory=factor(sapply(Rank, binderCategory), levels=c("NonBinder", "WeakBinder", "StrongBinder")))
  return(df.netmhc.all)
}

#' @export
#' @rdname DataPreparation_NetMHC
#' @name DataPreparation_NetMHC
NetMHC_Script <- function(netMHCFastaFileNames, outDir="./NetMHC/", outFileName="NetMHC_Script.txt"){
  NetMHC_Script_Single <- function(netMHCFastaFileName){
    l <- as.integer(gsub("mer", "", strsplit(netMHCFastaFileName, "_")[[1]][[3]]))
    o <- gsub(".fasta", ".xls", netMHCFastaFileName)
    paste0("../netMHC -a HLA-A0101,HLA-A0201,HLA-A0301,HLA-A2402,HLA-A2601,HLA-B0702,HLA-B0801,HLA-B2705,HLA-B3901,HLA-B4001,HLA-B5801,HLA-B1501 -f ", netMHCFastaFileName, " -p -l ", l, " -xls -xlsfile ", o)
  }
  readr::write_delim(data.frame("Script"=lapply(netMHCFastaFileNames, NetMHC_Script_Single)), path=file.path(outDir, outFileName), delim="\n", col_names=F)
}
