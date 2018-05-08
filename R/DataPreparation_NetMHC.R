#' Utility functions for NetMHC.
#'
#' @param peptideSet A set of peptide sequences.
#' @param outDir A path to the output directory.
#' @param netMHCOutputFileNames A set of tab-delimited NetMHC result files.
#' @param MHCType A character string. "MHCI" or "MHCII".
#' @export
#' @rdname DataPreparation_NetMHC
#' @name DataPreparation_NetMHC
NetMHC_PeptideExport <- function(peptideSet, outDir="./NetMHC/"){
  dir.create(outDir, showWarnings=F, recursive=T)
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
NetMHC_Import <- function(netMHCOutputFileNames, MHCType=c("MHCI", "MHCII")){
  if(MHCType=="MHCI"){
    hlaSuperTypeSet <- c("A01","A02","A03","A24","A26","B07","B08","B27","B39","B44","B58","B62")
    binderCategory <- function(rankPercent){
      if(rankPercent<=0.5) return("StrongBinder")
      if(rankPercent<=2) return("WeakBinder")
      return("NonBinder")
    }
    df.netmhc <- data.table::rbindlist(lapply(
      netMHCOutputFileNames,
      function(netMHCFileName){suppressWarnings(suppressMessages(readr::read_delim(netMHCFileName, "\t", escape_double=F, trim_ws=T, skip=1)))}
    )) %>%
      dplyr::select(Peptide, dplyr::matches("Rank"), -H_Avg_Ranks) %>%
      magrittr::set_colnames(c("Peptide", hlaSuperTypeSet)) %>%
      tidyr::gather(HLARestriction, Rank, -Peptide) %>%
      dplyr::mutate(BinderCategory=factor(sapply(Rank, binderCategory), levels=c("NonBinder", "WeakBinder", "StrongBinder")))
  }
  if(MHCType=="MHCII"){
    hlaSuperTypeSet <- c("DRB1","DRB3","DRB4","DRB5","DP","DQ")
    binderCategory <- function(rankPercent){
      if(rankPercent<=2) return("StrongBinder")
      if(rankPercent<=10) return("WeakBinder")
      return("NonBinder")
    }
    df.netmhc <- data.table::rbindlist(lapply(
      netMHCOutputFileNames,
      function(netMHCFileName){suppressWarnings(suppressMessages(readr::read_delim(netMHCFileName, "\t", escape_double=F, trim_ws=T, skip=1)))}
    )) %>%
      dplyr::select(Peptide, dplyr::matches("Rank")) %>%
      magrittr::set_colnames(c("Peptide", hlaSuperTypeSet)) %>%
      tidyr::gather(HLARestriction, Rank, -Peptide) %>%
      dplyr::mutate(BinderCategory=factor(sapply(Rank, binderCategory), levels=c("NonBinder", "WeakBinder", "StrongBinder")))
  }
  df.netmhc.rank <- df.netmhc %>% dplyr::select(-BinderCategory) %>% tidyr::spread(HLARestriction, Rank)
  df.netmhc.bind <- df.netmhc %>% dplyr::mutate(BinderCategory=as.numeric(BinderCategory)) %>% dplyr::select(-Rank) %>% tidyr::spread(HLARestriction, BinderCategory)
  list("DF"=df.netmhc, "DF_Rank"=df.netmhc.rank, "DF_Bind"=df.netmhc.bind)
}

