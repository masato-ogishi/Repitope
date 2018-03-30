#' Compile epitope datasets retrieved from IEDB and other sources.
#'
#' @param df An epitope dataframe.
#' @param IEDBAssayFileName A set of T cell assay results obtained from IEDB.
#' @param OtherFileNames Other sources. Must contain a "Dataset" column indicating the data source.
#' @param IEDBEpitopeFileName A set of epitopes obtained from IEDB.
#' @param IEDBEpitopeFileNames A set of epitope files from IEDB.
#' @importFrom dplyr %>%
#' @importFrom dplyr left_join
#' @importFrom dplyr bind_rows
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom dplyr mutate
#' @importFrom dplyr transmute
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr ungroup
#' @importFrom purrr flatten_chr
#' @importFrom readr read_csv
#' @importFrom magrittr set_colnames
#' @importFrom stringr str_replace
#' @importFrom stringr str_replace_all
#' @importFrom stringr str_detect
#' @importFrom stringr str_split
#' @importFrom scales percent
#' @export
#' @rdname DataPreparation_EpitopeDataset
#' @name DataPreparation_EpitopeDataset
Epitope_Import <- function(IEDBAssayFileName, OtherFileNames=NULL){
  df.epitope.iedb <- readr::read_csv(IEDBAssayFileName, skip=1) %>%
    magrittr::set_colnames(stringr::str_replace_all(colnames(.), " ", "")) %>%
    dplyr::transmute(
      Peptide=Description,
      Host=Name,
      Immunogenicity=ifelse(stringr::str_detect(QualitativeMeasure, "Positive"), "Positive", "Negative"),
      ImmunogenicityEvidence=AssayGroup,
      MHC=AlleleName,
      MHCEvidence=AlleleEvidenceCode,
      Dataset="IEDB"
    ) %>%
    dplyr::filter(ImmunogenicityEvidence %in% readr::read_csv(system.file("IEDB_ImmunogenicityEvidenceTable.csv", package="Repitope"))[[1]])
  if(!is.null(OtherFileNames)){
    df.epitope.others <- dplyr::bind_rows(lapply(OtherFileNames, readr::read_csv))
    df.epitope.iedb <- dplyr::bind_rows(df.epitope.iedb, df.epitope.others)
  }
  df.epitope.iedb <- df.epitope.iedb %>%
    dplyr::filter(nchar(Peptide) %in% 8:11) %>%         ## Length check.
    dplyr::filter(Peptide %in% sequenceFilter(Peptide)) ## Sequences with non-standard characters were discarded.
  return(df.epitope.iedb)
}

#' @export
#' @rdname DataPreparation_EpitopeDataset
#' @name DataPreparation_EpitopeDataset
Epitope_SolveImmunogenicityContradiction <- function(df){
  df <- df %>%
    dplyr::group_by(Peptide) %>%
    dplyr::summarise(
      Immunogenicity=paste0(sort(Immunogenicity), collapse="|"),
      ImmunogenicityEvidence=paste0(sort(unique(ImmunogenicityEvidence)), collapse="|"),
      MHC=stringr::str_replace_all(paste0(sort(unique(MHC)), collapse="|"), ", ", "|"),
      MHCEvidence=paste0(sort(unique(MHCEvidence)), collapse="|"),
      Host=paste0(sort(unique(Host)), collapse="|"),
      Dataset=paste0(sort(unique(Dataset)), collapse="|")
    ) %>%
    dplyr::mutate(Immunogenicity_Contradiction=stringr::str_detect(Immunogenicity, stringr::fixed("Negative|Positive")))

  cont <- table(df$"Immunogenicity_Contradiction")
  cat(cont[[2]], "/", nrow(df), " (", scales::percent(cont[[2]]/nrow(df)), ") peptides have contradicting annotations.\n", sep="")

  chooseOneImmunogenicity_Definitive <- function(ImmunogenicityString){
    paste0(unique(purrr::flatten_chr(stringr::str_split(ImmunogenicityString, stringr::fixed("|")))), collapse="|")
  }
  chooseOneImmunogenicity_Consensus <- function(ImmunogenicityString){
    which.max(table(purrr::flatten_chr(stringr::str_split(ImmunogenicityString, stringr::fixed("|")))))
  }
  df <- df %>%
    dplyr::mutate(
      #Immunogenicity_Definitive=sapply(Immunogenicity, chooseOneImmunogenicity_Definitive)
      Immunogenicity_Loose=ifelse(stringr::str_detect(Immunogenicity, "Positive"), "Positive", "Negative")
      #Immunogenicity_Strict=ifelse(stringr::str_detect(Immunogenicity, "Negative"), "Negative", "Positive")
      #Immunogenicity_Consensus=sapply(Immunogenicity, chooseOneImmunogenicity_Consensus)
    ) %>%
    dplyr::transmute(
      Peptide,
      Immunogenicity=Immunogenicity_Loose,  ## Seems the best
      ImmunogenicityEvidence,
      MHC,
      MHCEvidence,
      Host,
      Dataset
    ) %>%
    dplyr::filter(Immunogenicity %in% c("Positive", "Negative")) %>%
    dplyr::mutate(Immunogenicity=factor(Immunogenicity, levels=c("Positive", "Negative")))

  imm <- table(df$"Immunogenicity")
  cat(imm[[1]], "/", nrow(df), " (", scales::percent(imm[[1]]/nrow(df)), ") peptides are immunogenic in at least one of the observations.\n", sep="")

  return(df)
}

#' @export
#' @rdname DataPreparation_EpitopeDataset
#' @name DataPreparation_EpitopeDataset
Epitope_Add_Disease <- function(df, IEDBEpitopeFileName=system.file("IEDB_Epitope_Healthy.csv.gz", package="Repitope")){
  df.meta <- suppressWarnings(readr::read_csv(IEDBEpitopeFileName, skip=1))
  pept.meta <- df.meta$Description
  df$"Disease" <- "Disease"
  df$"Disease"[which(df$"Peptide" %in% pept.meta)] <- "Healthy"
  return(df)
}

#' @export
#' @rdname DataPreparation_EpitopeDataset
#' @name DataPreparation_EpitopeDataset
Epitope_Add_HLASerotypes <- function(df, IEDBEpitopeFileNames=system.file("IEDB_Epitope_Serotype_HLA-A01.csv.gz", package="Repitope")){
  hlaSerotypeDF <- function(IEDBEpitopeFileName){
    df.meta <- suppressWarnings(readr::read_csv(IEDBEpitopeFileName, skip=1))
    pept.meta <- df.meta$Description
    sero <- grep("HLA-", strsplit(basename(IEDBEpitopeFileName), "_")[[1]], value=T)
    sero <- stringr::str_replace(strsplit(sero, ".", fixed=T)[[1]][[1]], "HLA-", "")
    df.meta <- data.frame("Peptide"=pept.meta, "HLASerotype"=sero, stringsAsFactors=F)
    return(df.meta)
  }
  df.meta <- dplyr::bind_rows(lapply(IEDBEpitopeFileNames, hlaSerotypeDF)) %>%
    dplyr::group_by(Peptide) %>%
    dplyr::summarise(HLASerotype=paste0(sort(unique(HLASerotype)), collapse="|")) %>%
    dplyr::filter(nchar(Peptide) %in% 8:11) %>%
    dplyr::filter(Peptide %in% sequenceFilter(Peptide))
  df <- dplyr::left_join(df, df.meta, by="Peptide")
  return(df)
}
