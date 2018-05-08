#' Compile epitope datasets retrieved from IEDB and other sources.
#'
#' \code{Epitope_Import} imports IEDB and other epitope files. Contradicting annotations on immunogenicity are internally resolved. \cr
#' \code{compressedToLongFormat} converts a compresed column into a long-format column. The compressed strings should be separated by "|".\cr
#' \code{compressedToDummyDF} converts a compresed column into a dummy variable dataframe. The compressed strings should be separated by "|".\cr
#'
#' @param IEDBAssayFileName A set of T cell assay results obtained from IEDB.
#' @param OtherFileNames Other epitope source files. Must contain a "Dataset" column indicating the data source.
#' @param peptideLengthSet A set of peptide lengths to be incorporated.
#' @param df A dataframe/datatable of epitope dataset.
#' @param compressedColumnName A string indicating the names of the compressed column to be converted into a long format.
#' @importFrom dplyr %>%
#' @export
#' @rdname DataPreparation_EpitopeDataset
#' @name DataPreparation_EpitopeDataset
Epitope_Import <- function(IEDBAssayFileName=NULL, OtherFileNames=NULL, peptideLengthSet=8:11){
  if(!is.null(IEDBAssayFileName)){
    immEvi <- suppressMessages(readr::read_csv(system.file("IEDB_ImmunogenicityEvidenceTable.csv", package="Repitope"))[[1]])
    df <- suppressWarnings(suppressMessages(readr::read_csv(IEDBAssayFileName, skip=1))) %>%
      magrittr::set_colnames(stringr::str_replace_all(colnames(.), " ", "")) %>%
      dplyr::transmute(
        Peptide=Description,
        Immunogenicity=factor(ifelse(stringr::str_detect(QualitativeMeasure, "Positive"), "Positive", "Negative"), levels=c("Positive", "Negative")),
        ImmunogenicityEvidence=AssayGroup,
        MHC=AlleleName,
        MHCEvidence=AlleleEvidenceCode,
        Host=Name,
        Organism=ImmunogenOrganismSpecies,
        Dataset="IEDB"
      ) %>%
      dplyr::filter(ImmunogenicityEvidence %in% immEvi)
  }else{
    df <- NULL
  }
  if(!is.null(OtherFileNames)){
    df.others <- dplyr::bind_rows(lapply(OtherFileNames, function(f){suppressWarnings(suppressMessages(readr::read_csv(f)))})) %>%
      dplyr::transmute(
        Peptide,
        Immunogenicity,
        ImmunogenicityEvidence=NA,
        MHC,
        MHCEvidence=NA,
        Host=NA,
        Organism=NA,
        Dataset
      )
  }else{
    df.others <- NULL
  }
  df <- suppressWarnings(suppressMessages(dplyr::bind_rows(df, df.others))) %>%
    dplyr::filter(nchar(Peptide) %in% peptideLengthSet) %>%   ## Length check.
    dplyr::filter(Peptide %in% sequenceFilter(Peptide)) %>%   ## Sequences with non-standard characters were discarded.
    dplyr::mutate(Immunogenicity=factor(Immunogenicity, levels=c("Positive", "Negative"))) %>%
    dplyr::group_by(Peptide) %>%
    dplyr::summarise(
      Immunogenicity=paste0(sort(Immunogenicity), collapse="|"),
      ImmunogenicityEvidence=paste0(sort(unique(ImmunogenicityEvidence)), collapse="|"),
      MHC=stringr::str_replace_all(paste0(sort(unique(MHC)), collapse="|"), ", ", "|"),
      MHCEvidence=paste0(sort(unique(MHCEvidence)), collapse="|"),
      Host=paste0(sort(unique(Host)), collapse="|"),
      Organism=paste0(sort(unique(Organism)), collapse="|"),
      Dataset=paste0(sort(unique(Dataset)), collapse="|")
    )

  resolveImmunogenicityContradiction <- function(df){
    df <- df %>%
      dplyr::mutate(Immunogenicity_Contradiction=stringr::str_detect(Immunogenicity, stringr::fixed("Positive|Negative")))

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
        Immunogenicity_Loose=dplyr::if_else(stringr::str_detect(Immunogenicity, "Positive"), "Positive", "Negative")
        #Immunogenicity_Strict=dplyr::if_else(stringr::str_detect(Immunogenicity, "Negative"), "Negative", "Positive")
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
  df <- resolveImmunogenicityContradiction(df)

  return(data.table::as.data.table(df))
}

#' @export
#' @rdname DataPreparation_EpitopeDataset
#' @name DataPreparation_EpitopeDataset
compressedToLongFormat <- function(df, compressedColumnName){
  cmp <- stringr::str_split(df[[compressedColumnName]], stringr::fixed("|"))
  l <- sapply(cmp, length)
  df_long <- lapply(1:nrow(df), function(i){replicate(l[i], zoo::coredata(df[i,]), simplify=F)})
  df_long <- purrr::flatten_dfr(df_long)
  df_long[[compressedColumnName]] <- unlist(cmp)
  return(data.table::as.data.table(df_long))
}

#' @export
#' @rdname DataPreparation_EpitopeDataset
#' @name DataPreparation_EpitopeDataset
compressedToDummyDF <- function(df, compressedColumnName){
  df_dummy <- compressedToLongFormat(df, compressedColumnName) %>%
    dplyr::mutate(Value=1) %>%
    tidyr::drop_na(compressedColumnName) %>%
    tidyr::spread(compressedColumnName, Value, fill=0)
  left_join_0 <- function(x, y, fill=0){
    z <- dplyr::left_join(x, y)
    tmp <- setdiff(names(z), names(x))
    z <- tidyr::replace_na(z, setNames(as.list(rep(fill, length(tmp))), tmp))
    return(z)
  }
  df_dummy <- left_join_0(df, df_dummy, fill=0)
  return(data.table::as.data.table(df_dummy))
}

