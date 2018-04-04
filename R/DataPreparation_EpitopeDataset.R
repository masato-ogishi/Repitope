#' Compile epitope datasets retrieved from IEDB and other sources.
#'
#' \code{Epitope_Import} imports IEDB and other epitope files.\cr
#' \code{Epitope_SolveImmunogenicityContradiction} solves contradicted annotations on immunogenicity.\cr
#' \code{Epitope_Add_Disease} adds metadata of being related to diseases or not obtained from the IEDB database.\cr
#' \code{Epitope_Add_HLASerotypes} adds HLA serotype metadata obtained from the IEDB database.\cr
#' \code{compressedToLongFormat} converts a compresed column into a long-format column. The compressed strings should be separated by "|".\cr
#' \code{compressedToDummyDF} converts a compresed column into a dummy variable dataframe. The compressed strings should be separated by "|".\cr
#' \code{Epitope_Convert_HLASupertypes} converts the compressed HLA genotype information into HLA supertypes. The definition of HLA supertypes is derived from the Additional File 1 of Sidney et al., 2008.\cr
#'
#' @param df An epitope dataframe.
#' @param compressedColumnName A string indicating the names of the compressed column to be converted into a long format.
#' @param IEDBAssayFileName A set of T cell assay results obtained from IEDB.
#' @param OtherFileNames Other epitope source files. Must contain a "Dataset" column indicating the data source.
#' @param IEDBEpitopeFileName A set of epitopes obtained from IEDB.
#' @param IEDBEpitopeFileNames A set of epitope files from IEDB.
#' @importFrom dplyr %>%
#' @importFrom dplyr left_join
#' @importFrom dplyr bind_rows
#' @importFrom dplyr bind_cols
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom dplyr mutate
#' @importFrom dplyr transmute
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr ungroup
#' @importFrom dplyr first
#' @importFrom dplyr last
#' @importFrom tidyr drop_na
#' @importFrom tidyr spread
#' @importFrom purrr flatten_chr
#' @importFrom purrr flatten_dfr
#' @importFrom readr read_csv
#' @importFrom magrittr set_colnames
#' @importFrom stringr str_replace
#' @importFrom stringr str_replace_all
#' @importFrom stringr str_detect
#' @importFrom stringr str_split
#' @importFrom stringr fixed
#' @importFrom scales percent
#' @importFrom zoo coredata
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

#' @export
#' @rdname DataPreparation_EpitopeDataset
#' @name DataPreparation_EpitopeDataset
compressedToLongFormat <- function(df, compressedColumnName){
  cmp <- stringr::str_split(df[[compressedColumnName]], stringr::fixed("|"))
  l <- sapply(cmp, length)
  df_long <- lapply(1:nrow(df), function(i){replicate(l[i], zoo::coredata(df[i,]), simplify=F)})
  df_long <- purrr::flatten_dfr(df_long)
  df_long[[compressedColumnName]] <- unlist(cmp)
  return(df_long)
}

#' @export
#' @rdname DataPreparation_EpitopeDataset
#' @name DataPreparation_EpitopeDataset
compressedToDummyDF <- function(df, compressedColumnName){
  df_dummy <- compressedToLongFormat(df, compressedColumnName) %>%
    dplyr::mutate(Value=1) %>%
    tidyr::drop_na(compressedColumnName) %>%
    tidyr::spread(compressedColumnName, Value, fill=0)
  dplyr::left_join(df, df_dummy)
}

#' @export
#' @rdname DataPreparation_EpitopeDataset
#' @name DataPreparation_EpitopeDataset
Epitope_Convert_HLASupertypes <- function(df, compressedColumnName="MHC"){
  # Stanadrdize HLA strings
  HLAStrings <- df[[compressedColumnName]]
  hla <- HLAStrings %>%
    stringr::str_replace_all("HLA-", "") %>%
    stringr::str_replace_all("HLA ", "") %>%
    stringr::str_replace_all(":", "") %>%
    stringr::str_split(stringr::fixed("|"))
  hlaStandardize <- function(HLAStrings){
    hla.a <- stringr::str_detect(HLAStrings, "^A.+")
    hla.b <- stringr::str_detect(HLAStrings, "^B.+")
    hla <- hla.a | hla.b
    HLAStrings[which(hla==F)] <- NA
    HLAStrings <- unique(unlist(lapply(stringr::str_split(HLAStrings, " "), dplyr::first)))
    HLAStrings <- stringr::str_replace(HLAStrings, stringr::fixed("*"), "")
    HLAStrings <- stringr::str_replace(HLAStrings, "^A1$", "A01")
    HLAStrings <- stringr::str_replace(HLAStrings, "^A2$", "A02")
    HLAStrings <- stringr::str_replace(HLAStrings, "^A3$", "A03")
    HLAStrings <- stringr::str_replace(HLAStrings, "^B7$", "B07")
    HLAStrings <- stringr::str_replace(HLAStrings, "^B8$", "B08")
    return(HLAStrings)
  }
  cat("Standardizing HLA genotype strings...\n")
  hla <- pbapply::pblapply(hla, hlaStandardize)

  # Match supertypes
  HLA_ST <- suppressMessages(readr::read_csv(system.file("HLASupertypeConversion.csv", package="Repitope")))  ## Sidney et al., 2008. Additional File 1.
  HLA_ST <- HLA_ST[1:2]
  colnames(HLA_ST) <- c("Allele", "Supertype")
  HLA_ST$Supertype[which(HLA_ST$Supertype=="Unclassified")] <- NA
  HLA_ST$Supertype <- stringr::str_replace_all(HLA_ST$Supertype, " ", "|")
  HLA_ST$Allele <- stringr::str_replace_all(HLA_ST$Allele, stringr::fixed("*"), "")
  hlaSupertype_CompleteMatch <- function(HLAStrings){
    s <- paste0(paste0("^", HLAStrings, "$"), collapse="|")
    s <- dplyr::filter(HLA_ST, stringr::str_detect(Allele, s))$Supertype
    s <- sort(unique(unlist(stringr::str_split(s, stringr::fixed("|")))))
    return(s)
  }
  hlaSupertype_TwoDigitMatch <- function(HLAStrings){
    s <- paste0(paste0("^", HLAStrings, ".+"), collapse="|")
    s <- dplyr::filter(HLA_ST, stringr::str_detect(Allele, s))$Supertype
    s <- unlist(stringr::str_split(s, stringr::fixed("|")))
    s <- dplyr::last(names(sort(table(s))))
    return(s)
  }
  cat("Matching HLA genotypes to supertypes...\n")
  hla_comp <- pbapply::pblapply(hla, hlaSupertype_CompleteMatch)
  hla_tg <- pbapply::pblapply(hla, hlaSupertype_TwoDigitMatch)
  hla <- mapply(c, hla_comp, hla_tg, SIMPLIFY=F)
  hla <- lapply(hla, unique)
  hla[which(sapply(hla, length)==0)] <- NA
  hla <- unlist(lapply(hla, paste0, collapse="|"))
  return(dplyr::mutate(df, HLASupertype=hla))
}
