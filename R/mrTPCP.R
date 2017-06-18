# Perform modified repertoire-wide TCR-peptide contact profile (mrTPCP) analysis
#' Calculate modified repertoire-wide TCR-peptide contact profile (mrTPCP) variables with a fixed window size
#' @param peptList A character vector containing input peptide sequences.
#' @param rept A character vector containing input TCR CDR3 sequences. A default is NULL; in this case, a reference CD8 TCR repertoire used in the original paper will be used. Alternatively, users can provide a character vector containing a CDR3 sequence repertoire of interest.
#' @param wind A size of the window.
#' @param aa.index.id A character vector of AAIndex names to be used. Can also be NULL or "all"; in those cases, all AAIndex scales will be used.
#' @importFrom readr read_csv
#' @importFrom dplyr %>%
#' @importFrom dplyr filter
#' @importFrom stringr str_sub
#' @importFrom dplyr mutate
#' @importFrom stringr str_split
#' @importFrom stringr str_detect
#' @importFrom dplyr bind_rows
#' @importFrom dplyr bind_cols
#' @importFrom tidyr drop_na
#' @importFrom dplyr as_data_frame
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise_at
#' @importFrom dplyr ungroup
#' @importFrom dplyr select
#' @importFrom tidyr spread
#' @importFrom tidyr gather
#' @importFrom magrittr set_colnames
mrTPCP.fixedWindow <- function(peptList, rept=NULL, wind=4, aa.index.id=c("MIYS990106")){
  # Functions for string fragmentation
  slidingFragment <- function(seq, wind){
    aa.list <- stringr::str_split(seq, "") %>% unlist()
    len <- nchar(seq)
    if(len < wind) return(NA)
    sapply(1:(len-wind+1), function(i){paste0(aa.list[rep(i:len, length=wind)], collapse="")})
  }
  slidingFragment.Batch <- function(seqList, wind){
    lapply(seqList, function(s){slidingFragment(s, wind)})
  }

  # Import AACP AAIndex matrix
  AACP.AAIndex.DF <- readr::read_csv(system.file("AACP_AAIndex_Matrix.csv", package="Repitope"), na="")
  if(is.null(aa.index.id)||aa.index.id=="all"){
    AACP.AAIndex.DF <- AACP.AAIndex.Matrix.Format(AACP.AAIndex.DF)
  } else {
    AACP.AAIndex.DF <- AACP.AAIndex.DF %>%
      dplyr::filter(AAIndexID %in% aa.index.id) %>%
      AACP.AAIndex.Matrix.Format()
  }
  AAPairsList <- colnames(AACP.AAIndex.DF)
  AACPValues <- setNames(as.vector(as.matrix(AACP.AAIndex.DF)), AAPairsList)

  # Compile TCR repertoire into a set of fragments
  if(is.null(rept)){
    rept.seq <- read.delim(system.file("TCRRepertoire_CD8.txt", package="Repitope"), header=T, sep="\t")[["cdr3aa"]] %>%
      as.character() %>% toupper()
  } else {
    rept.seq <- rept %>% as.character() %>% toupper()
  }
  rept.frag <- unique(unlist(slidingFragment.Batch(rept.seq[nchar(rept.seq)>wind], wind)))

  # Length check
  if(length(unique(nchar(peptList)))>=2){
    print("The lengths of the input peptides are not the same!")
    return(NULL)
  }
  pept.length <- nchar(peptList[[1]])
  n.windows <- pept.length - wind + 1
  if(n.windows<=0){
    print("The window size is larger than the peptide length!")
    return(NULL)
  }

  # Calculate AACP values at the level of peptide fragments
  AACPByWindow.Fragment <- function(pept.frag, window.id){
    aaPairsMat <- mapply(function(x,y){paste0(x,y)},
                         x=strsplit(pept.frag, ""),
                         y=strsplit(rept.frag, ""))
    aacpMat <- matrix(AACPValues[aaPairsMat], nrow=wind)
    aacpSums <- colSums(aacpMat) # summed for each TCR fragment
    TCRFragmentIDs <- paste0("TCRFragment_", formatC(1:length(rept.frag), flag="0", width=(floor(log10(length(rept.frag)))+1)))
    aacpSumsDF <- data.frame("AACPValue"=aacpSums, "TCRFragmentID"=TCRFragmentIDs, "WindowID"=window.id)
    return(aacpSumsDF)
  }

  # Calculate AACP values at the level of peptides
  AACPByWindow.Peptide <- function(pept){
    Q10 <- function(v){quantile(v, probs=0.1, names=F)}
    Q20 <- function(v){quantile(v, probs=0.2, names=F)}
    Q30 <- function(v){quantile(v, probs=0.3, names=F)}
    Q40 <- function(v){quantile(v, probs=0.4, names=F)}
    Q60 <- function(v){quantile(v, probs=0.6, names=F)}
    Q70 <- function(v){quantile(v, probs=0.7, names=F)}
    Q80 <- function(v){quantile(v, probs=0.8, names=F)}
    Q90 <- function(v){quantile(v, probs=0.9, names=F)}
    mrtpcp_df <- suppressWarnings(dplyr::bind_rows(
      mapply(AACPByWindow.Fragment, slidingFragment(pept, wind), paste0("Wind", 1:n.windows), SIMPLIFY=F)
    ))
    mrtpcp_df.ByTCRFragment <- mrtpcp_df %>%
      dplyr::group_by(TCRFragmentID) %>%
      dplyr::summarise_at(vars(AACPValue), funs(min, max, median)) %>%
      magrittr::set_colnames(c("TCRFragmentID", paste0("AACPByTCRFragment_", c("min", "max", "median")))) %>%
      dplyr::ungroup() %>%
      tidyr::gather(Variable, AACPValue, -TCRFragmentID) %>%
      dplyr::group_by(Variable) %>%
      dplyr::summarise_at(vars(AACPValue), funs(min, max, median, Q10, Q20, Q30, Q40, Q60, Q70, Q80, Q90)) %>%
      tidyr::gather(Stat, AACPValue, -Variable) %>%
      dplyr::mutate(Variable=paste0(Variable, "_", Stat)) %>%
      dplyr::select(Variable, AACPValue) %>%
      tidyr::spread(Variable, AACPValue)
    mrtpcp_df.ByWindow <- mrtpcp_df %>%
      dplyr::group_by(WindowID) %>%
      dplyr::summarise_at(vars(AACPValue), funs(min, max, median, Q10, Q20, Q30, Q40, Q60, Q70, Q80, Q90)) %>%
      tidyr::gather(Variable, AACPValue, -WindowID) %>%
      dplyr::mutate(Variable=paste0(WindowID, "_", Variable)) %>%
      dplyr::select(Variable, AACPValue) %>%
      tidyr::spread(Variable, AACPValue)
    mrtpcp_df.ByWindow <- mrtpcp_df.ByWindow %>% magrittr::set_colnames(paste0("AACPByWindow_", colnames(mrtpcp_df.ByWindow)))
    mrtpcp_df <- dplyr::bind_cols(mrtpcp_df.ByTCRFragment, mrtpcp_df.ByWindow)
    mrtpcp_df <- mrtpcp_df %>% magrittr::set_colnames(paste0(colnames(mrtpcp_df), "_", wind))
    return(mrtpcp_df)
  }
  mrtpcp_df <- suppressWarnings(dplyr::bind_rows(lapply(peptList, AACPByWindow.Peptide)))
  return(mrtpcp_df)
}

#' Calculate modified repertoire-wide TCR-peptide contact profile (mrTPCP) variables
#' @param peptList A character vector containing input peptide sequences.
#' @param repertoire A character vector containing input TCR CDR3 sequences. A default is NULL; in this case, a reference CD8 TCR repertoire used in the original paper will be used. Alternatively, users can provide a character vector containing a CDR3 sequence repertoire of interest.
#' @param winds A vector of the sizes of window.
#' @param aa.index.id A character vector of AAIndex names to be used. Can also be NULL or "all"; in those cases, all AAIndex scales will be used.
#' @importFrom dplyr bind_cols
#' @importFrom dplyr as_data_frame
#' @export
mrTPCP <- function(peptList, repertoire=NULL, winds=4:5, aa.index.id=c("MIYS990106")){
  # Length check
  peptList_th <- peptList[nchar(peptList)>=max(winds)]

  # A batch mode
  mrtpcp_df <- data.frame(
    "Peptide"=peptList_th,
    dplyr::bind_cols(
      lapply(winds, function(windowSize){mrTPCP.fixedWindow(peptList_th, repertoire, windowSize, aa.index.id)})
    )
  ) %>% dplyr::as_data_frame()
  return(mrtpcp_df)
}
