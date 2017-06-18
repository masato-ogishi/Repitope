# Perform repertoire-wide TCR-peptide contact profile (rTPCP) analysis
#' Calculate repertoire-wide TCR-peptide contact profile (rTPCP) variables with a fixed window size
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
rTPCP.fixedWindow <- function(peptList, rept=NULL, wind=4, aa.index.id=c("MIYS990106")){
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

  # Import TCR repertoire
  if(is.null(rept)){
    rept.seq <- read.delim(system.file("TCRRepertoire_CD8.txt", package="Repitope"), header=T, sep="\t")[["cdr3aa"]] %>%
      as.character() %>% toupper()
  } else {
    rept.seq <- rept %>% as.character() %>% toupper()
  }

  # Length check
  peptList_th <- peptList[nchar(peptList)>=wind]
  rept_th <- rept.seq[nchar(rept.seq)>=wind]

  # Generate fragment pairs
  peptList.id <- paste0("Pept", formatC(1:length(peptList_th), flag="0", width=(floor(log10(length(peptList_th)))+1)))
  rept.id <- paste0("Rept", formatC(1:length(rept_th), flag="0", width=(floor(log10(length(rept_th)))+1)))
  l_max <- max(nchar(peptList_th), nchar(rept_th))
  peptList.2 <- substr(paste0(peptList_th, paste0(rep("X", l_max), collapse="")), 1, l_max)
  rept.2 <- substr(paste0(rept_th, paste0(rep("X", l_max), collapse="")), 1, l_max)
  peptList.3 <- t(sapply(seq(1, l_max-wind+1), function(i){stringr::str_sub(peptList.2, i, i+wind-1)}))
  rept.3 <- t(sapply(seq(1, l_max-wind+1), function(i){stringr::str_sub(rept.2, i, i+wind-1)}))

  fragpairs_df <- expand.grid(peptList.3, rept.3) %>%
    dplyr::mutate(PeptideID=rep(unlist(lapply(peptList.id, function(x){rep(x, l_max-wind+1)})), length(rept_th)*(l_max-wind+1))) %>%
    dplyr::mutate(RepertoireID=unlist(lapply(rept.id, function(x){rep(x, length(peptList_th)*(l_max-wind+1)*(l_max-wind+1))})))

  pairs_mat <- t(mapply(paste0, stringr::str_split(fragpairs_df$"Var1",""), stringr::str_split(fragpairs_df$"Var2","")))
  pairs_mat[stringr::str_detect(pairs_mat, "X")] <- NA

  # Calculate summed contact potentials for each of the fragment pairs
  fragpairs_df <- lapply(
    as.list(as.data.frame(t(AACP.AAIndex.DF))),
    function(aacpValues){
      m <- pairs_mat
      for(i in 1:length(AAPairsList)){
        m[m==AAPairsList[[i]]] <- aacpValues[[i]]
      }
      m <- matrix(as.numeric(m), ncol=wind)
      return(data.frame(fragpairs_df, m, "Fragment.AACP.Sum"=rowSums(m)))
    }
  )
  fragpairs_df <- data.frame(
    dplyr::bind_rows(fragpairs_df),
    "AAIndexID"=as.vector(
      mapply(rep, names(fragpairs_df), sapply(fragpairs_df, nrow))
    )
  ) %>% tidyr::drop_na()

  # Calculate representative statistics
  rangeDiff <- function(v){diff(range(v))}
  secondMax <- function(v){sort(v, partial=2)[length(v)-1]}
  secondMin <- function(v){sort(v, partial=2)[2]}
  deltaMax <- function(v){diff(v)[length(v)-1]}
  deltaMin <- function(v){diff(v)[1]}

  rtpcp_df <- fragpairs_df %>%
    dplyr::group_by(PeptideID, RepertoireID, AAIndexID) %>%
    dplyr::summarise_at(vars(Fragment.AACP.Sum), funs(mean, median, max, secondMax, deltaMax, min, secondMin, deltaMin, rangeDiff)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(PeptideID, AAIndexID) %>%
    dplyr::summarise_at(vars(mean, median, max, secondMax, deltaMax, min, secondMin, deltaMin, rangeDiff), funs(mean, median, max, min)) %>%
    dplyr::ungroup() %>%
    tidyr::gather(Variable, Value, -PeptideID, -AAIndexID) %>%
    dplyr::mutate(Variable=paste0(AAIndexID,"_",Variable)) %>%
    dplyr::select(-AAIndexID) %>%
    tidyr::spread(Variable, Value) %>%
    dplyr::select(-PeptideID)
  colnames(rtpcp_df) <- paste0(colnames(rtpcp_df), "_", wind)

  return(rtpcp_df)
}

#' Calculate repertoire-wide TCR-peptide contact profile (rTPCP) variables
#' @param peptList A character vector containing input peptide sequences
#' @param repertoire A character vector containing input TCR CDR3 sequences
#' @param winds A vector of the sizes of window
#' @param aa.index.id A character vector of AAIndex names to be used. Can also be NULL or "all", and all AAIndex scales will be used.
#' @importFrom dplyr bind_cols
#' @importFrom dplyr as_data_frame
#' @export
rTPCP <- function(peptList, repertoire=NULL, winds=4:5, aa.index.id=c("MIYS990106")){
  # Length check
  peptList_th <- peptList[nchar(peptList)>=max(winds)]

  # A batch mode
  rtpcp_df <- data.frame(
    "Peptide"=peptList_th,
    dplyr::bind_cols(
      lapply(winds, function(windowSize){rTPCP.fixedWindow(peptList_th, repertoire, windowSize, aa.index.id)})
    )
  ) %>% dplyr::as_data_frame()
  return(rtpcp_df)
}
