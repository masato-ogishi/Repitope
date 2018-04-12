#' Similarity network analysis.
#'
#' @param peptideSet A set of peptide sequences.
#' @param longerPeptideSet A set of peptide sequences. Should be one amino acid longer than \code{shorterPeptideSet}.
#' @param shorterPeptideSet A set of peptide sequences. Should be one amino acid shorter than \code{longerPeptideSet}.
#' @param numSet An attribute for the vertices.
#' @param directed Should the network be converted from undirected to directed? Directions are determined by the \code{numSet} provided.
#' @param weighted Should the network be converted to weihted? Edge weights are determined by the \code{numSet} provided.
#' @param forceMutType Should mutational types be annotated? A little bit time-consuming.
#' @importFrom dplyr %>%
#' @importFrom dplyr rename
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom dplyr if_else
#' @importFrom dplyr bind_rows
#' @importFrom dplyr inner_join
#' @importFrom DescTools Sort
#' @importFrom stringr str_sub
#' @importFrom stringr str_split
#' @importFrom stringi stri_reverse
#' @importFrom S4Vectors nchar
#' @importFrom Biostrings AAStringSet
#' @importFrom Biostrings pairwiseAlignment
#' @importFrom Biostrings mismatchTable
#' @importFrom stringdist stringdist
#' @importFrom Matrix sparseMatrix
#' @importFrom Matrix t
#' @importFrom igraph V
#' @importFrom igraph E
#' @importFrom igraph set_vertex_attr
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom igraph graph_from_data_frame
#' @importFrom igraph simplify
#' @importFrom igraph as_edgelist
#' @importFrom pbapply timerProgressBar
#' @importFrom pbapply setTimerProgressBar
#' @importFrom pbapply pblapply
#' @export
#' @rdname Utility_SimilarityNetwork
#' @name Utility_SimilarityNetwork
distMat_Auto <- function(peptideSet){
  # Input check...
  if(class(peptideSet)!="AAStringSet"){
    peptideSet <- Biostrings::AAStringSet(peptideSet)
  }
  if(length(unique(S4Vectors::nchar(peptideSet)))>=2){
    message("Input sequences must be of the same length!")
    return(NULL)
  }

  sequenceSet <- as.character(peptideSet)
  DistMat.auto <- Matrix::sparseMatrix(c(), c(), x=F, dims=c(length(peptideSet), length(peptideSet)))
  pbQ <- length(peptideSet)>1000
  if(pbQ){
    pb <- pbapply::timerProgressBar(min=1, max=length(peptideSet), char="+", style=1)
    for(i in 1:length(peptideSet)){
      pbapply::setTimerProgressBar(pb, i)
      DistMat.auto[i, 1:i] <- stringdist::stringdist(sequenceSet[[i]], peptideSet[1:i], method="hamming")<=1
    }
    cat("\n")
  }else{
    for(i in 1:length(peptideSet)){
      DistMat.auto[i, 1:i] <- stringdist::stringdist(sequenceSet[[i]], peptideSet[1:i], method="hamming")<=1
    }
  }
  diag(DistMat.auto) <- F
  return(DistMat.auto)  # An m x m triangular matrix: m = length(peptideSet)
}
#' @export
#' @rdname Utility_SimilarityNetwork
#' @name Utility_SimilarityNetwork
distMat_Juxtaposed <- function(longerPeptideSet, shorterPeptideSet){
  # Input check...
  if(class(longerPeptideSet)!="AAStringSet"){
    longerPeptideSet <- Biostrings::AAStringSet(longerPeptideSet)
  }
  if(class(shorterPeptideSet)!="AAStringSet"){
    shorterPeptideSet <- Biostrings::AAStringSet(shorterPeptideSet)
  }
  longerSeq <- as.character(longerPeptideSet)
  shorterSeq <- as.character(shorterPeptideSet)
  leng_longer <- unique(nchar(longerSeq))
  leng_shorter <- unique(nchar(shorterSeq))
  if(length(leng_longer)>=2){
    message("Input longer sequences must be of the same length!")
    return(NULL)
  }
  if(length(leng_shorter)>=2){
    message("Input shorter sequences must be of the same length!")
    return(NULL)
  }
  if(leng_longer<=leng_shorter){
    message("The longer sequences must be one amino acid longer than the shorter sequences!")
    return(NULL)
  }
  if(leng_longer-leng_shorter>=2){
    message("The length difference between longer and shorter sequences must be one amino acid!")
    return(NULL)
  }

  longerSeq.degenerate <- unlist(lapply(1:leng_longer, function(i){seq <- longerSeq; stringr::str_sub(seq, i, i) <- ""; return(seq)}))
  longerSeq.degenerate <- matrix(longerSeq.degenerate, ncol=leng_longer)
  DistMat.juxt <- Matrix::sparseMatrix(c(), c(), x=F, dims=c(length(shorterPeptideSet), length(longerPeptideSet)))
  pbQ <- length(shorterPeptideSet)>1000
  if(pbQ){
    pb <- pbapply::timerProgressBar(min=1, max=length(shorterPeptideSet), char="+", style=1)
    for(i in 1:length(shorterPeptideSet)){
      pbapply::setTimerProgressBar(pb, i)
      DistMat.juxt[i,] <- rowSums(longerSeq.degenerate==shorterSeq[[i]])!=0
    }
    cat("\n")
  }else{
    for(i in 1:length(shorterPeptideSet)){
      DistMat.juxt[i,] <- rowSums(longerSeq.degenerate==shorterSeq[[i]])!=0
    }
  }
  DistMat.juxt <- Matrix::t(DistMat.juxt)
  return(DistMat.juxt) # An m x n matrix: m = length(longerPeptideSet), n = length(shorterPeptideSet)
}
#' @export
#' @rdname Utility_SimilarityNetwork
#' @name Utility_SimilarityNetwork
singleAASimilarityNetwork <- function(peptideSet, numSet=NULL, directed=T, weighted=T, forceMutType=F){
  # Internally used workflows
  net_main <- function(peptideSet){
    ## Input check...
    if(class(peptideSet)!="AAStringSet"){
      peptideSet <- Biostrings::AAStringSet(peptideSet)
    }
    if(length(unique(peptideSet))!=length(peptideSet)){
      message("Duplicates are removed from the input sequences.")
      peptideSet <- unique(peptideSet)
    }

    ## Serialized adjacency matrix calculation
    aaStringSetList <- split(peptideSet, S4Vectors::nchar(peptideSet))
    if(length(names(aaStringSetList))>=2){
      sequenceLengthPairGrid <- suppressWarnings(dplyr::bind_rows(
        data.frame(V1=names(aaStringSetList), V2=names(aaStringSetList), Type="Auto"),
        data.frame(as.data.frame(t(combn(names(aaStringSetList), 2))), Type="Juxtaposed")
      )) %>% dplyr::rename(V2="V1", V1="V2") %>% dplyr::select(V1, V2, Type)
    }else{
      sequenceLengthPairGrid <- data.frame(V1=names(aaStringSetList), V2=names(aaStringSetList), Type="Auto")
    }
    sequenceLengthPairGrid <- sequenceLengthPairGrid %>%
      dplyr::filter((as.numeric(V1)-as.numeric(V2)) %in% c(0, 1)) ## V1>=V2
    sequenceLengthPairN <- nrow(sequenceLengthPairGrid)

    ends <- as.numeric(cumsum(table(S4Vectors::nchar(peptideSet))))
    starts <- (c(0, ends)+1)[1:length(ends)]
    positionGrid <- as.data.frame(t(data.frame(starts, ends)))
    colnames(positionGrid) <- names(aaStringSetList)

    DistMat <- Matrix::sparseMatrix(c(), c(), x=F, dims=rep(length(peptideSet), 2))

    message(paste0("Number of sequence length pairs = ", sequenceLengthPairN))
    for(i in which(sequenceLengthPairGrid$"Type"=="Auto")){
      message(paste0("Pair ", i, "/", sequenceLengthPairN, " | Auto: sequence length = ", sequenceLengthPairGrid[i,1]))
      s1 <- aaStringSetList[[sequenceLengthPairGrid[i,1]]]
      p1 <- positionGrid[[sequenceLengthPairGrid[i,1]]]
      p1 <- seq(p1[1], p1[2])
      DistMat[p1, p1] <- distMat_Auto(s1)
    }
    for(i in which(sequenceLengthPairGrid$"Type"=="Juxtaposed")){
      message(paste0("Pair ", i, "/", sequenceLengthPairN, " | Juxtaposed: sequence length pair = ", sequenceLengthPairGrid[i,1], " and ", sequenceLengthPairGrid[i,2]))
      s1 <- aaStringSetList[[sequenceLengthPairGrid[i,1]]]
      p1 <- positionGrid[[sequenceLengthPairGrid[i,1]]]
      p1 <- seq(p1[1], p1[2])
      s2 <- aaStringSetList[[sequenceLengthPairGrid[i,2]]]
      p2 <- positionGrid[[sequenceLengthPairGrid[i,2]]]
      p2 <- seq(p2[1], p2[2])
      DistMat[p1, p2] <- distMat_Juxtaposed(s1, s2)
    }
    DistMat <- DistMat|Matrix::t(DistMat) ## Symmetricalization

    ## Similarity network
    simNet <- DistMat %>%
      igraph::graph_from_adjacency_matrix(mode="undirected", weighted=NULL, diag=F) %>%
      igraph::set_vertex_attr("label", value=as.character(unlist(lapply(aaStringSetList, as.character)))) %>%
      igraph::simplify()
    return(simNet)
  }
  net_pairs_DF <- function(simpleSimNet, forceMutType=F){
    ## Get peptide pairs
    df <- as.data.frame(igraph::as_edgelist(simpleSimNet, names=T))
    colnames(df) <- c("Node1","Node2")
    df[["AASeq1"]] <- igraph::V(simpleSimNet)$label[df$"Node1"]
    df[["AASeq2"]] <- igraph::V(simpleSimNet)$label[df$"Node2"]
    df <- dplyr::select(df, -Node1, -Node2)

    ## Annotate mutational types and patterns
    mutType <- function(pept1, pept2){
      if(nchar(pept1)==nchar(pept2)){
        pwal <- Biostrings::pairwiseAlignment(pattern=pept1, subject=pept2, type="global")
        sbst <- Biostrings::mismatchTable(pwal)
        sbst <- paste0(c(as.character(sbst$PatternSubstring), sbst$PatternStart, as.character(sbst$SubjectSubstring)), collapse="_")
        return(sbst)
      }
      if(nchar(pept1)<nchar(pept2)){
        shorterseq <- pept1
        longerSeq <- pept2
        leng_longer <- nchar(longerSeq)
        longerSeq.degenerate <- unlist(lapply(1:leng_longer, function(i){seq <- longerSeq; stringr::str_sub(seq, i, i) <- ""; return(seq)}))
        longerSeq.degenerate <- matrix(longerSeq.degenerate, ncol=leng_longer)
        del.pos <- which(longerSeq.degenerate==shorterseq)
        sbst <- sapply(del.pos, function(p){paste0(c("-", p, substr(longerSeq, p, p)), collapse="_")})
        sbst <- paste0(sbst, collapse="|")
        return(sbst)
      }
      if(nchar(pept1)>nchar(pept2)){
        shorterseq <- pept2
        longerSeq <- pept1
        leng_longer <- nchar(longerSeq)
        longerSeq.degenerate <- unlist(lapply(1:leng_longer, function(i){seq <- longerSeq; stringr::str_sub(seq, i, i) <- ""; return(seq)}))
        longerSeq.degenerate <- matrix(longerSeq.degenerate, ncol=leng_longer)
        del.pos <- which(longerSeq.degenerate==shorterseq)
        sbst <- sapply(del.pos, function(p){paste0(c(substr(longerSeq, p, p), p, "-"), collapse="_")})
        sbst <- paste0(sbst, collapse="|")
        return(sbst)
      }
    }
    mutPattern <- function(pept1, pept2){
      if(nchar(pept1)==nchar(pept2)) return("Substitution")
      if(nchar(pept1)<nchar(pept2)) return("Insertion")
      if(nchar(pept1)>nchar(pept2)) return("Deletion")
    }
    message("Annotating mutational types and patterns...")
    if(forceMutType==T|nrow(df)<=50000){
      df[["MutType"]] <- unlist(pbapply::pblapply(1:nrow(df), function(i){mutType(df$"AASeq1"[[i]], df$"AASeq2"[[i]])})) ## a bit slow...
    }else{
      df[["MutType"]] <- NA
    }
    df[["MutPattern"]] <- unlist(pbapply::pblapply(1:nrow(df), function(i){mutType(df$"AASeq1"[[i]], df$"AASeq2"[[i]])}))

    ## Output
    return(df)
  }

  # A basic, non-directional, non-weighted network
  simNet <- net_main(peptideSet)
  simNet_PairsDF <- net_pairs_DF(simNet, forceMutType=forceMutType)
  if(is.null(numSet)) return(list("SimilarityNetwork"=simNet, "PairsDF"=simNet_PairsDF))

  # A directional, weighted network
  df_num <- data.frame("Peptide"=as.character(peptideSet), "Score"=numSet, stringsAsFactors=F)
  simNet_PairsDF_DW <- simNet_PairsDF %>%
    dplyr::inner_join(df_num, by=c("AASeq1"="Peptide")) %>%
    dplyr::inner_join(df_num, by=c("AASeq2"="Peptide")) %>%
    dplyr::transmute(AASeq1, AASeq2, Score1=Score.x, Score2=Score.y, MutType, MutPattern)
  mutType_Rev <- function(mutType){
    paste0(lapply(strsplit(mutType, "|", fixed=T)[[1]], function(mut){stringi::stri_reverse(mut)}), collapse="|")
  }
  simNet_PairsDF_DW <- dplyr::bind_rows(
    simNet_PairsDF_DW,
    simNet_PairsDF_DW %>%
      dplyr::rename(AASeq1=AASeq2, AASeq2=AASeq1, Score1=Score2, Score2=Score1) %>%
      dplyr::mutate(MutType=sapply(MutType, mutType_Rev))
  ) %>%
    dplyr::filter(Score1>=Score2) %>%
    dplyr::transmute(AASeq1, AASeq2, Score1, Score2, DeltaScore=Score1-Score2, MutType, MutPattern)
  simNet_DW <- igraph::graph_from_data_frame(simNet_PairsDF_DW, directed=directed)
  simNet_DW <- igraph::set_vertex_attr(simNet_DW, name="label", value=igraph::V(simNet_DW)$name)
  if(weighted==T) igraph::E(simNet_DW)$weight <- igraph::E(simNet_DW)$DeltaScore
  return(list("SimilarityNetwork"=simNet, "PairsDF"=simNet_PairsDF,
              "SimilarityNetwork_DW"=simNet_DW, "PairsDF_DW"=simNet_PairsDF_DW))
}
