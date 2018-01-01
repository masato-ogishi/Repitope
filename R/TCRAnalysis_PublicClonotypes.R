#' Public XCR clonotypes in the similarity network
#'
#' @param aaStringSet A set of XCR clonotype sequences (converted as an AAStringSet object).
#' @param longerAAStringSet A set of XCR clonotype sequences. Should be one amino acid longer than \code{shorterAAStringSet}.
#' @param shorterAAStringSet A set of XCR clonotype sequences. Should be one amino acid shorter than \code{longerAAStringSet}.
#' @param numSet An attribute for the verteces.
#' @param directed Should the network be converted from undirected to directed? Directions are determined by the \code{numSet} provided.
#' @param weighted Should the network be converted to weihted? Edge weights are determined by the \code{numSet} provided.
#' @param simNet A similarity network object.
#' @param sizeThresholdSet A set of integers specifying the maximum numbers of clonotypes retained as public clonotypes. Based on the degree distribution, appropriate degree threshold will be internally calculated.
#' @param outputFileHeader A path for saving the outputs from the public clonotype analysis.
#' @importFrom dplyr %>%
#' @importFrom dplyr rename
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom dplyr if_else
#' @importFrom dplyr bind_rows
#' @importFrom dplyr inner_join
#' @importFrom DescTools Sort
#' @importFrom stringr str_sub
#' @importFrom S4Vectors nchar
#' @importFrom Biostrings AAStringSet
#' @importFrom Biostrings pairwiseAlignment
#' @importFrom Biostrings mismatchTable
#' @importFrom stringdist stringdist
#' @importFrom Matrix sparseMatrix
#' @importFrom Matrix t
#' @importFrom igraph set_vertex_attr
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom igraph graph_from_data_frame
#' @importFrom igraph simplify
#' @importFrom igraph as_edgelist
#' @importFrom igraph decompose
#' @importFrom igraph V
#' @importFrom igraph E
#' @importFrom igraph write_graph
#' @importFrom igraph degree
#' @importFrom igraph subgraph
#' @importFrom pbapply timerProgressBar
#' @importFrom pbapply setTimerProgressBar
#' @importFrom pbapply pblapply
#' @export
#' @rdname TCRAnalysis_PublicClonotypes
#' @name TCRAnalysis_PublicClonotypes
distMat_Auto <- function(aaStringSet){
  # Input check...
  if(class(aaStringSet)!="AAStringSet"){
    message("Input sequences are converted to AAStringSet object.")
    aaStringSet <- Biostrings::AAStringSet(aaStringSet)
  }
  if(length(unique(S4Vectors::nchar(aaStringSet)))>=2){
    message("Input sequences must be of the same length!")
    return(NULL)
  }

  sequenceSet <- as.character(aaStringSet)
  DistMat.auto <- Matrix::sparseMatrix(c(), c(), x=F, dims=c(length(aaStringSet), length(aaStringSet)))
  pbQ <- length(aaStringSet)>1000
  if(pbQ){
    pb <- pbapply::timerProgressBar(min=1, max=length(aaStringSet), char="+", style=1)
    for(i in 1:length(aaStringSet)){
      pbapply::setTimerProgressBar(pb, i)
      DistMat.auto[i, 1:i] <- stringdist::stringdist(sequenceSet[[i]], aaStringSet[1:i], method="hamming")<=1
    }
    cat("\n")
  }else{
    for(i in 1:length(aaStringSet)){
      DistMat.auto[i, 1:i] <- stringdist::stringdist(sequenceSet[[i]], aaStringSet[1:i], method="hamming")<=1
    }
  }
  diag(DistMat.auto) <- F
  return(DistMat.auto)  # An m x m triangular matrix: m = length(aaStringSet)
}
#' @export
#' @rdname TCRAnalysis_PublicClonotypes
#' @name TCRAnalysis_PublicClonotypes
distMat_Juxtaposed <- function(longerAAStringSet, shorterAAStringSet){
  # Input check...
  if(class(longerAAStringSet)!="AAStringSet"){
    message("Input sequences are converted to AAStringSet object.")
    longerAAStringSet <- Biostrings::AAStringSet(longerAAStringSet)
  }
  if(class(shorterAAStringSet)!="AAStringSet"){
    message("Input sequences are converted to AAStringSet object.")
    shorterAAStringSet <- Biostrings::AAStringSet(shorterAAStringSet)
  }
  longerSeq <- as.character(longerAAStringSet)
  shorterSeq <- as.character(shorterAAStringSet)
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
  DistMat.juxt <- Matrix::sparseMatrix(c(), c(), x=F, dims=c(length(shorterAAStringSet), length(longerAAStringSet)))
  pbQ <- length(shorterAAStringSet)>1000
  if(pbQ){
    pb <- pbapply::timerProgressBar(min=1, max=length(shorterAAStringSet), char="+", style=1)
    for(i in 1:length(shorterAAStringSet)){
      pbapply::setTimerProgressBar(pb, i)
      DistMat.juxt[i,] <- rowSums(longerSeq.degenerate==shorterSeq[[i]])!=0
    }
    cat("\n")
  }else{
    for(i in 1:length(shorterAAStringSet)){
      DistMat.juxt[i,] <- rowSums(longerSeq.degenerate==shorterSeq[[i]])!=0
    }
  }
  DistMat.juxt <- Matrix::t(DistMat.juxt)
  return(DistMat.juxt) # An m x n matrix: m = length(longerAAStringSet), n = length(shorterAAStringSet)
}
#' @export
#' @rdname TCRAnalysis_PublicClonotypes
#' @name TCRAnalysis_PublicClonotypes
singleAASimilarityNetwork <- function(aaStringSet, numSet=NULL, directed=T, weighted=T){
  # Internally used workflows
  net_main <- function(aaStringSet){
    ## Input check...
    if(class(aaStringSet)!="AAStringSet"){
      message("Input sequences are converted to AAStringSet object.")
      aaStringSet <- Biostrings::AAStringSet(aaStringSet)
    }
    if(length(unique(aaStringSet))!=length(aaStringSet)){
      message("Duplicates are removed from the input sequences.")
      aaStringSet <- unique(aaStringSet)
    }

    ## Serialized adjacency matrix calculation
    aaStringSetList <- split(aaStringSet, S4Vectors::nchar(aaStringSet))
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

    ends <- as.numeric(cumsum(table(S4Vectors::nchar(aaStringSet))))
    starts <- (c(0, ends)+1)[1:length(ends)]
    positionGrid <- as.data.frame(t(data.frame(starts, ends)))
    colnames(positionGrid) <- names(aaStringSetList)

    DistMat <- Matrix::sparseMatrix(c(), c(), x=F, dims=rep(length(aaStringSet), 2))

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
  net_pairs_DF <- function(simpleSimNet){
    ## Get peptide pairs
    df <- as.data.frame(igraph::as_edgelist(simpleSimNet, names=T))
    colnames(df) <- c("Node1","Node2")
    df[["AASeq1"]] <- igraph::V(simpleSimNet)$label[df$"Node1"]
    df[["AASeq2"]] <- igraph::V(simpleSimNet)$label[df$"Node2"]

    ## Annotate mutational types and patterns
    mutType <- function(pept1, pept2){
      if(nchar(pept1)==nchar(pept2)){
        pwal <- Biostrings::pairwiseAlignment(pattern=pept1, subject=pept2, type="global")
        sbst <- Biostrings::mismatchTable(pwal)
        sbst <- paste0(sbst$PatternSubstring, sbst$PatternStart, sbst$SubjectSubstring, collapse="")
        return(sbst)
      }
      if(nchar(pept1)<nchar(pept2)){
        shorterseq <- pept1
        longerSeq <- pept2
        leng_longer <- nchar(longerSeq)
        longerSeq.degenerate <- unlist(lapply(1:leng_longer, function(i){seq <- longerSeq; stringr::str_sub(seq, i, i) <- ""; return(seq)}))
        longerSeq.degenerate <- matrix(longerSeq.degenerate, ncol=leng_longer)
        del.pos <- which(longerSeq.degenerate==shorterseq)
        sbst <- sapply(del.pos, function(p){paste0("-", p, substr(pept2, p, p), collapse="")})
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
        sbst <- sapply(del.pos, function(p){paste0(substr(pept2, p, p), p, "-", collapse="")})
        sbst <- paste0(sbst, collapse="|")
        return(sbst)
      }
    }
    message("Annotating mutational types and patterns...")
    df[["MutType"]] <- unlist(pbapply::pblapply(1:nrow(df), function(i){mutType(df$"AASeq1"[i], df$"AASeq2"[i])}))
    df[["MutPattern"]] <- dplyr::if_else(nchar(df$"AASeq1")==nchar(df$"AASeq2"), "Substitution", "Indel")

    ## Output
    return(df)
  }

  # A basic, non-directional, non-weighted network
  simNet <- net_main(aaStringSet)
  simNet_PairsDF <- net_pairs_DF(simNet)
  if(is.null(numSet)){
    return(list("SimilarityNetwork"=simNet, "PairsDF"=simNet_PairsDF))
  }

  # A directional, weighted network
  df_num <- data.frame("Peptide"=as.character(aaStringSet), "Score"=numSet, stringsAsFactors=F)
  simNet_PairsDF_Directional <- simNet_PairsDF %>%
    dplyr::inner_join(df_num, by=c("AASeq1"="Peptide")) %>%
    dplyr::inner_join(df_num, by=c("AASeq2"="Peptide")) %>%
    dplyr::transmute(AASeq1, AASeq2, Score1=Score.x, Score2=Score.y, MutType, MutPattern)
  simNet_PairsDF_Directional <- dplyr::bind_rows(
    simNet_PairsDF_Directional,
    dplyr::transmute(simNet_PairsDF_Directional, AASeq1=AASeq2, AASeq2=AASeq1, Score1=Score2, Score2=Score1, MutType, MutPattern)
  )
  simNet_Directional <- simNet_PairsDF_Directional %>%
    dplyr::transmute(AASeq1, AASeq2, DeltaScore=Score1-Score2, MutType, MutPattern) %>%
    igraph::graph_from_data_frame(directed=T)
  simNet_Directional <- simNet_Directional %>%
    igraph::set_vertex_attr("label", value=igraph::V(simNet_Directional)$name)
  if(weighted==T) igraph::E(simNet_Directional)$weight <- igraph::E(simNet_Directional)$DeltaScore
  return(list("SimilarityNetwork"=simNet, "PairsDF"=simNet_PairsDF,
              "SimilarityNetwork_Directional"=simNet_Directional, "PairsDF_Directional"=simNet_PairsDF_Directional))
}

#' @export
#' @rdname TCRAnalysis_PublicClonotypes
#' @name TCRAnalysis_PublicClonotypes
publicClonotypeAnalysis <- function(simNet, sizeThresholdSet=10000, outputFileHeader="./Graph/Graph_"){
  # The largest connected component
  gList <- igraph::decompose(simNet)
  gSizeList <- sapply(gList, function(g){length(igraph::V(g))})
  message("Node size of the largest connected component:", max(gSizeList))
  igraph::write_graph(gList[[which(gSizeList==max(gSizeList))]],
                      file=paste0(outputFileHeader, "LargestConnected.graphml"),
                      format="graphml")

  # Split by degree
  dg <- igraph::degree(simNet)
  dg.df <- DescTools::Sort(data.frame("AASeq"=igraph::V(simNet)$"label", "Degree"=dg), ord="Degree", decreasing=T)
  message("Node with the maximum degree: ", igraph::V(simNet)$"label"[dg==max(dg)], ", degree = ", max(dg))
  subgraphByDegree <- function(th){
    minDegree <- dg.df[th,][["Degree"]]
    gSub <- suppressWarnings(igraph::subgraph(simNet, which(dg>=minDegree)))
    message(
      "Upper limit of the number of nodes = ", th, "\n",
      "The minimum degree threshold = ", minDegree, "\n",
      "The number of nodes selected = ", length(igraph::V(gSub))
    )
    sink(paste0(outputFileHeader, "PubCloneAnalysis_Log", th, ".txt"))
    cat("Public clonotype analysis.\n", sep="")
    cat(
      "Upper limit of the number of nodes = ", th, "\n",
      "The minimum degree threshold = ", minDegree, "\n",
      "The number of nodes selected = ", length(igraph::V(gSub)), sep=""
    )
    sink()
    igraph::write_graph(gSub,
                        file=paste0(outputFileHeader, "Degree", minDegree, ".graphml"),
                        format="graphml")
    clones <- igraph::V(gSub)$"label"
    write.table(clones,
                file=paste0(outputFileHeader, "PubClones", th, ".txt"),
                quote=F, sep="\t", row.names=F, col.names=F)
    return(clones)
  }
  cloneList <- lapply(sizeThresholdSet, subgraphByDegree)
  return(cloneList)
}
