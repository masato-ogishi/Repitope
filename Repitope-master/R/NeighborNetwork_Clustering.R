#' Neighbor network clustering analysis.
#'
#' \code{neighborNetwork_ConnectedSubGraph} and \code{neighborNetwork_ConnectedSubGraphDF} extract the minimum connected subgraph focusing on a target peptide.\cr
#' \code{neighborNetwork_Cluster} and \code{neighborNetwork_Cluster_Batch} conduct a network clustering analysis using a walktrap algorithm.\cr
#' \code{neighborNetwork_Cluster_FeatureDF} computes per-cluster and per-peptide features, including "escape potential".
#'
#' @param neighborNetResult The result returned from \code{neighborNetwork}.
#' @param neighborNetClusterResult The result returned from \code{neighborNetwork_Cluster_Batch}.
#' @param peptide A target peptide sequence.
#' @param graph A directed and weighted neighbor network of the target peptide.
#' @param metadataDF A dataframe which has "Peptide" and "ImmunogenicityScore" columns. Optionally, an "Immunogenicity" column would be integrated if present.
#' @param seed A random seed.
#' @param plot Logical. Whether the network cluster plot shuould be generated.
#' @param coreN The number of threads for parallelization.
#' @export
#' @rdname NeighborNetwork_Clustering
#' @name NeighborNetwork_Clustering
neighborNetwork_ConnectedSubGraphDF <- function(neighborNetResult, peptide){
  dplyr::bind_rows(
    neighborNetResult$"PairDF_DW" %>%
      dplyr::filter(AASeq1==peptide),
    neighborNetResult$"PairDF_DW" %>%
      dplyr::filter(AASeq2==peptide) %>%
      dplyr::mutate(ScoreRatio=1/ScoreRatio)
  ) %>%
    tidyr::separate(col="MutType", into=c("AA1","Pos","AA2"))
}

#' @export
#' @rdname NeighborNetwork_Clustering
#' @name NeighborNetwork_Clustering
neighborNetwork_ConnectedSubGraph <- function(neighborNetResult, peptide){
  connectedPeptides <- c(
    peptide,
    dplyr::filter(neighborNetResult$"PairDF_DW", AASeq1==peptide)$"AASeq2",
    dplyr::filter(neighborNetResult$"PairDF_DW", AASeq2==peptide)$"AASeq1"
  )
  connectedSubgraph <- igraph::induced_subgraph(
    neighborNetResult$"NeighborNetwork_DW",
    which(igraph::V(neighborNetResult$"NeighborNetwork_DW")$label %in% connectedPeptides)
  )
  return(connectedSubgraph)
}

#' @export
#' @rdname NeighborNetwork_Clustering
#' @name NeighborNetwork_Clustering
neighborNetwork_Cluster <- function(peptide, graph, metadataDF, seed=12345, plot=T){
  ## Peptide labels
  peptideLabels <- peptide
  peptideSet <- igraph::V(graph)$"name"
  igraph::V(graph)$label[!peptideSet %in% peptideLabels] <- ""

  ## Metadata
  metadataDF <- dplyr::filter(metadataDF, Peptide %in% peptideSet)
  df_meta <- dplyr::select(metadataDF, Peptide, ImmunogenicityScore)
  if("Immunogenicity" %in% colnames(metadataDF)){
    df_meta$Immunogenicity <- metadataDF$Immunogenicity
    igraph::V(graph)$Immunogenicity <- df_meta$Immunogenicity
  }
  igraph::V(graph)$ImmunogenicityScore <- df_meta$ImmunogenicityScore

  ## Clusters
  set.seed(seed)
  df_meta <- dplyr::left_join(
    dplyr::tibble("Peptide"=peptideSet,
                  "Target"=dplyr::if_else(peptideSet==peptide, "Target", "Neighbor"),
                  "ClusterID"=paste0("Cluster", igraph::cluster_walktrap(graph)$membership)),
    df_meta,
    by="Peptide"
  )

  ## No plot ver.
  if(plot!=T) return(df_meta)

  ## Coordinates
  set.seed(seed)
  l <- igraph::layout_nicely(graph)
  df_meta <- dplyr::left_join(
    magrittr::set_colnames(cbind(dplyr::tibble("Peptide"=peptideSet), as.data.frame(l)), c("Peptide","x","y")),
    df_meta,
    by="Peptide"
  )

  ## Consensus per cluster
  clusteredPeptides <- df_meta %>%
    dplyr::arrange(ClusterID) %>%
    dplyr::mutate(ClusterID=as.character(ClusterID)) %>%
    data.table::as.data.table() %>%
    split(by="ClusterID") %>%
    lapply(function(d){d[["Peptide"]]})
  consensusSequence <- function(sequenceSet){
    if(length(sequenceSet)==1) return(sequenceSet)
    sink(tempfile())
    s <- msa::msaConsensusSequence(msa::msaClustalW(sequenceSet, type="protein"), type="Biostrings", ambiguityMap="X", threshold=0.5)
    sink()
    return(s)
  }
  clusterConsensusSeqs <- sapply(clusteredPeptides, consensusSequence)
  clusterConsensusSeqs <- stringr::str_replace_all(clusterConsensusSeqs, stringr::fixed("-"), "X")

  ## Graph plot
  clusterGraphPlot <- function(g, meta, layout, seed=12345){
    ### Vertex annotations
    colPal <- function(x){
      pal <- colorRamp(append(ggsci::pal_d3()(2), "white", after=1), space="rgb")
      cols <- pal(x)
      apply(cols, 1, function(x){rgb(x[1], x[2], x[3], maxColorValue=255)})
    }
    clusterColors <- meta %>%
      dplyr::group_by(ClusterID) %>%
      dplyr::summarise(ImmunogenicityColor=colPal(mean(ImmunogenicityScore)))
    clusterColors <- scales::alpha(clusterColors$"ImmunogenicityColor", alpha=0.75)
    vertexColors <- colPal(meta$"ImmunogenicityScore")
    if(is.null(meta$"Immunogenicity")){
      vertexShapes <- "circle"
    }else{
      vertexShapes <- dplyr::if_else(meta$"Immunogenicity"=="Positive", "circle", "square")
    }
    vertexLabels <- igraph::V(g)$"label"

    ### Cluster labels
    clusterCentroids <- meta %>%
      dplyr::group_by(ClusterID) %>%
      dplyr::summarise(x=mean(x), y=mean(y))
    clusterCentroids$x <- scales::rescale(clusterCentroids$x, to=c(-1, 1))
    clusterCentroids$y <- scales::rescale(clusterCentroids$y, to=c(-1, 1))

    ### Graph
    try(dev.off(), silent=T)
    set.seed(seed)
    plot(
      igraph::cluster_walktrap(g), g,
      layout=layout,
      col=vertexColors,
      mark.border="black",
      mark.col=clusterColors,
      vertex.size=3,
      vertex.shape=vertexShapes,
      vertex.label=vertexLabels,
      vertex.label.color="black",
      vertex.label.cex=1.25,
      vertex.label.dist=0.5,
      vertex.label.family="sans",
      vertex.color=vertexColors,
      vertex.frame.color="black",
      edge.width=0.5,
      edge.arrow.size=0.25,
      edge.arrow.width=0.1,
      edge.color=scales::alpha("gray50", 0.75)
    )
    for(i in 1:length(clusterCentroids$ClusterID)){
      text(
        rep(clusterCentroids$x[[i]], 2),
        clusterCentroids$y[[i]]+c(0, -0.1),
        c(clusterCentroids$ClusterID[[i]], clusterConsensusSeqs[[i]]),
        pos=3, cex=1.25, family="sans"
      )
    }
    return(recordPlot())
  }
  neighborPlot <- clusterGraphPlot(
    g=graph,
    meta=df_meta,
    layout=l,
    seed=seed
  )
  return(list(
    "SummaryDF"=df_meta,
    "NeighborPlot"=neighborPlot
  ))
}

#' @export
#' @rdname NeighborNetwork_Clustering
#' @name NeighborNetwork_Clustering
neighborNetwork_Cluster_Batch <- function(neighborNetResult, metadataDF, seed=12345, coreN=parallel::detectCores(logical=F)){
  ## Extract peptide sequences
  peptideSet <- igraph::V(neighborNetResult$"NeighborNetwork_DW")$"name"

  ## Cluster analysis
  message("Clustering neighbor network...")
  cl <- parallel::makeCluster(coreN, type="PSOCK")
  snow::clusterSetupRNGstream(cl, seed=rep(seed, 6))
  doSNOW::registerDoSNOW(cl)
  res <- foreach::foreach(pept=peptideSet, .inorder=F)%dopar%{
    graph <- Repitope::neighborNetwork_ConnectedSubGraph(neighborNetResult, pept)
    df_meta <- Repitope::neighborNetwork_Cluster(pept, graph, metadataDF, seed, plot=F)
    pos <- grep("Target", df_meta$"Target")
    clust <- df_meta$"ClusterID"[[pos]]
    score <- df_meta$"ImmunogenicityScore"[[pos]]
    list("Peptide"=pept, "Score"=score, "ClusterID"=clust, "SummaryDF"=df_meta)
  }
  parallel::stopCluster(cl)
  gc();gc()
  return(res)
}

#' @export
#' @rdname NeighborNetwork_Clustering
#' @name NeighborNetwork_Clustering
neighborNetwork_Cluster_FeatureDF <- function(neighborNetClusterResult, coreN=parallel::detectCores(logical=F)){
  message("Computing cluster-based metrics...")
  cl <- parallel::makeCluster(coreN, type="PSOCK")
  doSNOW::registerDoSNOW(cl)
  dt_feat <- foreach::foreach(i=1:length(neighborNetClusterResult), .inorder=F, .packages=c("dplyr","data.table"))%dopar%{
    summaryDF <- neighborNetClusterResult[[i]][["SummaryDF"]]
    summaryDF <- dplyr::arrange(summaryDF, dplyr::desc(Target))
    dt <- data.table::as.data.table(summaryDF[1,])
    aveDF <- summaryDF %>%
      dplyr::group_by(ClusterID) %>%
      dplyr::summarise(ImmunogenicityScore=mean(ImmunogenicityScore))
    dt[,ImmunogenicityScore_Cluster_Average:=dplyr::filter(aveDF, ClusterID==dt$"ClusterID"[1])$"ImmunogenicityScore"]
    aveDF <- dplyr::filter(aveDF, ClusterID!=dt$"ClusterID"[1])
    dt[,ImmunogenicityScore_Cluster_Diff_Max:=max(aveDF$"ImmunogenicityScore") - ImmunogenicityScore_Cluster_Average]
    dt[,ImmunogenicityScore_Cluster_Diff_Min:=ImmunogenicityScore_Cluster_Average - min(aveDF$"ImmunogenicityScore")]  ## Considered to be an "escape potential"
    dt[,ImmunogenicityScore_Diff_Max:=max(summaryDF$"ImmunogenicityScore") - ImmunogenicityScore]
    dt[,ImmunogenicityScore_Diff_Min:=ImmunogenicityScore - min(summaryDF$"ImmunogenicityScore")]
    dt
  } %>%
    data.table::rbindlist()
  parallel::stopCluster(cl)
  gc();gc()
  return(dt_feat)
}
