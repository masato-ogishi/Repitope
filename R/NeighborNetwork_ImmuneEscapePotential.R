#' Immune escape potentials derived from neighbor network clustering analysis.
#'
#' @param peptide The target peptide sequence.
#' @param neighborNetResult The result returned from \code{neighborNetwork}.
#' @param metadataDF A dataframe which has "Peptide", "Immunogenicity", and "ImmunogenicityScore" columns.
#' @param seed A random seed.
#' @importFrom dplyr %>%
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom dplyr tibble
#' @importFrom dplyr left_join
#' @importFrom dplyr arrange
#' @importFrom dplyr mutate
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr if_else
#' @importFrom magrittr set_colnames
#' @importFrom data.table as.data.table
#' @importFrom igraph induced_subgraph
#' @importFrom igraph V
#' @importFrom igraph E
#' @importFrom igraph layout_nicely
#' @importFrom igraph cluster_walktrap
#' @importFrom scales alpha
#' @importFrom scales rescale
#' @importFrom ggseqlogo geom_logo
#' @importFrom ggseqlogo theme_logo
#' @importFrom ggpubr rremove
#' @import ggplot2
#' @export
#' @rdname NeighborNetwork_ImmuneEscapePotential
#' @name NeighborNetwork_ImmuneEscapePotential
neighborNetwork_ClusterPlot <- function(peptide, neighborNetResult, metadataDF, seed=12345){
  ## Subgraph
  connectedPeptides <- c(
    peptide,
    dplyr::filter(neighborNetResult$"PairDF_DW", AASeq1==peptide)$"AASeq2",
    dplyr::filter(neighborNetResult$"PairDF_DW", AASeq2==peptide)$"AASeq1"
  )
  connectedSubgraph <- igraph::induced_subgraph(
    neighborNetResult$"NeighborNetwork_DW",
    which(igraph::V(neighborNetResult$"NeighborNetwork_DW")$label %in% connectedPeptides)
  )

  ## Metadata
  peptideLabels <- peptide
  peptideSet <- igraph::V(connectedSubgraph)$"name"
  df_meta <- dplyr::filter(metadataDF, Peptide %in% peptideSet) %>%
    dplyr::select(Peptide, Immunogenicity, ImmunogenicityScore)
  igraph::V(connectedSubgraph)$label[!peptideSet %in% peptideLabels] <- ""
  igraph::V(connectedSubgraph)$Immunogenicity <- df_meta$Immunogenicity
  igraph::V(connectedSubgraph)$ImmunogenicityScore <- df_meta$ImmunogenicityScore
  igraph::E(connectedSubgraph)$weight[igraph::E(connectedSubgraph)$weight==0] <- 0.001

  ## Coordinates
  set.seed(seed)
  l <- igraph::layout_nicely(connectedSubgraph)
  df_meta <- dplyr::left_join(
    magrittr::set_colnames(cbind(dplyr::tibble("Peptide"=peptideSet), as.data.frame(l)), c("Peptide","x","y")),
    df_meta,
    by="Peptide"
  )

  ## Clusters
  set.seed(seed)
  df_meta <- dplyr::left_join(
    dplyr::tibble("Peptide"=peptideSet,
                  "ClusterID"=paste0("Cluster", igraph::cluster_walktrap(connectedSubgraph)$membership)),
    df_meta,
    by="Peptide"
  )
  clusteredPeptides <- df_meta %>%
    dplyr::arrange(ClusterID) %>%
    dplyr::mutate(ClusterID=as.character(ClusterID)) %>%
    data.table::as.data.table() %>%
    split(by="ClusterID") %>%
    lapply(function(d){d[["Peptide"]]})
  seqLogoPlot <- ggplot() +
    ggseqlogo::geom_logo(clusteredPeptides, seq_type="aa") +
    facet_wrap(~seq_group, nrow=length(clusteredPeptides)) +
    ggseqlogo::theme_logo() +
    ggpubr::rremove("y.text") +
    ggpubr::rremove("y.title") +
    theme(legend.position="right", legend.direction="vertical")

  ## Graph plot
  clusterGraphPlot <- function(g, meta, layout, seed=12345){
    ### Colors
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

    ### Labels
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
      vertex.shape=dplyr::if_else(meta$Immunogenicity=="Positive", "circle", "square"),
      vertex.label=igraph::V(g)$"label",
      vertex.label.color="black",
      vertex.label.cex=1.0,
      vertex.label.dist=0.5,
      vertex.label.family="sans",
      vertex.color=vertexColors,
      vertex.frame.color="black",
      edge.width=0.1,
      edge.arrow.size=0.1,
      edge.arrow.width=0.1,
      edge.color=scales::alpha("gray50", 0.75)
    )
    text(
      clusterCentroids$x,
      clusterCentroids$y,
      clusterCentroids$ClusterID,
      cex=1.25,
      family="sans"
    )
    return(recordPlot())
  }
  neighborPlot <- clusterGraphPlot(
    g=connectedSubgraph,
    meta=df_meta,
    layout=l,
    seed=seed
  )
  return(list(
    "NeighborNetwork"=connectedSubgraph,
    "SummaryDF"=df_meta,
    "NeighborPlot"=neighborPlot,
    "SeqLogoPlot"=seqLogoPlot
  ))
}
