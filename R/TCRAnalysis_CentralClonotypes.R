#' "Central" clonotypes in the similarity repertoire network
#'
#' @param simNet A similarity network object.
#' @param sizeSet A set of integers specifying the minimum numbers of clonotypes retained as central clonotypes. Based on the degree distribution, appropriate degree threshold will be internally calculated.
#' @param outputFileHeader A path for saving the outputs from the central clonotype analysis.
#' @importFrom dplyr %>%
#' @importFrom igraph decompose
#' @importFrom igraph V
#' @importFrom igraph E
#' @importFrom igraph write_graph
#' @importFrom igraph degree
#' @importFrom igraph subgraph
#' @importFrom DescTools Sort
#' @export
#' @rdname TCRAnalysis_CentralClonotypes
#' @name TCRAnalysis_CentralClonotypes
centralClonotypeAnalysis <- function(simNet, sizeSet=10000, outputFileHeader="./Graph/Graph_"){
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
    cat("Central clonotype analysis.\n", sep="")
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
  cloneList <- lapply(sizeSet, subgraphByDegree)
  return(cloneList)
}
