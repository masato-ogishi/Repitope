#' Public TCR clonotypes contributing as central nodes in a similarity network
#' 
#' \code{similarityNetworkCentralNodes} is a working function to extract central nodes from TCR CDR3 similarity network, where any pair of clonotypes with only one amino acid difference (levenstein distance of one) is considered linked.\cr
#' \code{similarityNetworkCentralNodes.Large} randomely splits the input TCR set, performs central node analysis, extracts important TCRs in each similarity networks, and finally performs central node analysis again with the set of important TCRs.\cr
#' \code{publicTCRAnalysis} is a wrapper function.
#' 
#' @param TCRSet A set of TCR sequences.
#' @param largeTCRSet A large set of TCR sequences.
#' @param combinedTCRDF A dataframe of TCR clonotypes. Must contain columns labeled as "aaSeqCDR3" and "N". Column "N" represents the number of datasets in which the CDR3 sequence were observed.
#' @param n.data An integer representing the threshold of publicity. TCRs observed in at least \code{n.data} repertoire datasets are considered public.
#' @param n.max.clonotype An integer representing the maximum number of clonotypes for construction of similarity network.
#' @param n.max.iteration An integer representing the maximum number of iterations.
#' @param seed A random seed.
#' @param outDir A directory in which the public clonotype analysis results are saved.
#' @param header A filename header for public clonotype analysis results.
#' @importFrom stringdist stringdistmatrix
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom igraph set_vertex_attr
#' @importFrom igraph V
#' @importFrom igraph cluster_fast_greedy
#' @importFrom igraph eigen_centrality
#' @importFrom dplyr %>%
#' @importFrom dplyr group_by
#' @importFrom dplyr top_n
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom seqinr write.fasta
#' @importFrom parallel splitIndices
#' @importFrom tibble data_frame
#' @export
#' @rdname publicTCRAnalysis
#' @name publicTCRAnalysis
similarityNetworkCentralNodes <- function(TCRSet){
  TCRSet <- sequenceFilter(TCRSet)
  tcr_distMat <- stringdist::stringdistmatrix(TCRSet, method="lv") %>% as.matrix()
  tcr_distMat[(tcr_distMat!=1)] <- 0
  tcr_g <- igraph::graph_from_adjacency_matrix(tcr_distMat, mode="undirected", weighted=NULL, diag=F) %>%
    igraph::set_vertex_attr("label", value=TCRSet)
  tcr_clust_df <- tibble::data_frame("TCR"=igraph::V(tcr_g)$label,
                                     "Cluster"=igraph::cluster_fast_greedy(tcr_g)$membership,
                                     "EigenCentrality"=igraph::eigen_centrality(tcr_g)$vector) %>%
    dplyr::group_by(Cluster) %>% dplyr::top_n(1, EigenCentrality) ## Multiple TCRs having the same centrality will be retained.
  rm(tcr_distMat);gc();gc()
  return(list("Graph"=tcr_g, "CentralityDF"=tcr_clust_df))
}
#' @export
#' @rdname publicTCRAnalysis
#' @name publicTCRAnalysis
similarityNetworkCentralNodes.Large <- function(largeTCRSet, n.max.clonotype=10000, n.max.iteration=20, seed=12345){
  # In case where the input TCR set is not so large.
  if(length(largeTCRSet)<=n.max.clonotype){
    message("The size of the TCR set is not as large as the specified upper limit. No random splitting was performed.")
    return(similarityNetworkCentralNodes(largeTCRSet))
  }

  # A working function of central node extraction
  centralNodes <- function(largeTCRSet, n.max.clonotype, seed, itn=0){
    # Upper limit of iterations
    if(itn>=n.max.iteration){
      message("Iteration limit reached. Final similarity network analysis is attempted...")
      TCRSet_Central_Final <- try(similarityNetworkCentralNodes(largeTCRSet))
      if(class(TCRSet_Central_Final)=="try-error"){
        message("Final similarity network analysis was not successful. Clonotype set is returned.")
        return(largeTCRSet)
      }else{
        message(paste0("Final similarity network analysis is finished. The number of total iterations = ", itn))
        return(TCRSet_Central_Final)
      }
    }
    
    # Declare the next iteration loop
    message(paste0("Iteration: ", itn+1))
    
    # Randomly Split the large TCRSet
    set.seed(seed)
    largeTCRSet.reordered <- sample(largeTCRSet, size=length(largeTCRSet), replace=F)
    N <- ceiling(length(largeTCRSet)/n.max.clonotype)
    TCRSet_Index_List <- parallel::splitIndices(length(largeTCRSet), N)

    # Perform central node analysis to extract important TCRs in each sub-network
    TCRSet_Central_List <- list()
    for(i in 1:N){
      message(paste0("Similarity network analysis: Chunk ", i, "/", N, " is started..."))
      TCRSet_Central_List[[i]] <- similarityNetworkCentralNodes(largeTCRSet.reordered[TCRSet_Index_List[[i]]])$"CentralityDF"$"TCR"
      message(paste0("Similarity network analysis: Chunk ", i, "/", N, " is finished."))
    }
    TCRSet_Central <- unlist(TCRSet_Central_List)

    # Perform central node analysis again with the important TCRs
    if(length(TCRSet_Central)<=n.max.clonotype){
      message("Iteration finished. Final similarity network analysis is started...")
      TCRSet_Central_Final <- similarityNetworkCentralNodes(TCRSet_Central)
      message(paste0("Final similarity network analysis is finished. The number of total iterations = ", itn+1))
      return(TCRSet_Central_Final)
    }else{
      return(list("TCRSet_Central"=TCRSet_Central, "IterationN"=itn+1))
    }
  }

  # Iterative loop of central node extraction
  TCRSet_Central <- list("TCRSet_Central"=largeTCRSet, "IterationN"=0)
  while(all.equal(names(TCRSet_Central),c("TCRSet_Central", "IterationN"))){
    TCRSet_Central <- centralNodes(TCRSet_Central$"TCRSet_Central", n.max.clonotype, seed, TCRSet_Central$"IterationN")
    gc();gc()
    if(all.equal(names(TCRSet_Central),c("Graph", "CentralityDF"))) break
    if(is.character(TCRSet_Central)) break
  }
  return(TCRSet_Central)
}
#' @export
#' @rdname publicTCRAnalysis
#' @name publicTCRAnalysis
publicTCRAnalysis <- function(combinedTCRDF, n.data=2, n.max.clonotype=10000, n.max.iteration=20, 
                              seed=12345, 
                              outDir, header="TCR_"){
  dir.create(outDir, showWarnings=F, recursive=T)
  out1 <- file.path(outDir, paste0(header, "NData", n.data,".csv"))
  out2 <- file.path(outDir, paste0(header, "NData", n.data,".fasta"))
  out3 <- file.path(outDir, paste0(header, "Central_NData", n.data, "_NMaxClono", n.max.clonotype, "_NMaxIte", n.max.iteration, "_Seed", seed, ".csv"))
  out4 <- file.path(outDir, paste0(header, "Central_NData", n.data, "_NMaxClono", n.max.clonotype, "_NMaxIte", n.max.iteration, "_Seed", seed, ".rds"))
  
  # Show the number of clonotypes appeared n.data or more datasets
  TCRSet_N_DF <- cumsum(rev(table(combinedTCRDF$"N")))
  TCRSet_N_DF <- tibble::data_frame("N_Data"=as.numeric(names(TCRSet_N_DF)), "N_Clonotype"=TCRSet_N_DF) %>%
    dplyr::filter(N_Data==n.data)
  print(TCRSet_N_DF)
  TCRSet_N <- dplyr::filter(combinedTCRDF, N>=n.data)
  write.csv(TCRSet_N, out1, row.names=F)
  seqinr::write.fasta(as.list(TCRSet_N$"aaSeqCDR3"), TCRSet_N$"aaSeqCDR3", out2)
  gc();gc()
  
  # Public clonotype analysis
  TCRSet_Central <- similarityNetworkCentralNodes.Large(TCRSet_N$"aaSeqCDR3", n.max.clonotype, n.max.iteration, seed)
  write.csv(TCRSet_Central$"CentralityDF", out3, row.names=F)
  saveRDS(TCRSet_Central, out4)
  gc();gc()
  return(TCRSet_Central)
}
