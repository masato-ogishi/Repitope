#' Neighbor network analysis.
#'
#' @param peptideSet A set of peptide sequences.
#' @param numSet An attribute for the vertices.
#' @param directed Should the network be converted from undirected to directed? Directions are determined using the \code{numSet} provided.
#' @param weighted Should the network be converted to weihted? Edge weights are determined using the \code{numSet} provided.
#' @param coreN The number of cores to be used for parallelization.
#' @export
#' @rdname NeighborNetwork
#' @name NeighborNetwork
neighborNetwork <- function(
  peptideSet,
  numSet=NULL,
  directed=F,
  weighted=F,
  coreN=parallel::detectCores(logical=F)
){
  # Start calculation
  set.seed(12345)
  time.start <- proc.time()

  # Explore neighbor pairs
  cl <- parallel::makeCluster(coreN, type="PSOCK")
  doSNOW::registerDoSNOW(cl)
  peptideList <- split(peptideSet, nchar(peptideSet))
  cat("Searching neighbors by single-aa substitutions...\n")
  dt_neighbor_sub <- foreach::foreach(i=1:length(peptideList), .inorder=F)%dopar%{
    l <- nchar(peptideList[[i]][1])
    dt_edges <- apply(data.table::CJ(1:l, Biostrings::AA_STANDARD), 1,
                      function(v){
                        aa_pos <- as.numeric(v[1])
                        aa_sub <- v[2]
                        s <- peptideList[[i]]
                        s_prime <- s
                        stringr::str_sub(s_prime, start=aa_pos, end=aa_pos) <- aa_sub
                        pairPos <- setdiff(which(s_prime %in% s), which(s_prime==s))
                        dt_edges <- data.table::data.table(
                          "AASeq1"=s,
                          "AASeq2"=s_prime,
                          "AA1"=stringr::str_sub(s, start=aa_pos, end=aa_pos),
                          "MutPosition"=aa_pos,
                          "AA2"=aa_sub,
                          "MutPattern"="Substitution"
                        )
                        dt_edges <- dt_edges[pairPos,]
                        dt_edges[,MutType:=paste0(AA1, "_", MutPosition, "_", AA2)][,AA1:=NULL][,MutPosition:=NULL][,AA2:=NULL]
                        return(dt_edges)
                      })
    data.table::rbindlist(dt_edges)
  } %>%
    data.table::rbindlist() %>%
    unique(fromLast=F, by=c("AASeq1","AASeq2","MutPattern"))
  cat("Searching neighbors by single-aa insertions...\n")
  dt_neighbor_ins <- foreach::foreach(i=1:(length(peptideList)-1), .inorder=F)%dopar%{
    l <- nchar(peptideList[[i]][1])
    dt_edges <- apply(data.table::CJ(1:l, Biostrings::AA_STANDARD), 1,
                      function(v){
                        aa_pos <- as.numeric(v[1])
                        aa_sub <- v[2]
                        s <- peptideList[[i]]
                        s_prime <- s
                        s_long <- peptideList[[i+1]]
                        stringi::stri_sub(s_prime, aa_pos, aa_pos-1) <- aa_sub
                        dt_edges <- data.table::data.table(
                          "AASeq1"=s,
                          "AASeq2"=s_prime,
                          "AA1"="-",
                          "MutPosition"=aa_pos,
                          "AA2"=aa_sub,
                          "MutPattern"="Insertion"
                        )
                        dt_edges <- dt_edges[AASeq2 %in% s_long,]
                        dt_edges[,MutType:=paste0(AA1, "_", MutPosition, "_", AA2)][,AA1:=NULL][,MutPosition:=NULL][,AA2:=NULL]
                        return(dt_edges)
                      })
    data.table::rbindlist(dt_edges)
  } %>%
    data.table::rbindlist() %>%
    unique(fromLast=F, by=c("AASeq1","AASeq2","MutPattern"))
  cat("Searching neighbors by single-aa deletions...\n")
  dt_neighbor_del <- foreach::foreach(i=2:length(peptideList), .inorder=F)%dopar%{
    l <- nchar(peptideList[[i]][1])
    dt_edges <- lapply(1:l, function(aa_pos){
      s <- peptideList[[i]]
      s_prime <- s
      s_short <- peptideList[[i-1]]
      stringi::stri_sub(s_prime, aa_pos, aa_pos) <- ""
      dt_edges <- data.table::data.table(
        "AASeq1"=s,
        "AASeq2"=s_prime,
        "AA1"=stringr::str_sub(s, start=aa_pos, end=aa_pos),
        "MutPosition"=aa_pos,
        "AA2"="-",
        "MutPattern"="Deletion"
      )
      dt_edges <- dt_edges[AASeq2 %in% s_short,]
      dt_edges[,MutType:=paste0(AA1, "_", MutPosition, "_", AA2)][,AA1:=NULL][,MutPosition:=NULL][,AA2:=NULL]
      return(dt_edges)
    })
    data.table::rbindlist(dt_edges)
  } %>%
    data.table::rbindlist() %>%
    unique(fromLast=F, by=c("AASeq1","AASeq2","MutPattern"))
  parallel::stopCluster(cl)
  gc();gc()
  cat("Merging...\n")
  dt_neighbor <- rbind(
    dt_neighbor_sub, dt_neighbor_ins, dt_neighbor_del
  )
  dt_neighbor[,MutPattern:=factor(MutPattern, levels=c("Substitution","Insertion","Deletion"))]
  data.table::setorder(dt_neighbor, AASeq1, AASeq2, MutPattern, MutType)

  # Construct a neighbor network object
  net_neighbor <- igraph::graph_from_data_frame(dt_neighbor, directed=F)
  net_neighbor <- igraph::set_vertex_attr(net_neighbor, name="label", value=igraph::V(net_neighbor)$name)
  if(is.null(numSet)){
    out <- list("NeighborNetwork"=net_neighbor, "PairDF"=dt_neighbor)
  }else{
    ## Weighted and directed network using a given set of scores
    ## Direction is determined such that the score always increases after mutation
    dt_num <- data.table::data.table("Peptide"=peptideSet, "Score"=numSet)
    dt_neighbor_dw <- data.table::copy(dt_neighbor)
    dt_neighbor_dw <- merge(dt_neighbor_dw, dt_num, by.x="AASeq1", by.y="Peptide")
    dt_neighbor_dw <- merge(dt_neighbor_dw, dt_num, by.x="AASeq2", by.y="Peptide", suffixes=c("1","2"))
    dt_neighbor_dw <- dt_neighbor_dw[Score2>=Score1,]
    dt_neighbor_dw[,ScoreRatio:=Score2/Score1]
    dt_neighbor_dw[,EdgeWeight:=Score1/Score2]
    dt_neighbor_dw[EdgeWeight<0.001, EdgeWeight:=0.001] ### set an arbitrary minimum threshold
    net_neighbor_dw <- igraph::graph_from_data_frame(dt_neighbor_dw, directed=directed)
    net_neighbor_dw <- igraph::set_vertex_attr(net_neighbor_dw, name="label", value=igraph::V(net_neighbor_dw)$name)
    if(weighted==T) igraph::E(net_neighbor_dw)$weight <- igraph::E(net_neighbor_dw)$EdgeWeight
    out <- list(
      "NeighborNetwork"=net_neighbor, "PairDF"=dt_neighbor,
      "NeighborNetwork_DW"=net_neighbor_dw, "PairDF_DW"=dt_neighbor_dw
    )
  }

  # Finish the timer
  time.end <- proc.time()
  message("Overall time required = ", format((time.end-time.start)[3], nsmall=2), "[sec]")

  # Output
  gc();gc()
  return(out)
}

