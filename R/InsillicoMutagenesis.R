#' In sillico mutagenesis analysis.
#'
#' \code{InSillicoMutagenesis} generates single-mutated peptide sequences.
#' \code{InsillicoMutagenesis_GenerateIntermediates} generates mutational intermediates between two peptide sequences.
#' \code{InsillicoMutagenesis_PairwiseSimilarityNetwork} generates single-mutated peptide sequences and construct similarity networks.
#'
#' @param peptideSet A set of peptide sequences.
#' @param peptidePair A pair of peptide sequences.
#' @param peptidePairDF A dataframe which has two columns of paired peptide sequences.
#' @param union Logical. Whether the output list should be combined.
#' @param coreN The number of threads.
#' @importFrom Biostrings AA_STANDARD
#' @importFrom stringdist stringdist
#' @importFrom Matrix t
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom igraph set_vertex_attr
#' @importFrom igraph union
#' @importFrom igraph simplify
#' @importFrom pbapply pblapply
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#' @importFrom parallel clusterExport
#' @export
#' @rdname InsillicoMutagenesis
#' @name InsillicoMutagenesis
InSillicoMutagenesis <- function(peptideSet, union=T, coreN=NULL){
  if(is.null(coreN)){
    cl <- NULL
  }else{
    cl <- parallel::makeCluster(coreN)
  }
  mut <- pbapply::pblapply(
    peptideSet,
    function(pept){
      apply(expand.grid(1:nchar(pept), Biostrings::AA_STANDARD), 1, function(v){substr(pept, v[1], v[1]) <- v[2];return(pept)})
    },
    cl=cl
  )
  if(!is.null(coreN)) parallel::stopCluster(cl=cl)
  gc();gc()
  if(union==T) mut <- unique(unlist(mut))
  return(mut)
}

#' @export
#' @rdname InsillicoMutagenesis
#' @name InsillicoMutagenesis
InsillicoMutagenesis_GenerateIntermediates <- function(peptidePair){
  # Peptides
  peptide.start <- peptidePair[[1]]
  peptide.goal <- peptidePair[[2]]

  # The minimum levenstein distance
  minDist <- stringdist::stringdist(peptide.start, peptide.goal, method="lv")
  cat("The minimum number of required substitution: ", minDist, "\n", sep="")

  # In sillico mutagenesis
  Nth_Distant_Peptides <- function(peptideSet, peptide.goal, nth){
    peptide.mut <- inSillicoMutagenesis(peptideSet)
    d <- stringdist::stringdist(peptide.mut, peptide.goal, method="lv")
    peptide.mut[which(d==nth)]
  }
  peptide.mut.list <- as.list(numeric(minDist))
  peptide.mut.list[[1]] <- peptide.start
  peptide.mut.list[[minDist+1]] <- peptide.goal
  if(minDist<=1) return(unlist(peptide.mut.list))
  for(i in seq(1, minDist-1)){
    peptide.mut.list[[i+1]] <- Nth_Distant_Peptides(peptide.mut.list[[i]], peptide.goal, nth=minDist-i)
  }
  return(unique(unlist(peptide.mut.list)))
}

#' @export
#' @rdname InsillicoMutagenesis
#' @name InsillicoMutagenesis
InsillicoMutagenesis_PairwiseSimilarityNetwork <- function(peptidePairDF, coreN=NULL){
  simnet_singlepair <- function(pept1, pept2){
    peptideSet.mut <- Repitope::InSillicoMutagenesis(c(pept1, pept2))
    d <- Repitope::distMat_Auto(peptideSet.mut)
    d <- d|Matrix::t(d)
    g <- igraph::graph_from_adjacency_matrix(d, mode="undirected", weighted=NULL, diag=F)
    g <- igraph::set_vertex_attr(g, "name", value=peptideSet.mut)
    return(g)
  }
  if(is.null(coreN)){
    cl <- NULL
  }else{
    cl <- parallel::makeCluster(coreN)
    parallel::clusterExport(cl=cl, list=c("peptidePairDF","simnet_singlepair"), envir=environment())
  }
  simNetList <- pbapply::pblapply(
    1:nrow(peptidePairDF),
    function(i){simnet_singlepair(peptidePairDF[i,1], peptidePairDF[i,2])},
    cl=cl
  )
  if(!is.null(coreN)) parallel::stopCluster(cl=cl)
  gc();gc()
  return(simNetList)
}
