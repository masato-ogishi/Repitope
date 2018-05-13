#' In sillico mutagenesis analysis.
#'
#' \code{InSillicoMutagenesis} generates single-mutated peptide sequences.
#' \code{InSillicoMutagenesis_GenerateIntermediates} generates mutational intermediates between two peptide sequences.
#'
#' @param peptideSet A set of peptide sequences.
#' @param peptidePair A pair of peptide sequences.
#' @param coreN The number of threads for parallelization.
#' @export
#' @rdname InSillicoMutagenesis
#' @name InSillicoMutagenesis
InSillicoMutagenesis <- function(peptideSet, coreN=parallel::detectCores()){
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
  mut <- unique(unlist(mut))
  return(mut)
}

#' @export
#' @rdname InSillicoMutagenesis
#' @name InSillicoMutagenesis
InSillicoMutagenesis_GenerateIntermediates <- function(peptidePair, coreN=parallel::detectCores()){
  # Peptides
  peptide.start <- peptidePair[[1]]
  peptide.goal <- peptidePair[[2]]

  # The minimum levenstein distance
  minDist <- stringdist::stringdist(peptide.start, peptide.goal, method="lv")
  cat("The minimum number of required substitution: ", minDist, "\n", sep="")

  # In sillico mutagenesis
  Nth_Distant_Peptides <- function(peptideSet, peptide.goal, nth){
    peptide.mut <- InSillicoMutagenesis(peptideSet, coreN=coreN)
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
  mut <- unique(unlist(peptide.mut.list))
  return(mut)
}

