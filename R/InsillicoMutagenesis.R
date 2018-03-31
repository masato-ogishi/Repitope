#' In sillico mutagenesis analysis.
#'
#' @param peptideSet A set of peptide sequences.
#' @param peptide.start A peptide sequence.
#' @param peptide.goal A peptide sequence.
#' @importFrom Biostrings AA_STANDARD
#' @importFrom stringdist stringdist
#' @importFrom pbapply pblapply
#' @export
#' @rdname InsillicoMutagenesis
#' @name InsillicoMutagenesis
InSillicoMutagenesis <- function(peptideSet){
  unique(unlist(lapply(
    peptideSet,
    function(pept){
      apply(expand.grid(1:nchar(pept), Biostrings::AA_STANDARD), 1, function(v){substr(pept, v[1], v[1]) <- v[2];return(pept)})
    }
  )))
}
InsillicoMutagenesis_GenerateIntermediates <- function(peptide.start, peptide.goal){
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
