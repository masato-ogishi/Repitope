#' Filter amino acid sequences to exclude those containing non-standard letters.
#' @param sequenceSet A set of amino acid sequences.
#' @importFrom stringr str_split
#' @importFrom Biostrings AA_STANDARD
#' @importFrom stringr str_detect
#' @importFrom stringr fixed
#' @export
#' @rdname Utility_SequenceIO
#' @name Utility_SequenceIO
sequenceFilter <- function(sequenceSet){
  s <- sequenceSet[!is.na(sequenceSet)]
  s <- toupper(s)
  letters <- unique(unlist(stringr::str_split(s, "")))
  letters.exclude <- setdiff(letters, Biostrings::AA_STANDARD)
  for(l in letters.exclude){
    s <- s[!stringr::str_detect(s, stringr::fixed(l))]
  }
  s
}
