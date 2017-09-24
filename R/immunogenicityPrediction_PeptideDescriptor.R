#' Descriptive statistics of peptide fragment descriptors.
#' @param peptideSet A set of peptide sequences.
#' @param fragLen The length of TCR sequence fragments to be matched against the peptide set.
#' @importFrom stringr str_sub
#' @importFrom matrixStats rowMins
#' @importFrom matrixStats rowMaxs
#' @importFrom matrixStats rowMeans2
#' @importFrom matrixStats rowMedians
#' @importFrom dplyr %>%
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom dplyr bind_rows
#' @importFrom tidyr gather
#' @importFrom tidyr spread
#' @importFrom tidyr unite
#' @importFrom parallel splitIndices
#' @importFrom parallel detectCores
#' @importFrom snowfall sfInit
#' @importFrom snowfall sfLibrary
#' @importFrom snowfall sfExport
#' @importFrom snowfall sfLapply
#' @importFrom snowfall sfStop
#' @import Peptides
#' @export
#' @rdname peptideDescriptor
#' @name peptideDescriptor
peptideDescriptorAnalysis <- function(peptideSet, fragLen=5){
  time.start <- proc.time()
  
  # Working functions
  peptideDescriptor.Batch <- function(peptide){
    unlist(list(
      Peptides::aIndex(peptide),
      Peptides::blosumIndices(peptide),
      Peptides::boman(peptide),
      Peptides::charge(peptide),
      Peptides::crucianiProperties(peptide),
      Peptides::fasgaiVectors(peptide),
      Peptides::hmoment(peptide),
      Peptides::hydrophobicity(peptide),
      Peptides::instaIndex(peptide),
      Peptides::kideraFactors(peptide),
      Peptides::mswhimScores(peptide),
      Peptides::pI(peptide),
      Peptides::protFP(peptide),
      Peptides::vhseScales(peptide),
      Peptides::zScales(peptide)
    ))
  }
  peptideDescriptor.NameSet <- c("AliphaticIndex","BLOSUM1","BLOSUM2","BLOSUM3","BLOSUM4","BLOSUM5","BLOSUM6","BLOSUM7","BLOSUM8","BLOSUM9","BLOSUM10","Boman","Charge","PP1","PP2","PP3","F1","F2","F3","F4","F5","F6","HydrophobicMoment","Hydrophobicity","Instability","KF1","KF2","KF3","KF4","KF5","KF6","KF7","KF8","KF9","KF10","MSWHIM1","MSWHIM2","MSWHIM3","pI","ProtFP1","ProtFP2","ProtFP3","ProtFP4","ProtFP5","ProtFP6","ProtFP7","ProtFP8","VHSE1","VHSE2","VHSE3","VHSE4","VHSE5","VHSE6","VHSE7","VHSE8","Z1","Z2","Z3","Z4","Z5")
  peptideDescriptor.FragStat.Single <- function(peptide, fragLen){
    f <- sapply(1:(max(nchar(peptide))-fragLen+1), function(i){stringr::str_sub(peptide, i, i+fragLen-1)})
    d <- sapply(f[nchar(f)==fragLen], function(s){peptideDescriptor.Batch(s)})
    data.frame("Peptide"=peptide,
               "AADescriptor"=peptideDescriptor.NameSet, 
               "Min"=matrixStats::rowMins(d), 
               "Max"=matrixStats::rowMaxs(d), 
               "Mean"=matrixStats::rowMeans2(d),
               "Median"=matrixStats::rowMedians(d)) 
  }
  
  # Parallelization
  snowfall::sfInit(parallel=T, cpus=parallel::detectCores(), type="SOCK")
  sink(tempfile())
  suppressMessages(snowfall::sfLibrary(ff))
  suppressMessages(snowfall::sfLibrary(Peptides))
  suppressMessages(snowfall::sfLibrary(matrixStats))
  sink()
  snowfall::sfExport("peptideDescriptor.FragStat.Single","peptideDescriptor.NameSet","peptideDescriptor.Batch")
  
  # Descriptive statistics of peptide fragment descriptors
  message(paste0("Parallelized fragment descriptor calculation was started. (Memory occupied = ", memory.size(), "[Mb])"))
  df_feature <- dplyr::bind_rows(snowfall::sfLapply(peptideSet, function(p){peptideDescriptor.FragStat.Single(p, fragLen)}))
  message(paste0("Parallelized fragment descriptor calculation was finished. (Memory occupied = ", memory.size(), "[Mb])"))
  gc();gc();
  df_feature <- df_feature %>%
    tidyr::gather(Stat, Value, -Peptide, -AADescriptor) %>%
    tidyr::unite(Feature, AADescriptor, Stat, sep="_") %>%
    dplyr::mutate(Feature=paste0(Feature, "_", fragLen)) %>%
    tidyr::spread(Feature, Value)
  snowfall::sfStop()
  
  # Output
  time.end <- proc.time()
  message(paste0("Overall time required = ", (time.end-time.start)[3], "[sec]"))
  return(df_feature)
}
