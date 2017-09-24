#' MiXCR workflow.
#' 
#' \code{MiXCRScript} generates MiXCR script.\cr
#' \code{MiXCRClonotypeCombine} reads and combines MiXCR-derived clonotypes. Note: Only TCR beta chains are retained.
#' 
#' @param dir A directory in which FASTQ files are stored.
#' @param mixcrPath A path to MiXCR JAR file.
#' @param javaMaxMemory A maximum memory size allocated to the Java virtual machine.
#' @param outDir A directory in which the clonotype summary table is saved.
#' @param rdsFileName A filename for the clonotype summary.
#' @param countSummarise A logical. If True, the number of datasets containing each CDR3 sequence is counted. If False, raw SRA accessions are provided as such.
#' @importFrom stringr str_split
#' @importFrom stringr str_detect
#' @importFrom dplyr %>%
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr distinct
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom dplyr mutate
#' @importFrom dplyr bind_rows
#' @importFrom tidyr separate
#' @importFrom readr read_csv
#' @importFrom magrittr set_colnames
#' @importFrom dplyr left_join
#' @importFrom dplyr one_of
#' @export
#' @rdname MiXCR
#' @name MiXCR
MiXCRScript <- function(dir="C:/SRA/output", mixcrPath="C:/MiXCR/", javaMaxMemory="6G"){
  # MiXCR script template [type 1: single-end, type 2: paired-end]
  mixcrJarPath <- normalizePath(file.path(mixcrPath, "mixcr.jar"))
  template <- function(accessionNumber, type){
    outDir <- suppressWarnings(normalizePath(file.path(dir, accessionNumber)))
    if(type==1){
      s <- paste0(c(
        paste0("mkdir ", outDir),
        paste0("java -jar -Xmx", javaMaxMemory, " ", mixcrJarPath, " align -r ", outDir, "\\log.txt ", outDir, "_pass_1.fastq.gz ", outDir, "\\alignments.vdjca"),
        paste0("java -jar -Xmx", javaMaxMemory, " ", mixcrJarPath, " assemble -r ", outDir, "\\log.txt ", outDir, "\\alignments.vdjca ", outDir, "\\clones.clns"),
        paste0("java -jar -Xmx", javaMaxMemory, " ", mixcrJarPath, " exportClones ", outDir, "\\clones.clns ", outDir, "\\clones.txt")),
        collapse = "\n"
      )
    }
    if(type==2){
      s <- paste0(c(
        paste0("mkdir ", outDir),
        paste0("java -jar -Xmx", javaMaxMemory, " ", mixcrJarPath, " align -r ", outDir, "\\log.txt ", outDir, "_pass_1.fastq.gz ", outDir, "_pass_2.fastq.gz ", outDir, "\\alignments.vdjca"),
        paste0("java -jar -Xmx", javaMaxMemory, " ", mixcrJarPath, " assemble -r ", outDir, "\\log.txt ", outDir, "\\alignments.vdjca ", outDir, "\\clones.clns"),
        paste0("java -jar -Xmx", javaMaxMemory, " ", mixcrJarPath, " exportClones ", outDir, "\\clones.clns ", outDir, "\\clones.txt")),
        collapse = "\n"
      )
    }
    return(s)
  }
  
  # MiXCR script
  f1 <- list.files(path=dir, pattern=".fastq.gz")
  if(!length(f1)==0){
    d1 <- stringr::str_split(f1, "_", simplify=T) %>% as.data.frame() %>% dplyr::group_by(V1) %>% dplyr::summarise(Count=n())
    mixcr.script <- paste0(mapply(template, d1$V1, d1$Count), collapse="\n")
    write.table(mixcr.script, file.path(dir, "MiXCR_script.txt"), row.names=F, col.names=F, quote=F)
  }else{
    message("No FastQ files.")
    return(NULL)
  }
  
  # MiXCR script [existing 'clones.txt' files omitted]
  f2 <- list.files(path=dir, pattern="clones.txt", full.names=T, recursive=T)
  f2 <- gsub(normalizePath(paste0(dir, "/"), "/"), "", f2, fixed=T)
  if(!length(f2)==0){
    d2 <- stringr::str_split(f2, "/", simplify=T) %>% as.data.frame() %>% dplyr::group_by(V1) %>% dplyr::distinct(V1)
    d3 <- d1 %>% dplyr::filter(V1 %in% setdiff(d1$V1, d2$V1))
    mixcr.script.rem <- paste0(mapply(template, d3$V1, d3$Count), collapse="\n")
    write.table(mixcr.script.rem, file.path(dir, "MiXCR_script_remaining.txt"), row.names=F, col.names=F, quote=F)
  }else{
    message("No existing MiXCR folders detected.")
    return(file.path(dir, "MiXCR_script.txt"))
  }
  return(list(file.path(dir, "MiXCR_script.txt"), file.path(dir, "MiXCR_script_remaining.txt")))
}
#' @export
#' @rdname MiXCR
#' @name MiXCR
MiXCRClonotypeCombine <- function(dir, outDir, rdsFileName="TCR_Combined.rds", countSummarise=F){
  tcrDF <- dplyr::bind_rows(
    lapply(list.files(path=dir, pattern="clones.txt", full.names=T, recursive=T), 
           function(f){
             read.delim(f) %>% 
               dplyr::filter(stringr::str_detect(allVHitsWithScore, "TRB")) %>% 
               dplyr::select(cloneFraction, aaSeqCDR3) %>% 
               dplyr::mutate(FilePath=f)
           })
  ) %>%
    dplyr::filter(aaSeqCDR3 %in% sequenceFilter(aaSeqCDR3)) %>%
    tidyr::separate(FilePath, c("D1", "D2", "D3", "Accession"), sep="/") %>%
    dplyr::select(aaSeqCDR3, Accession)
  if(countSummarise==T){
    tcrDF <- tcrDF %>% dplyr::group_by(aaSeqCDR3) %>% dplyr::summarise(N=n())
  }
  saveRDS(tcrDF, file.path(outDir, rdsFileName))
  return(tcrDF)
}
