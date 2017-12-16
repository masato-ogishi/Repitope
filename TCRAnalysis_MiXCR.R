#' MiXCR workflow.
#' 
#' \code{MiXCRScript} generates MiXCR script.\cr
#' \code{MiXCRClonotypeCombine} reads and combines MiXCR-derived clonotypes. Note: Only TCR beta chains are retained.\cr
#' 
#' @param dir A directory in which FASTQ files are stored.
#' @param mixcrPath A path to MiXCR JAR file.
#' @param javaMaxMemory A maximum memory size allocated to the Java virtual machine.
#' @param dirClones A directory in which MiXCR clonotype TXT files are stored.
#' @param clonesFileNames (Optional) Alternativey, MiXCR clonotype TXT filenames can be provided as a list.
#' @param countSummarise Logical. If True, the number of datasets containing each CDR3 sequence is counted. If False, raw SRA accessions are provided as such.
#' @param rdsFileName (Optional) A filename to which the clonotype summary table is saved as an RDS file.
#' @importFrom stringr str_split
#' @importFrom stringr str_detect
#' @importFrom stringr fixed
#' @importFrom dplyr %>%
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr distinct
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom dplyr mutate
#' @importFrom dplyr left_join
#' @importFrom dplyr one_of
#' @importFrom purrr is_empty
#' @importFrom magrittr set_colnames
#' @importFrom pbapply pblapply
#' @importFrom data.table fread
#' @importFrom data.table rbindlist
#' @importFrom Biostrings AA_STANDARD
#' @export
#' @rdname TCRAnalysis_MiXCR
#' @name TCRAnalysis_MiXCR
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
    if(purrr::is_empty(d3)){message("No MiXCR analysis remaining to be done.")}
    mixcr.script.rem <- paste0(mapply(template, d3$V1, d3$Count), collapse="\n")
    write.table(mixcr.script.rem, file.path(dir, "MiXCR_script_remaining.txt"), row.names=F, col.names=F, quote=F)
  }else{
    message("No existing MiXCR folders detected.")
    return(file.path(dir, "MiXCR_script.txt"))
  }
  return(list(file.path(dir, "MiXCR_script.txt"), file.path(dir, "MiXCR_script_remaining.txt")))
}
#' @export
#' @rdname TCRAnalysis_MiXCR
#' @name TCRAnalysis_MiXCR
MiXCRClonotypeCombine <- function(dirClones=NULL, clonesFileNames=NULL, countSummarise=F, rdsFileName=NULL){
  TRBClonotypeTableImport <- function(FILEPATH){
    DT <- suppressWarnings(data.table::fread(FILEPATH, sep="\t", data.table=T))
    DT <- subset(DT, stringr::str_detect(allVHitsWithScore, "TRB"))[,list(aaSeqCDR3),]
    s <- DT[["aaSeqCDR3"]]
    s <- s[!is.na(s)]
    s <- toupper(s)
    letters <- unique(unlist(stringr::str_split(s, "")))
    letters.exclude <- setdiff(letters, Biostrings::AA_STANDARD)
    for(l in letters.exclude){
      DT <- subset(DT, !stringr::str_detect(aaSeqCDR3, stringr::fixed(l)))
    }
    if(nrow(DT)==0){ return(DT) }
    DT[["Accession"]] <- unlist(stringr::str_split(FILEPATH, "/"))[[4]]
    return(DT)
  }
  if(!is.null(dirClones)){
    tcrFiles <- list.files(path=dirClones, pattern="clones.txt", full.names=T, recursive=T)
  }
  if(!is.null(clonesFileNames)){
    tcrFiles <- clonesFileNames
  }
  tcrDT <- pbapply::pblapply(tcrFiles, TRBClonotypeTableImport)
  tcrDT <- data.table::rbindlist(tcrDT[sapply(tcrDT, ncol)==2])
  if(countSummarise==T){
    tcrDT <- tcrDT[,length(unique(Accession)), by="aaSeqCDR3"]
    colnames(tcrDT) <- c("aaSeqCDR3", "N")
  }
  if(!is.null(rdsFileName)){
    saveRDS(tcrDT, rdsFileName)
  }
  return(tcrDT)
}
