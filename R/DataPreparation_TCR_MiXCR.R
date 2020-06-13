#' MiXCR workflow.
#'
#' \code{MiXCR_Script} generates MiXCR script.\cr
#' \code{MiXCR_CombineClonotypes} reads and combines MiXCR-derived clonotypes. Note: Only TCR beta chains are retained.\cr
#'
#' @param dir A directory in which FASTQ files are stored.
#' @param mixcrPath A path to MiXCR JAR file.
#' @param javaMaxMemory A maximum memory size allocated to the Java virtual machine.
#' @param rnaseq Logical. Specifies if the analysis takes RNASeq datasets.
#' @param species A character string. "hsa" or "mmu".
#' @param dirClones A directory in which MiXCR clonotype TXT files are stored.
#' @param clonesFileNames (Optional) Alternativey, MiXCR clonotype TXT filenames can be provided as a list.
#' @param countSummarise Logical. If True, the number of datasets containing each CDR3 sequence is counted. If False, raw SRA accessions are provided as such.
#' @param rdsFileName (Optional) A filename to which the clonotype summary table is saved as an RDS file.
#' @export
#' @rdname DataPreparation_TCR_MiXCR
#' @name DataPreparation_TCR_MiXCR
MiXCR_Script <- function(
  dir="C:/SRA/output",
  mixcrPath="C:/MiXCR/",
  javaMaxMemory="60G",
  rnaseq=F,
  species="hsa"
){
  # MiXCR script template [type 1: single-end, type 2: paired-end]
  mixcrJarPath <- normalizePath(file.path(mixcrPath, "mixcr.jar"))
  template <- function(accessionNumber, type, rnaseq){
    outDir <- suppressWarnings(normalizePath(file.path(dir, accessionNumber)))
    if(rnaseq==F){
      if(type==1){
        s <- paste0(c(
          paste0("mkdir ", outDir),
          paste0("java -jar -Xmx", javaMaxMemory, " ", mixcrJarPath, " align -s ", species, " -r ", outDir, "\\log.txt ", outDir, "_pass.fastq.gz ", outDir, "\\alignments.vdjca"),
          paste0("java -jar -Xmx", javaMaxMemory, " ", mixcrJarPath, " assemble -r ", outDir, "\\log.txt ", outDir, "\\alignments.vdjca ", outDir, "\\clones.clns"),
          paste0("java -jar -Xmx", javaMaxMemory, " ", mixcrJarPath, " exportClones ", outDir, "\\clones.clns ", outDir, "\\clones.txt")),
          collapse = "\n"
        )
      }
      if(type==2){
        s <- paste0(c(
          paste0("mkdir ", outDir),
          paste0("java -jar -Xmx", javaMaxMemory, " ", mixcrJarPath, " align -s ", species, " -r ", outDir, "\\log.txt ", outDir, "_pass_1.fastq.gz ", outDir, "_pass_2.fastq.gz ", outDir, "\\alignments.vdjca"),
          paste0("java -jar -Xmx", javaMaxMemory, " ", mixcrJarPath, " assemble -r ", outDir, "\\log.txt ", outDir, "\\alignments.vdjca ", outDir, "\\clones.clns"),
          paste0("java -jar -Xmx", javaMaxMemory, " ", mixcrJarPath, " exportClones ", outDir, "\\clones.clns ", outDir, "\\clones.txt")),
          collapse = "\n"
        )
      }
    }else{
      ## RNASeq-specific workflow
      ## https://mixcr.readthedocs.io/en/master/rnaseq.html
      if(type==1){
        s <- paste0(c(
          paste0("mkdir ", outDir),
          paste0("java -jar -Xmx", javaMaxMemory, " ", mixcrJarPath, " align -p rna-seq -s ", species, " -r ", outDir, "\\log.txt ", outDir, "_pass.fastq.gz ", outDir, "\\alignments.vdjca"),
          paste0("java -jar -Xmx", javaMaxMemory, " ", mixcrJarPath, " assemblePartial -r ", outDir, "\\log.txt ", outDir, "\\alignments.vdjca ", outDir, "\\alignments_rescued_1.vdjca"),
          paste0("java -jar -Xmx", javaMaxMemory, " ", mixcrJarPath, " assemblePartial -r ", outDir, "\\log.txt ", outDir, "\\alignments_rescued_1.vdjca ", outDir, "\\alignments_rescued_2.vdjca"),
          paste0("java -jar -Xmx", javaMaxMemory, " ", mixcrJarPath, " extend -r ", outDir, "\\log.txt ", outDir, "\\alignments_rescued_2.vdjca ", outDir, "\\alignments_extended.vdjca"),
          paste0("java -jar -Xmx", javaMaxMemory, " ", mixcrJarPath, " assemble -r ", outDir, "\\log.txt ", outDir, "\\alignments_extended.vdjca ", outDir, "\\clones.clns"),
          paste0("java -jar -Xmx", javaMaxMemory, " ", mixcrJarPath, " exportClones ", outDir, "\\clones.clns ", outDir, "\\clones.txt")),
          collapse = "\n"
        )
      }
      if(type==2){
        s <- paste0(c(
          paste0("mkdir ", outDir),
          paste0("java -jar -Xmx", javaMaxMemory, " ", mixcrJarPath, " align -p rna-seq -s ", species, " -r ", outDir, "\\log.txt ", outDir, "_pass_1.fastq.gz ", outDir, "_pass_2.fastq.gz ", outDir, "\\alignments.vdjca"),
          paste0("java -jar -Xmx", javaMaxMemory, " ", mixcrJarPath, " assemblePartial -r ", outDir, "\\log.txt ", outDir, "\\alignments.vdjca ", outDir, "\\alignments_rescued_1.vdjca"),
          paste0("java -jar -Xmx", javaMaxMemory, " ", mixcrJarPath, " assemblePartial -r ", outDir, "\\log.txt ", outDir, "\\alignments_rescued_1.vdjca ", outDir, "\\alignments_rescued_2.vdjca"),
          paste0("java -jar -Xmx", javaMaxMemory, " ", mixcrJarPath, " extend -r ", outDir, "\\log.txt ", outDir, "\\alignments_rescued_2.vdjca ", outDir, "\\alignments_extended.vdjca"),
          paste0("java -jar -Xmx", javaMaxMemory, " ", mixcrJarPath, " assemble -r ", outDir, "\\log.txt ", outDir, "\\alignments.vdjca ", outDir, "\\clones.clns"),
          paste0("java -jar -Xmx", javaMaxMemory, " ", mixcrJarPath, " exportClones ", outDir, "\\clones.clns ", outDir, "\\clones.txt")),
          collapse = "\n"
        )
      }
    }
    return(s)
  }

  # MiXCR script
  f1 <- list.files(path=dir, pattern=".fastq.gz")
  if(!length(f1)==0){
    d1 <- stringr::str_split(f1, "_", simplify=T) %>%
      as.data.frame() %>%
      dplyr::group_by(V1) %>%
      dplyr::summarise(Count=dplyr::n())
    mixcr.script <- paste0(mapply(template, d1$V1, d1$Count, MoreArgs=list(rnaseq=rnaseq)), collapse="\n")
    write.table(mixcr.script, file.path(dir, "MiXCR_script.txt"), row.names=F, col.names=F, quote=F)
  }else{
    message("No FastQ files.")
    return(NULL)
  }

  # MiXCR script [existing 'clones.txt' files omitted]
  f2 <- list.files(path=dir, pattern="clones.txt", full.names=T, recursive=T)
  f2 <- gsub(normalizePath(paste0(dir, "/"), "/"), "", f2, fixed=T)
  if(!length(f2)==0){
    d2 <- stringr::str_split(f2, "/", simplify=T) %>%
      as.data.frame() %>%
      dplyr::group_by(V1) %>%
      dplyr::distinct(V1)
    d3 <- d1 %>%
      dplyr::filter(V1 %in% setdiff(d1$V1, d2$V1))
    if(purrr::is_empty(d3)){message("No MiXCR analysis remaining to be done.")}
    mixcr.script.rem <- paste0(mapply(template, d3$V1, d3$Count, MoreArgs=list(rnaseq=rnaseq)), collapse="\n")
    write.table(mixcr.script.rem, file.path(dir, "MiXCR_script_remaining.txt"), row.names=F, col.names=F, quote=F)
  }else{
    message("No existing MiXCR folders detected.")
    return(file.path(dir, "MiXCR_script.txt"))
  }
  return(list(
    file.path(dir, "MiXCR_script.txt"),
    file.path(dir, "MiXCR_script_remaining.txt")
  ))
}

#' @export
#' @rdname DataPreparation_TCR_MiXCR
#' @name DataPreparation_TCR_MiXCR
MiXCR_CombineClonotypes <- function(
  dirClones=NULL,
  clonesFileNames=NULL,
  countSummarise=F,
  rdsFileName=NULL
){
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
