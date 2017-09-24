#' Downloading fastQ files from NCBI Sequence Read Archive (SRA).
#' 
#' \code{fastqDumpScript} generates fustq-dump script.
#' 
#' @param accessionSet A set of SRA Run Accession numbers.
#' @param outDir A directory for outputs.
#' @importFrom stringr str_split
#' @importFrom dplyr group_by
#' @importFrom dplyr distinct
#' @export
#' @rdname SRA
fastqDumpScript <- function(accessionSet, outDir="C:/SRA/output"){
  dir.create(outDir, showWarnings=F, recursive=T)
  fastq.dump.script <- paste0( 
    sapply(sort(accessionSet), function(a){paste0("C:\\SRA\\sratoolkit.2.8.1-3-win64\\bin\\fastq-dump.exe --outdir ", normalizePath(outDir), " --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-files --clip --accession ", a)}), 
    collapse="\n"
  )
  write.table(fastq.dump.script, file.path(outDir, "fastq-dump_script.txt"), row.names=F, col.names=F, quote=F)
  f <- list.files(path=outDir, pattern="fastq.gz", full.names=T, recursive=T)
  f <- gsub(normalizePath(paste0(outDir, "/"), "/"), "", f, fixed=T)
  if(!length(f)==0){
    d <- stringr::str_split(f, "_", simplify=T) %>% as.data.frame() %>% dplyr::group_by(V1) %>% dplyr::distinct(V1)
    fastq.dump.script.rem <- paste0( 
      sapply(sort(setdiff(accessionSet, d$V1)), function(a){paste0("C:\\SRA\\sratoolkit.2.8.1-3-win64\\bin\\fastq-dump.exe --outdir ", normalizePath(outDir), " --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-files --clip --accession ", a)}), 
      collapse="\n"
    )
    write.table(fastq.dump.script.rem, file.path(outDir, "fastq-dump_script_remaining.txt"), row.names=F, col.names=F, quote=F)
  }else{
    message("No existing FastQ files detected.")
    return(file.path(outDir, "fastq-dump_script.txt"))
  }
  return(list(file.path(outDir, "fastq-dump_script.txt"), file.path(outDir, "fastq-dump_script_remaining.txt")))
  ## https://edwards.sdsu.edu/research/fastq-dump/
}
