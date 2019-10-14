#' Parallelized correlation matrix computation.
#'
#' @param mat A matrix, a data.frame, or a data.table.
#' @param nblocks The number of sub-matrices.
#' @param coreN The number of threads for parallelization. Disable by setting \code{NULL}.
#' @export
#' @rdname Utility_ParCor
#' @name Utility_ParCor
parCor <- function(mat, nblocks=10, coreN=parallel:::detectCores(logical=F)){
  ## Add dummy variables if ncol(x) %% nblocks does not give a remainder of 0
  if (ncol(mat) %% nblocks != 0){
    DUMMY <- data.table::as.data.table(matrix(data=0, nrow=nrow(mat), ncol=nblocks-(ncol(mat) %% nblocks)))
    colnames(DUMMY) <- paste0("DUMMY_", 1:ncol(DUMMY))
    x <- cbind(data.table::as.data.table(mat), DUMMY)
  }else{
    x <- data.table::as.data.table(mat)
  }

  ## Split blocks
  NCOL <- ncol(x)
  SPLIT <- split(1:NCOL, rep(1:nblocks, each=NCOL/nblocks))
  COMBS <- data.table::CJ(1:length(SPLIT), 1:length(SPLIT))
  COMBS <- t(apply(COMBS, 1, sort))
  COMBS <- unique(COMBS)

  ## Parallel computation
  cat("Starting parallelized computation of sub-matrices.\n")
  cl <- parallel::makeCluster(coreN, type="SOCK")
  doSNOW::registerDoSNOW(cl)
  sink(tempfile())
  pb <- pbapply::timerProgressBar(max=nrow(COMBS), style=1)
  sink()
  opts <- list(progress=function(n){pbapply::setTimerProgressBar(pb, n)})
  results <- foreach::foreach(i=1:nrow(COMBS), .packages="data.table", .inorder=T, .options.snow=opts)%dopar%{
    COMB <- COMBS[i, ]
    G1 <- SPLIT[[COMB[1]]]
    G2 <- SPLIT[[COMB[2]]]
    cor(x[, ..G1], x[, ..G2])
  }
  cat("\n")
  close(pb)
  parallel::stopCluster(cl)
  gc();gc()

  ## Integrate sub-matrices into one correlation matrix
  corMAT <- matrix(0, ncol=NCOL, nrow=NCOL)
  for(i in 1:nrow(COMBS)){
    COMB <- COMBS[i, ]
    G1 <- SPLIT[[COMB[1]]]
    G2 <- SPLIT[[COMB[2]]]
    COR <- results[[i]]
    corMAT[G1, G2] <- COR
  }
  close(pb)
  corMAT <- Matrix::forceSymmetric(corMAT, uplo="U")
  corMAT <- as.matrix(corMAT)
  colnames(corMAT) <- colnames(x)
  rownames(corMAT) <- colnames(x)
  P <- setdiff(1:ncol(x), grep("^DUMMY_", colnames(x)))
  corMAT <- corMAT[P,P]

  ## Finish
  rm(list=setdiff(ls(), c("corMAT")))
  gc();gc()
  return(corMAT)
}
