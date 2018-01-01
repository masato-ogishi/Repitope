#' Parallelized correlation matrix calculation.
#'
#' @param mat A matrix. A data.frame or a data.table is internally converted into a matrix.
#' @param num_splits The number of chunk blocks. Default setting is the number of core (for parallelization).
#' @param verbose Either 0, 1, 2, or TRUE/FALSE.
#' @importFrom parallel detectCores
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterExport
#' @importFrom parallel stopCluster
#' @importFrom pbapply pblapply
#' @importFrom Matrix forceSymmetric
#' @importFrom bigmemory as.big.matrix
#' @importFrom bigmemory attach.big.matrix
#' @export
#' @rdname Utility_ParCor
#' @name Utility_ParCor
parCor <- function(mat, num_splits=parallel::detectCores(), verbose=2){
  # Input check
  if(!any(class(mat)=="matrix")){
    mat <- try(as.matrix(mat), silent=T)
    if(identical(class(mat), "try-error")){
      message("The input mat cannot be converted into matrix!")
      return(NULL)
    }
  }

  # Convert to big.matrix
  tmp.timestamp <- format(Sys.time(), "%Y.%b.%d.%H.%M.%OS3")
  tmp.dir <- tempdir()
  mat.binfile <- paste0("ParCor.", tmp.timestamp, ".bin")
  mat.descfile <- paste0("ParCor.", tmp.timestamp, ".desc")
  mat <- bigmemory::as.big.matrix(
    mat, type = "double", separated = F,
    backingpath = tmp.dir,
    backingfile = mat.binfile,
    descriptorfile = mat.descfile
  )

  # Internally used functions
  compute_split_size <- function(n, num_splits){
    split_size <- ceiling(n / num_splits)
    stopifnot(split_size > 1)
    return(split_size)
  }
  compute_splits_beg <- function(n, split_size){
    beg <- seq(1, n, by = split_size)
    beg <- beg[beg <= n]
    return(beg)
  }
  compute_splits_end <- function(n, split_size, beg){
    end <- beg - 1 + split_size
    end[length(end)] <- n
    return(end)
  }

  # Parallelized computation of correlation matrix

  ### variables
  n <- nrow(mat)
  p <- ncol(mat)

  ### split into batches
  if(verbose>0|verbose) cat(" - bigcrossprod: preparing batches of block-pairs from ", num_splits, "splits\n")
  split_size <- compute_split_size(p, num_splits)
  beg <- compute_splits_beg(p, split_size)
  end <- compute_splits_end(p, split_size, beg)
  num_splits <- length(beg)
  batches <- lapply(1:num_splits, function(i){
    list(batch = i, beg = beg[i], end = end[i])
  })
  grid <-  expand.grid(1:num_splits, 1:num_splits)
  colnames(grid) <- c("cell1", "cell2")
  grid <- subset(grid, cell1 <= cell2)
  num_batches <- nrow(grid)

  ### computes blocks of the cov. matrix
  if(verbose>0|verbose) cat(" - bigcrossprod: computing `crossprod` by blocks:", num_batches, "batches\n")
  cl <- parallel::makeCluster(min(num_splits, parallel::detectCores()), type="SOCK")
  parallel::clusterExport(
    cl,
    list("grid", "batches", "beg", "end", "mat.descfile", "tmp.dir"),
    envir=environment()
  )
  blocks <- pbapply::pblapply(1:num_batches, function(i){
    tmpmat <- bigmemory::attach.big.matrix(obj=mat.descfile, path=tmp.dir)
    cell1 <- grid[i, "cell1"]
    cell2 <- grid[i, "cell2"]
    ind1 <- with(batches[[cell1]], seq(beg, end))
    ind2 <- with(batches[[cell2]], seq(beg, end))
    if(cell1==cell2){
      cp <- crossprod(tmpmat[, ind1])
    }else{
      cp <- crossprod(tmpmat[, ind1], tmpmat[, ind2])
    }
    return(list(ind1 = ind1, ind2 = ind2, cp = cp))
  }, cl=cl)
  parallel::stopCluster(cl)
  gc();gc()

  ### build the output matrix block by block
  if(verbose>0|verbose) cat(" - bigcrossprod: merging the batches into", p, "x", p, "output matrix\n")
  outmat <- matrix(0, ncol = p, nrow = p)
  for(i in 1:length(blocks)){
    outmat[blocks[[i]]$ind1, blocks[[i]]$ind2] <- blocks[[i]]$cp
  }
  outmat <- Matrix::forceSymmetric(outmat, uplo = "U")
  outmat <- as.matrix(outmat)
  rownames(outmat) <- colnames(mat)
  colnames(outmat) <- colnames(mat)

  ### convert cov. matrix to correlation matrix
  outmat <- cov2cor(outmat)

  ### output
  file.remove(file.path(tmp.dir, mat.binfile))
  gc();gc()
  return(outmat)
}
