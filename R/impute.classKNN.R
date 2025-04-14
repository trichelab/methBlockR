#' impute missing sites via k-NN with user-defined classes for comparables
#'
#' @param x       a MethylationExperiment (or something with getBeta)
#' @param col     a factor with the class of each column, or its name in colData
#' @param NAs     maximum row NA fraction (0 to 1) allowed within a class (0.5)
#' @param maxK    maximum K (positive integer) nearest neighbors to use (99)
#' @param BPPARAM a BiocParallel::MulticoreParam() or similar for processing
#' @param ...     other arguments to pass along to impute::impute.knn
#'
#' @details       in principle, this could be combined with impute.GimmeCpG
#'                in practice, can iterate this with coarsening 'col' factors
#'
#' @return        x, but with a new assay named ImputedBeta 
#' 
#' @seealso       impute::impute.knn
#'
#' @importFrom    impute impute.knn
#' 
#' @export
impute.classKNN <- function(x, col, NAs=.5, maxK=99, BPPARAM=NULL, ...){

  stopifnot(is(x, "SummarizedExperiment") & ("Beta" %in% assayNames(x)))
  if (length(col) != ncol(x)) {
    if (col %in% names(colData(x))) col <- colData(x)[, col]
    else stop("col does not match ncol(x)")
  }

  # Delayed-unfriendly
  B0 <- lapply(split(t(getBeta(x)), col), t)
  if (is.null(BPPARAM)) BPPARAM <- SerialParam(progressbar=TRUE)
  B <- bplapply(B0, .imputeWithin, NAs=NAs, maxK=maxK, BPPARAM=BPPARAM, ...)
  assay(x, withDimnames=FALSE, "ImputedBeta") <- do.call(cbind, B)[,colnames(x)]
  return(x)

}


# helper fn
.imputeWithin <- function(x, NAs=0.5, maxK=99, ...) {

  k <- pmin(maxK, pmax(1, round(ncol(x)/2)))
  rr <- which(missingness(x, 1, "ImputedBeta") <= maxNA)
  x[rr, ] <- fexpit(suppressMessages(impute.knn(flogit(x[rr, ]))$data))
  return(x)

}

