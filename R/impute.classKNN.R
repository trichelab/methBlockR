#' impute missing sites via k-NN with user-defined classes for comparables
#'
#' @param x       a MethylationExperiment (or something with getBeta)
#' @param col     a factor with the class of each column, or its name in colData
#' @param maxK    maximum K (positive integer) nearest neighbors to use (10)
#' @param BPPARAM a BiocParallel::MulticoreParam() or similar for processing
#' @param ...     other arguments to pass along to impute::impute.knn
#'
#' @details       in principle, this could be combined with impute.GimmeCpG
#'                in practice, can iterate this with coarsening 'col' factors
#'
#' @return        x, but with an assay named RawBetas holding the originals
#' 
#' @seealso       impute::impute.knn
#'
#' @importFrom    impute impute.knn
#' 
#' @export
impute.classKNN <- function(x, col, maxK=10, BPPARAM=NULL, ...){

  stopifnot(is(x, "SummarizedExperiment") & ("Beta" %in% assayNames(x)))
  if (length(col) != ncol(x)) {
    if (col %in% names(colData(x))) col <- colData(x)[, col]
    else stop("col does not match ncol(x)")
  }

  if (is.null(BPPARAM)) BPPARAM <- SerialParam(progressbar=TRUE)
  assay(x, withDimnames=FALSE, "RawBeta") <- assay(x, "Beta")
  do.call(cbind, bplapply(.splitByCol(x, col),
                          .impute, maxK=maxK, ..., BPPARAM=BPPARAM))

}


# helper fn
.splitByCol <- function(x, col) { 
 
  ss <- levels(factor(col))
  names(ss) <- ss 
  res <- lapply(ss, function(s) x[, which(col == s), drop=FALSE])
  # has to be done immediately before handoff to .impute()
  for (s in names(res)) attr(res[[s]], "grouping") <- s
  names(res) <- ss
  return(res)

}


# helper fn
.impute <- function(x, maxK=10, maxNA=0.8, ...) {

  missed1 <- missingness(x, 1)
  missed2 <- missingness(x[which(missed1 < maxNA), ], 2)
  if (max(missed1) > 0) {
    ri <- which(missed1 < maxNA)
    rj <- which(missed2 < maxNA)
    k <- min(maxK, max(1, length(rj) - 1))
    message("Imputing group '", attr(x, "grouping"), "'...")
    Brr <- fexpit(suppressMessages(impute.knn(flogit(getBeta(x[ri, rj])))$data))
    assay(x[ri, rj], "Beta") <- Brr
  }
  message("Finished imputing group '", attr(x, "grouping"), "'.")
  return(x)

}
