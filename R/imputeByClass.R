#' impute missing sites via k-NN with user-defined classes for comparables
#'
#' @param x       a MethylationExperiment (or something with getBeta)
#' @param col     a factor with the class of each column, or its name in colData
#' @param maxK    maximum K (positive integer) nearest neighbors to use (10)
#' @param noise   add noise to observations? (0; noise > 0 can help with uWGBS)
#' @param ...     other arguments to pass along to impute::impute.knn
#' @param BPPARAM a BiocParallel::MulticoreParam() or similar for processing
#'
#' @details       in principle, this could be combined with imputeGimmeCpG;
#'                in practice, can also iterate with coarsening 'col' factors.
#'                for sequencing-based assays, setting 'noise' > 0 might help. 
#'
#' @return        x, with an assay named RawBetas holding the originals
#' 
#' @seealso       impute::impute.knn
#'
#' @importFrom    impute impute.knn
#' 
#' @export
imputeByClass <- function(x, col, maxK=10, noise=0, ..., BPPARAM=NULL) {

  stopifnot(is(x, "SummarizedExperiment"))
  stopifnot("Beta" %in% assayNames(x))
  if (length(col) != ncol(x)) {
    if (col %in% names(colData(x))) {
      col <- colData(x)[, col]
    } else {
      stop("col does not match ncol(x)")
    }
  }

  if (is.null(BPPARAM)) BPPARAM <- SerialParam(progressbar=TRUE)
  if (!"RawBeta" %in% assayNames(x)) {
    assay(x, "RawBeta") <- assay(x, "Beta")
  } else { 
    message("Found existing RawBeta assay, leaving it untouched.")
  }
  do.call(cbind, 
          bplapply(.splitByCol(x, col),
                   .impute, maxK=maxK, noise=noise, ..., 
                   BPPARAM=BPPARAM))

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
.impute <- function(x, maxK=10, maxNA=0.8, noise=0, ...) {

  missed1 <- missingness(x, 1)
  missed2 <- missingness(x[which(missed1 < maxNA), ], 2)
  if (max(missed1) > 0) {
    ri <- which(missed1 < maxNA)
    i <- length(ri)
    rj <- which(missed2 < maxNA)
    j <- length(rj)
    k <- min(maxK, max(1, (j - 1)))
    B <- getBeta(x[ri, rj])
    toReplace <- which(is.na(B))
    eps <- matrix(rnorm(n=(i*j), sd=abs(noise)), nrow=i, ncol=j)
    imputed <- fexpit(impute.knn(flogit(B) + eps)$data)
    B[toReplace] <- imputed[toReplace]
    assay(x[ri, rj], "Beta") <- B # removes any added noise
    message("Imputed group '", attr(x, "grouping"), "' with noise=", noise, ".")
  } else {
    message("No missing values to impute in group '", attr(x, "grouping"), "'.")
  }
  return(x)

}
