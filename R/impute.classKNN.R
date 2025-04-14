#' impute missing sites via k-NN with user-defined classes for comparables
#'
#' @param x       a MethylationExperiment (or something with getBeta)
#' @param col     a factor with the class of each column, or its name in colData
#' @param ...     other arguments to pass along to impute::impute.knn
#'
#' @details       in principle, this can be combined with impute.GimmeCpG
#'
#' @return        x, but with a new assay named ImputedBeta 
#' 
#' @seealso       impute::impute.knn
#'
#' @importFrom    impute impute.knn
#' 
#' @export
impute.classKNN <- function(x, col, ...) { 

  if (!is(x, "SummarizedExperiment") | !("Beta" %in% assayNames(x))) { 
    stop("Could not find beta values to impute")
  }
  if (length(col) != ncol(x)) {
    if (col %in% names(colData(x))) { 
      col <- colData(x)[, col]
    } else {
      stop("col does not match ncol(x)")
    }
  }

  tbl <- table(col)
  assay(x, "ImputedBeta") <- getBeta(x)
  k <- pmax(1, (as.integer(tbl) - 1))
  names(k) <- names(tbl) 

  # replace with bplapply?
  for (i in names(k)) {
    nn <- k[i]
    j <- which(col == i)
    message("Imputing missing values for class ", i, "...") 
    assay(x[, j], "ImputedBeta") <- 
      fexpit(impute.knn(flogit(assay(x[, j], "ImputedBeta")), k=nn, ...)$data)
  }

  return(x)

}
