#' impute missing sites via k-NN with user-defined classes for comparables
#'
#' @param X       a MethylationExperiment (or something with getBeta)
#' @param INDEX   a factor representing the class of each observation
#' @param ...     other arguments to pass along to impute::impute.knn
#'
#' @details       in principle, this can be combined with impute.GimmeCpG
#'
#' @return        x, but with a new assay Imputed that has been imputed.
#' 
#' @seealso       impute::impute.knn
#'
#' @importFrom    impute impute.knn
#' 
#' @export
impute.classKNN <- function(X, INDEX, ...) { 

  if (!is(X, "SummarizedExperiment") | !("Beta" %in% assayNames(X))) { 
    stop("Could not find beta values to impute")
  }

  tbl <- table(INDEX) 
  k <- max(1, as.integer(tbl) - 1)
  names(k) <- names(tbl) 
  assay(X, "Imputed", withDim=FALSE) <- .imputeByClass(X, INDEX, k)
  return(X) 

}


# helper fn
.imputeByClass <- function(X, INDEX, k, ...) { 

  B <- getBeta(X)
  for (i in names(k)) { 
    nn <- k[i]
    j <- which(INDEX == i)
    message("Imputing missing values for class ", i, "...") 
    B[, j] <- fexpit(impute.knn(flogit(B[, j]), k=nn, ...)$data)
  }
  return(B)

}
