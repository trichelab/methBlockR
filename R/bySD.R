#' bySD
#'
#' This function selects the k most variable rows of x and returns their values.
#'
#' @param x       a numeric matrix, possibly with NAs 
#' @param k       how many rows to return (500)
#' @param minmean minimum mean value (-Inf)
#' @param maxmean maximum mean value (+Inf) 
#'
#' @return        the most variable _k_ rows of _x_ (with suitable means)
#' 
#' @importFrom    matrixStats rowMeans2 rowSds
#' 
#' @export
bySD <- function(x, k=500, minmean=-Inf, maxmean=Inf) { 
  means <- rowMeans2(x, na.rm=TRUE)
  eligible <- which(means < maxmean & means > minmean)
  sds <- rowSds(x[eligible, ], na.rm=TRUE)
  k <- min(length(eligible), k)
  x[rev(order(sds))[seq_len(k)], ]
}
