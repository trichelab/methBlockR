#' bySD
#'
#' This function selects the k most variable rows of x and returns their values.
#'
#' @param x       a numeric matrix, possibly with NAs 
#' @param k       maximum number of rows to return (100)
#' @param minmean minimum mean value (0.025)
#' @param maxmean maximum mean value (0.975)
#'
#' @return        the most variable rows of x (with suitable means)
#' 
#' @importFrom    matrixStats rowMeans2 rowSds
#' 
#' @export
#' 
bySD <- function(x, k=100, minmean=0.025, maxmean=0.975) { 

  x <- filterByMean(x, min=minmean, max=maxmean)
  k <- min(nrow(x), k)
  sds <- rowSds(x, na.rm=TRUE)
  x[rev(order(sds))[seq_len(k)], ]

}
