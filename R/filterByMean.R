#' filterByMean
#' 
#' does what it says
#' 
#' @param x     a matrix to filter
#' @param min   minimum allowable mean (0.025)
#' @param max   maximum allowable mean (0.975)
#' @param INDEX which dimension to filter along (1)
#' 
#' @return      only the permissible (rows/columns) of X
#'
#' @details     this could probably be adapted to operate on INDEX=1:2
#'
#' @importFrom  matrixStats rowMeans2 colMeans2
#'
#' @export
#'
filterByMean <- function(x, min=0.025, max=0.975, INDEX=1) {

  stopifnot(INDEX %in% 1:2)
  means <- switch(INDEX, rowMeans2(x, na.rm=TRUE), colMeans2(x, na.rm=TRUE))
  eligible <- which(means <= maxmean & means >= minmean)
  rows <- ifelse(INDEX == 1, eligible, seq_len(nrow(x)))
  cols <- ifelse(INDEX == 2, eligible, seq_len(ncol(x)))
  x[rows, cols]

}
