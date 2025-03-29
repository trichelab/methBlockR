#' byExtremality 
#'
#' For array-based DNA methylation, particularly rates across regions,
#' we can do better (a lot better) than MAD. An ideal class-splitting Beta-
#' distributed feature will have Bernoulli variance (== all 0s and 1s): 
#' 
#' max(SD(X_j)) if X_j ~ Beta(a, b) <= max(SD(X_j)) if X_j ~ Bernoulli(a/(a+b))
#'
#' for X having a known mean and SD, solvable for a + b by MoM. Now define
#'
#' extremality(X_j) = sd(X_j) / bernoulliSD(mean(X_j))
#'
#' This function selects the k most extremal rows of x and returns their values.
#'
#' @param     x   a matrix of beta values (proportions), by block or by locus
#' @param     k   how many rows to return (500)
#' @return    the most extremal _k_ rows of _x_
#' 
#' @import    matrixStats
#' 
#' @export
byExtremality <- function(x, k=500) {
  k <- min(nrow(x), k)
  extremality <- .extremality(x)
  x[rev(order(extremality))[seq_len(k)], ]
}

# helper fn
.extremality <- function(x) {
  means <- rowMeans2(x, na.rm=TRUE)
  bernoulliSd <- sqrt(means * (1 - means))
  actualSd <- rowSds(x, na.rm=TRUE)
  pmin(1, actualSd / bernoulliSd)
}
