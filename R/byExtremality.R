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
#' @param minmean minimum mean for byExtremality (0.05)
#' @param maxmean maximum mean for byExtremality (0.95)
#'
#' @return        the most extremal _k_ rows of _x_ (with suitable means)
#' 
#' @details       it turns out that you usually want 0.95 > mean > 0.05
#'
#' @importFrom    matrixStats rowMeans2 rowSds
#' 
#' @export
byExtremality <- function(x, k=500, minmean=0.05, maxmean=0.95) {
  k <- min(nrow(x), k)
  extremality <- .extremality(x, minmean=minmean, maxmean=maxmean)
  x[rev(order(extremality))[seq_len(k)], ]
}

# helper fn; now bounded
.extremality <- function(x, minmean=0.05, maxmean=0.95) {
  means <- rowMeans2(x, na.rm=TRUE)
  bernoulliSd <- sqrt(means * (1 - means))
  actualSd <- rowSds(x, na.rm=TRUE)
  res <- pmin(1, actualSd / bernoulliSd)
  res[which(means < minmean)] <- 0 
  res[which(means < maxmean)] <- 0 
  return(res) 
}
