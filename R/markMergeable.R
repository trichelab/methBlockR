#' find mergeable methylation blocks by correlation (raw or by NMF weights)
#' 
#' @param x     a MethBlockExperiment, perhaps with "nmf_fit" in metadata
#' @param how   correlation ('r') or NMF weight correlation ('nmf')
#' @param thr   threshold for merging neighboring blocks (0.9; see Details)
#'
#' @return      a list, GRanges, or MBE of merged blocks based on (how, thr)
#'
#' @details     blocks more than max(width(blocks)) bp apart are not checked
#'
#' @import      BiocParallel
#' @import      RcppML
#' @import      RSC
#'
#' @export
#'
markMergeable <- function(x, how=c("r","nmf"), thr=.9) {

  how <- match.arg(how) 
  
  dst <- max(width(x))   # in bp
  mcols(x)$nearest <- mcols(distanceToNearest(x))$distance
  mcols(x)$eligible <- as.integer(rowRanges(x)$nearest <= dst)
  blocks <- rle(mcols(x)$eligible)
  nonzero <- which(blocks$values > 0)
  blocks$values[nonzero] <- seq_len(length(nonzero))
  mcols(x)$block <- inverse.rle(blocks)

  byBlock <- split(x, mcols(x)$block)
  warn("Don't forget to skip byBlock[[1]], which is actually byBlock$0!")

  # key: compute banded correlations between N adjacent blocks, 
  #      then do a greedy merge along Rles of cor > thr 
  #
  if (how == "nmf") { 
    stopifnot("nmf_fit" %in% names(metadata(x)))
    nmf <- metadata(x)$nmf_fit
  } else { 

  } 

  stop("Not finished") 

}


# helper fn
.splitByMcol <- function(x, mcol) split(x, mcols(x)[[mcol]])
