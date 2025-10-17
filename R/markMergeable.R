#' find mergeable methylation blocks by correlation (raw or by NMF weights)
#' 
#' @param x     a MethBlockExperiment or similar, perhaps "nmf_fit" in metadata
#' @param how   correlation ('corr') or NMF weight correlation ('nmf')
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
markMergeable <- function(x, how=c("corr","nmf"), thr=.9) {

  x <- sort(x) 
  how <- match.arg(how) 

  dst <- max(width(x))   # in bp
  mcols(x)$nearest <- mcols(distanceToNearest(x))$distance
  mcols(x)$eligible <- as.integer(rowRanges(x)$nearest <= dst)
 
  # need to split by chrom first
  byChr <- lapply(.splitByChrom(x), .nameBlocks)

  byBlock <- split(x, mcols(x)$block)
  message("Don't forget to skip byBlock[[1]], which is actually byBlock$0!")

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


# helper fn 
.splitByChrom <- function(x) split(x, seqnames(x))


# helper fn; assumes by chr
.nameBlocks <- function(x) { 

  chr <- unique(seqnames(x))
  stopifnot(length(chr) == 1)

  blocks <- rle(mcols(x)$eligible)
  nonzero <- which(blocks$values > 0)
  blocks$values[nonzero] <- paste0(chr, "_block", seq_len(length(nonzero)))
  mcols(x)$block <- inverse.rle(blocks)

  return(x) 

}
