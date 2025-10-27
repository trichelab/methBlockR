#' find mergeable methylation blocks by correlation (raw or by NMF weights)
#' 
#' @param x     a MethBlockExperiment or similar, perhaps "nmf_fit" in metadata
#' @param how   correlation ('corr') or NMF weight correlation ('nmf')
#' @param dst   how many base pairs (max) to look in each direction? (1000)
#' @param p     threshold p(cor) for merging neighboring blocks (0.01; Details)
#' @param minr  minimum acceptable correlation (regardless of p) to merge (0.5)
#'
#' @details     Given n pairs of observations, a statistical test of correlation
#'              can be performed against a t distribution with n-2 df, where 
#'              the test statistic is t* = (r * sqrt(n-2)) / sqrt(1 - (r**2)),
#'              and the one-side p-value is pt(t*, df=n-2, lower.tail=FALSE). 
#'              Since we don't care about strong negative correlation, we only
#'              test for positive correlations unlikely to be due to chance,
#'              so we find t* = qt(1 - (p/dstblocks(x)), df=ncol(x) - 2),
#'              where dstblocks(x) is the number of blocks in metablocks,
#'              then solve for the minimum correlation r generating t* with n-2,
#'              and that r is our lower bound to merge blocks within metablocks.
#'
#' @return      an Rle of mergeable block candidates based on (how, dst, thr)
#'
#' @import      BiocParallel
#' @import      RcppML
#'
#' @export
#'
markMergeable <- function(x, how=c("corr","nmf"), dst=1000, p=0.01, minr=0.5) {

  x <- sort(sortSeqlevels(x))
  mcols(x)$nearest <- mcols(distanceToNearest(x))$distance
  mcols(x)$eligible <- as.integer(rowRanges(x)$nearest <= dst)
 
  # need to split by chrom first
  x <- do.call(rbind, lapply(.splitByChrom(x), .nameMetaBlocks))
  mcols(x)$metablock <- Rle(factor(mcols(x)$metablock))
  noblock <- which(levels(mcols(x)$metablock) == "0")
  levels(mcols(x)$metablock)[noblock] <- NA_character_

  message("Don't forget to skip 'metablock' 1, which is actually 0 (noblock)")

  # significance based thresholding (only overrides minr if greater than minr) 
  r_test <- seq(0.01, 0.99, 0.01)
  p_star <- p / sum(!is.na(mcols(x)$metablock))
  t_star <- qt(1 - p_star, df = ncol(x) - 2)
  r_stat <- .oneSidedCor(r=r_test, n=ncol(x))
  r_thr <- r_test[which(r_stat > t_star)[1]]
  if (r_thr > minr) minr <- r_thr
  message("Threshold for merging set at ", minr)

  # key: compute correlations between adjacent blocks, 
  #      then do a greedy merge along Rles if cor > thr 
  #
  how <- match.arg(how) 
  if (how == "nmf") { 
    topifnot("nmf_fit" %in% names(metadata(x)))
    message("this may be wrong") 
    wts <- metadata(x)$nmf_fit
  } else {
         
  }

  

  stop("Not finished") 

}


# helper fn
.splitByMcol <- function(x, mcol) split(x, mcols(x)[[mcol]])


# helper fn 
.splitByChrom <- function(x) split(x, seqnames(x))


# helper fn; assumes by chr
.nameMetaBlocks <- function(x) { 

  chr <- unique(seqnames(x))
  stopifnot(length(chr) == 1)
  stopifnot("eligible" %in% names(mcols(x)))
  message("Labeling metablocks on ", chr, "... ", appendLF=FALSE) 
  blox <- rle(mcols(x)$eligible)
  nonzero <- which(blox$values > 0)
  blox$values[nonzero] <- paste0(chr, "_metablock", seq_len(length(nonzero)))
  mcols(x)$metablock <- inverse.rle(blox)
  message("done.")
  return(x) 

}


# helper fn
.oneSidedCor <- function(r, n) (r * sqrt(n-2)) / sqrt(1 - (r**2))
