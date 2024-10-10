#' assign loci to regions (blocks) using standard GenomicRanges tools
#' 
#' This is a somewhat trivial function that I got tired of writing over & over
#'
#' @param loci      a GRanges of (usually smaller) regions
#' @param blocks    a GRanges of (usually larger) regions 
#' @param fill      assign unmatched loci to their own blocks? (FALSE)
#' 
#' @return          the loci but with mcols(loci)$block filled with block names
#' 
#' @details
#' If you have specified a genome() for loci and blocks, which you ought to, 
#' this function will raise an error if the genomes or seqlevels do not match.
#' 
#' @seealso         asMethBlocks
#' 
#' @import          GenomicRanges 
#'
#' @export
#'
assignBlocks <- function(loci, blocks, fill=FALSE) { 

  ol <- findOverlaps(loci, blocks)
  loci$block <- NA
  if (length(ol) == 0) {
    warning("No overlaps found between loci and blocks!")
    return(loci)
  }

  unstrand <- function(x) { strand(x) <- "*"; return(x) }
  if (fill) loci$block <- as.character(unstrand(loci))
  loci[queryHits(ol)]$block <- as.character(blocks[subjectHits(ol)])
  loci$block <- factor(loci$block) 
  return(loci)

}
