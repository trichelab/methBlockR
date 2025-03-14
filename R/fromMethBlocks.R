#' extract methylation rates at probes corresponding to methylation blocks and
#' expand the matrix to match the original input dimensions, for (e.g.) data 
#' transfer or comparison across platforms and technologies (arrays, ONT, etc.)
#' 
#' see documentation for asMethBlocks and toMethBlocks for more details 
#' 
#' @param x         a GenomicRatioSet or something like it with chr12:34-56 rows
#' @param loci      probe coordinates (or a GRanges) to unpack blocks onto 
#'
#' @return          an object reinflated to locus level using rates
#' 
#' @seealso         switchMethBlocksGenome
#' @seealso         toMethBlocks
#' @seealso         asMethBlocks
#' 
#' @import          SummarizedExperiment
#' @import          GenomicRanges
#'
#' @export
fromMethBlocks <- function(x, loci) {

  g0 <- unique(genome(x))
  if (is.na(g0)) {
    warning("No genome specified for `x`. You're on thin ice here...")
    g0 <- "hg19"
  }
  if (is(loci, "GRanges")) {
    g1 <- unique(genome(loci))
    if (is.na(g1)) {
      warning("No genome specified for `loci`. Make sure it matches!")
      genome(loci) <- g0
    } else {
      stopifnot(g1 == g0)
    }
  } else {
    loci <- as(loci, "GRanges")
    genome(loci) <- g0
  }

  # find the loci that overlap the blocks in the dataset 
  message("Mapping methylation blocks to probes in ", g, "...")
  newBetas <- matrix(NA_real_, nrow=length(loci), ncol=ncol(x))
  rownames(newBetas) <- names(loci)
  colnames(newBetas) <- colnames(x)
  ol <- findOverlaps(loci, granges(x))
  newBetas[queryHits(ol), ] <- getBeta(x)[subjectHits(ol), ]
  SummarizedExperiment(assays=SimpleList(Beta=newBetas),
                       rowRanges=loci,
                       colData=colData(x),
                       metadata=metadata(x))

}
