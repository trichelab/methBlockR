#' Collapse a GenomicRatioSet or similar over methylation blocks (per Kaplan)
#' 
#' By default, data(methBlocks) is used for summarization, singleton probes are
#' included, methylation block rownames are returned against hg19, and neither 
#' singletons (!proper) nor unstable mappings (!stable) are filtered out. This
#' version of asMethBlocks supports "custom" as a target; be careful with that.
#'
#' @param x         a SummarizedExperiment-derived object with CpG loci as rows
#' @param g         genome to use for methylation block coordinates (hg19)
#' @param custom    if genome == "custom", a GRanges (cannot be NULL) of blocks
#'
#' @return          an object with same colData but with new rowRanges & assays
#' 
#' @details
#' If you have specified a genome() for x, which you very much should, this 
#' function will raise an error if the genomes or seqlevels do not match. A copy
#' of the subset of methBlocks used for mapping will be returned in metadata(x).
#' 
#' @seealso         switchMethBlocksGenome
#' 
#' @import          BiocParallel 
#' @import          GenomicRanges 
#'
#' @export
asMethBlocks <- function(x, g=c("hg19","hg38","mm10","custom"), custom=NULL) {

  g <- match.arg(g)
  stopifnot(is(x, "SummarizedExperiment"))
  if (g == "custom" & (is.null(custom) | !is(custom, "GRanges"))) {
    stop("You must provide a GRanges of blocks if using custom blocks.")
  } else {
    # FIXME: this is gross 
    if (exists("methBlocks")) {
      message("You seem to have defined `methBlocks` already. Using that.")
    } else { 
      data("methBlocks", package="methBlockR") 
    }
  }

  message("Checking for masked (NA) probes...")
  N <- ncol(x)
  keepCpGs <- rownames(x)[rowSums(is.na(getBeta(x))) < N]
  M <- length(keepCpGs)
  dropCpGs <- setdiff(rownames(x), keepCpGs)
  M_NA <- length(dropCpGs)
  message("Retained ", M, " CpGs, dropped ", M_NA, " with ", N, "/", N, " NAs.")

  message("Mapping probes to methylation blocks in ", g, " genome...")
  methBlocks <- methBlocks[keepCpGs, ]
  methBlocks <- subset(methBlocks, !is.na(methBlocks[[g]]))
  keepCpGs <- intersect(keepCpGs, rownames(methBlocks))
  M <- length(keepCpGs)
  message("Retained ", M, " probes mapped to ", g, ".")

  # preparation for aggregation
  tbl <- table(methBlocks[[g]])
  blocks <- names(tbl)[tbl > 1] # blocks with more than one probe
  byBlock <- subset(methBlocks, methBlocks[[g]] %in% blocks)[, g, drop=FALSE]

  # this will become rowRanges(amb)
  mbgr <- as(methBlocks[[g]], "GRanges")
  genome(mbgr) <- g 
  mbgr$probes <- countOverlaps(mbgr, granges(x)[keepCpGs])
  mbgr$singleton <- width(mbgr) == 1
  mbgr <- unique(mbgr)
  names(mbgr) <- as.character(mbgr)

  message("Computing per-block average methylation across ", M, " probes...")
  placeholder <- rownames(subset(methBlocks, !duplicated(methBlocks[[g]])))
  blockBetas <- getBeta(x[placeholder, ])
  rownames(blockBetas) <- names(mbgr) 

  # only recompute regions with more than one probe
  byBlock[[g]] <- factor(byBlock[[g]])
  toSquash <- split.data.frame(getBeta(x[rownames(byBlock), ]), byBlock[[g]])
  res <- do.call(rbind, lapply(toSquash, colMeans, na.rm=TRUE)) # bplapply fails
  blockBetas[rownames(res), ] <- res
 
  if (!is(x, "MethylationExperiment")) message("Reconstructing ",class(x),"...")
  amb <- x[seq_len(nrow(blockBetas)), ]
  rownames(amb) <- rownames(blockBetas)
  assay(amb, "Beta") <- blockBetas
  rowRanges(amb) <- mbgr

  message("Done.")
  return(amb)

}
