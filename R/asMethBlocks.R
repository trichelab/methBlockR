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

  if (is(x, "GenomicRatioSet") & any(width(rowRanges(x))>50)) { 
    message("Converting a blocked GenomicRatioSet to a MethBlockExperiment...")
    # just convert to MBE 
    g <- unique(genome(x))
    if (!g %in% names(mcols(x))) {
      data("methBlocks", package="methBlockR")
      genomes <- c("hg19", "hg38") 
      gg <- setdiff(genomes, g) 
      mcols(x)[, g] <- methBlocks[match(rownames(fsb), methBlocks[[g]]), g]
      mcols(x)[, gg] <- methBlocks[match(rownames(fsb), methBlocks[[g]]), gg]
    }
    return(as(as(x, "SingleCellExperiment"), "MethBlockExperiment"))
  }

  g <- match.arg(g)
  if (is(x, "MethylationExperiment")) { 
    # coerce the object so that we can pass it through the constructor
    assay(x, "CN") <- NULL # drop intensities since we move to rates
    return(MethBlockExperiment(as(x, "SingleCellExperiment"), g=g))
  }

  stopifnot(is(x, "SummarizedExperiment"))
  if (g == "custom" & (is.null(custom) | !is(custom, "GRanges"))) {
    stop("You must provide a GRanges of blocks if using custom blocks.")
  } else {
    data("methBlocks")
  }

  message("Checking for masked (NA) probes...")
  N <- ncol(x)
  keepCpGs <- rownames(x)[rowSums(is.na(assay(x, "Beta"))) < N]
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
  mbgr <- unique(as(methBlocks[[g]], "GRanges"))
  names(mbgr) <- as.character(mbgr)
  genome(mbgr) <- g 

  # rowData for the results 
  mbgr$probes <- tbl[names(mbgr)]
  mbgr$singleton <- width(mbgr) == 1
  data(mbHg19Hg38) # need to expand this to mm10, T2T, HPRC, etc. 
  ord <- match(names(mbgr), mbHg19Hg38[[g]])
  mcols(mbgr)[[g]] <- mbHg19Hg38[ord, g]
  for (h in setdiff(names(mbHg19Hg38), g)) mcols(mbgr)[[h]] <- mbHg19Hg38[ord,h]

  message("Computing per-block average methylation across ", M, " probes...")
  placeholder <- rownames(subset(methBlocks, !duplicated(methBlocks[[g]])))
  blockBetas <- assay(x[placeholder, ], "Beta")
  rownames(blockBetas) <- names(mbgr) 

  # only recompute regions with more than one probe
  byBlock[[g]] <- factor(byBlock[[g]])
  toSquash <- split.data.frame(assay(x[rownames(byBlock), ], "Beta"), 
                               byBlock[[g]])
  res <- do.call(rbind, lapply(toSquash, colMeans, na.rm=TRUE)) # bplapply fails
  blockBetas[rownames(res), ] <- res
 
  amb <- x[seq_len(nrow(blockBetas)), ]
  rownames(amb) <- rownames(blockBetas)
  assay(amb, "Beta") <- blockBetas
  rowRanges(amb) <- mbgr
  message("Done.")
  return(amb)

}
