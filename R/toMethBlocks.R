#' extract the probes corresponding to methylation blocks and collapse them
#' 
#' see documentation for asMethBlocks and fromMethBlocks for more details 
#' 
#' @param x         a GenomicRatioSet or something like it with cgXXXXXXXXX rows
#' @param y         block coordinates (or a GRanges) to summarize x over 
#' @param g         genome to use for block coordinates if none found (hg19)
#'
#' @return          an object or matrix summarized down to the methBlocks in y
#' 
#' @seealso         switchMethBlocksGenome
#' @seealso         fromMethBlocks
#' @seealso         asMethBlocks
#' 
#' @import          BiocParallel 
#'
#' @export
toMethBlocks <- function(x, y, g=c("hg19","hg38"), methBlocks=NULL) {

  g <- match.arg(g)
  if (is(y, "GRanges")) {
    g1 <- unique(genome(y))
    if (!is.na(g1) & g1 %in% c("hg19", "hg38")) g <- g1
    else genome(y) <- g
  } else {
    y <- as(y, "GRanges")
    genome(y) <- g
  }

  if (!is.null("methBlocks")) {
    message("You seem to have defined `methBlocks` already. Using that.")
  } else { 
    data("methBlocks", package="miser") 
  }

  # check for 1:1 overlap:
  if (all(as(y, "character") %in% methBlocks[[g]])) { 
    CpGs <- rownames(subset(methBlocks, methBlocks[[g]] %in% as(y,"character")))
    CpGs <- intersect(CpGs, rownames(x))
  } else { 
    stop("Currently, only 1:1 methBlock:BED naming is supported.")
  }
  x <- x[CpGs, ] 

  # check for masked probes 
  message("Checking for masked (NA) probes...")
  N <- ncol(x)
  keepCpGs <- rownames(x)[rowSums(is.na(x)) < N]
  M <- length(keepCpGs)
  dropCpGs <- setdiff(rownames(x), keepCpGs)
  M_NA <- length(dropCpGs)
  x <- x[keepCpGs, ] 
  M <- length(keepCpGs)
  CpGs <- intersect(CpGs, keepCpGs)
  message("Retained ", M, " CpGs, dropped ", M_NA, " with ", N, "/", N, " NAs.")

  message("Mapping probes to methylation blocks in ", g, "...")
  mb <- methBlocks[CpGs, g, drop=FALSE]
  mb <- subset(mb, !is.na(mb[[g]]))
  mb[[g]] <- factor(mb[[g]]) 

  # preparation for aggregation
  tbl <- table(mb[[g]])
  blocks <- names(tbl)[tbl > 1] # blocks with more than one probe
  byBlock <- subset(mb, mb[[g]] %in% blocks)[, g, drop=FALSE]

  message("Computing per-block average methylation across ", M, " probes...")
  placeholder <- rownames(subset(mb, !duplicated(mb[[g]])))
  blockBetas <- x[placeholder, ]
  rownames(blockBetas) <- unique(mb[[g]])
  byBlock[[g]] <- factor(byBlock[[g]])
  toSquash <- split.data.frame(x[rownames(byBlock), ], byBlock[[g]])
  res <- do.call(rbind, lapply(toSquash, colMeans, na.rm=TRUE)) 
  blockBetas[rownames(res), ] <- res
  message("Done.")
  return(blockBetas)

}
