#' block-randomized sampling to ensure somewhat representative comparisons 
#'
#' @param x       a MethylationExperiment or comparable rectangular object
#' @param col     a factor for each row/column, or its name in colData
#' @param maxN    maximum size per-subset? (Inf; use the smallest subset size)
#' @param SNPcall recall SNPs if present? (FALSE)
#' @param INDEX   block-downsample on columns (2) or rows (1)? (2) 
#'
#' @details       use set.seed beforehand if you want reproducibility.
#'
#' @return        a subset of x, downsampled to match the smallest N(col)
#' 
#' @export
blockDownsample <- function(x, col, maxN=Inf, SNPcall=FALSE, INDEX=2) {

  if (length(col) != dim(x)[INDEX]) { 
    if (is(x, "SummarizedExperiment") & INDEX == 2) {
      if (col %in% names(colData(x))) {
        col <- colData(x)[, col]
      }
    } else {
      stop("col does not match dim(x)[", INDEX, "]")
    }
  }

  ss <- levels(factor(col))
  names(ss) <- ss
  n <- min(maxN, min(table(col)))
  idx <- sort(do.call(c, lapply(ss, function(s) sample(which(col == s), n))))
  res <- switch(INDEX, x[idx, ], x[, idx])

  if (is(x, "SummarizedExperiment") & INDEX == 2) {
    mdn <- names(metadata(res))
    tbl <- table(mdn) 
    if ("SNPs" %in% names(tbl)) { 
      if (all(colnames(res) %in% colnames(metadata(res)$SNPs) == ncol(x))) { 
        metadata(res)[["SNPs"]] <- metadata(res)[["SNPs"]][, colnames(res)]
        if ("SNPcalls" %in% names(tbl)) {
          metadata(res)[which(mdn == "SNPcalls")] <- NULL
          if (SNPcall) {
            message("Recalling SNPs for the merged dataset.") 
            metadata(res)[["SNPcalls"]] <- SNPcalls(res)
          }
        }
      } else {  
        message("Missing columns in metadata(<",class(res),">)$SNPs; skipped")
      }
    }
  }

  return(res)

}
