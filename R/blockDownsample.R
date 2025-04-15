#' block-randomized sampling to ensure somewhat representative comparisons 
#'
#' @param x       a MethylationExperiment or comparable object
#' @param col     a factor for each column, or its name in colData
#' @param SNPcall recall SNPs if present? (TRUE)
#'
#' @details       use set.seed beforehand if you want reproducibility.
#'
#' @return        a subset of x, downsampled to match the smallest N(col)
#' 
#' @export
blockDownsample <- function(x, col, SNPcall=TRUE) {

  if (length(col) != ncol(x)) {
    if (col %in% names(colData(x))) col <- colData(x)[, col]
    else stop("col does not match ncol(x)")
  }

  ss <- levels(factor(col))
  names(ss) <- ss
  n <- min(table(col))
  idx <- sort(do.call(c, lapply(ss, function(s) sample(which(col == s), n))))
  res <- x[, idx]
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

  return(res)

}
