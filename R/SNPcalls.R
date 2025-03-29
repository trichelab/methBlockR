#' fit a 1-to-3-component mixture model to flogit(SNPbetas) for each sample
#' 
#' @param SNPs      a SummarizedExperiment or a matrix of SNP beta values 
#'
#' @return          0/1/2 calls for each SNP in each sample
#' 
#' @details         mixR::mixfit(flogit(x), ncomp=1:3) is run for each sample. 
#' 
#' @import          mixR
#' 
#' @export
SNPcalls <- function(SNPs) {

  if (is(SNPs, "SummarizedExperiment")) {
    NC <- ncol(SNPs)
    if ("SNPs" %in% names(metadata(SNPs))) {
      SNPs <- metadata(SNPs)$SNPs[, colnames(SNPs)]
    } else {
      SNPs <- getBeta(SNPs)[grep("^rs", rownames(SNPs)), colnames(SNPs)]
    }
    stopifnot(ncol(SNPs) == NC)
  }

  # Could be parallelized if truly absurd sample sizes are being fit. 
  apply(SNPs, 2, .SNPcall)

}


# helper for heavy lifting
.SNPcall <- function(x) {

  calls <- rep(NA_integer_, length(x))
  names(calls) <- names(x)
  y <- x
  yy <- which(!is.na(y))
  x <- as.numeric(flogit(y[yy]))
  res <- suppressMessages(try(mixR::select(x=x, ncomp=1:3), silent=TRUE))
  if (!inherits(res, "try-error")) {
    k <- res$ncomp[which(res$best == "*")]
    SNPfit <- mixR::mixfit(x=x, ncomp=k)
    SNPprob <- SNPfit$comp.prob
    SNPclass <- apply(SNPprob, 1, which.max) - 1
    if (length(SNPclass) == length(x)) {
      names(SNPclass) <- names(x)
      calls[yy] <- SNPclass
    }
  }
  calls

}
