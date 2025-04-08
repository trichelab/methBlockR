#' fit a mixture model to flogit(SNPbetas) for each sample
#' 
#' @param SNPs      a SummarizedExperiment or a matrix of SNP beta values 
#' @param BPPARAM   a BiocParallelParam() object, or (default) SerialParam()
#'
#' @return          0/1/2 calls for each SNP in each sample
#' 
#' @details         mixR::mixfit(flogit(x), ncomp=1:3) is run for each sample. 
#'                  if the mixture fit fails, non-NA loci will be called as 0.
#'                  Inexplicably, the default SerialParam() is usually fastest.
#' 
#' @import          BiocParallel
#' @import          mixR
#' 
#' @export
SNPcalls <- function(SNPs, BPPARAM=NULL, ...) {

  if (is(SNPs, "SummarizedExperiment")) {
    NC <- ncol(SNPs)
    if ("SNPs" %in% names(metadata(SNPs))) {
      SNPs <- metadata(SNPs)$SNPs[, colnames(SNPs)]
    } else {
      SNPs <- getBeta(SNPs)[grep("^rs", rownames(SNPs)), colnames(SNPs)]
    }
    stopifnot(ncol(SNPs) == NC)
  }
  
  rsids <- rownames(SNPs)
  byspecimen <- split(t(SNPs), colnames(SNPs))
  for (i in names(byspecimen)) attr(byspecimen[[i]], "ID") <- i
  if (is.null(BPPARAM)) BPPARAM <- SerialParam(progressbar = TRUE)
  # could move QC evaluations here and add to attrs(calls) 
  message("Calling SNPs for ", ncol(SNPs), " samples...") 
  res <- bplapply(byspecimen, .SNPcall, BPPARAM=BPPARAM)
  calls <- do.call(cbind, res)
  rownames(calls) <- rsids
  qcNames <- c("fail", "note", "pass")
  callRanges <- apply(colRanges(calls, na.rm=TRUE), 1, diff)
  qcResult <- factor(qcNames[callRanges + 1])
  names(qcResult) <- colnames(calls)
  attr(calls, "qcResult") <- qcResult
  attr(calls, "qcColors") <- c(fail="darkred", note="orange", pass="green")
  return(calls)

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
    ncomp <- res$ncomp[which(res$best == "*")]
    SNPfit <- suppressMessages(try(mixR::mixfit(x=x, ncomp=ncomp), silent=TRUE))
    if (!inherits(SNPfit, "try-error")) {
      SNPclass <- try(apply(SNPfit$comp.prob, 1, which.max), silent=TRUE)
      if (inherits(SNPclass, "try-error")) {
        warning("SNP calling failed for ", attr(y, "ID"))
        calls[yy] <- 0
      } else { 
        if (length(SNPclass) == length(yy)) {
          names(SNPclass) <- names(x)[yy]
          calls[yy] <- SNPclass - 1
        } else { 
          calls[yy] <- 0
        }
      }
    } else {
      calls[yy] <- 0
    }
  } else { 
    calls[yy] <- 0
  }
  calls

}
