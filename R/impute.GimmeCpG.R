#' impute missing sites via GimmeCpG (binned within 35bp, weighted out to 1kbp) 
#'
#' @param x       a MethylationExperiment (or something with getBeta and coords)
#' @param g       a genome (hg19 is the default, or use coords from granges(x)) 
#' @param BPPARAM a BiocParallelParam object (default SerialParam(progr=TRUE)) 
#' 
#' @return        x, but with a new assay Imputed that has been imputed.
#' 
#' @details       this function has not yet been tested enough for prime time.
#'
#' @export
impute.GimmeCpG <- function(x, g=NULL, BPPARAM=SerialParam(progressbar=TRUE)) {

  if (!is(x, "SummarizedExperiment") | !("Beta" %in% assayNames(x))) { 
    stop("Could not find beta values to impute")
  }

  if (is.null(g) & is.na(unique(genome(x)))) { 
    warning("No genome found, assuming hg19... this is risky") 
    g <- "hg19"
  } else if (!is.na(unique(genome(x)))) { 
    g <- unique(genome(x))
  }

  stop("Not completed yet")
  # need to add in +/- 1kb hits for imputable CpGs aren't typically masked

}


# helper fn; use precomputed whenever possible; this should be done in C++
.dw <- function(b1, d1, b2, d2) ((b1*abs(d1))+(b2*abs(d2)))/(abs(d1)+abs(d2))
