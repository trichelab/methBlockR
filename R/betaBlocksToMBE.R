#' read beta blocks (wgbstools beta_to_table output) into a MethBlockExperiment
#' 
#' @param tbx   either a TabixFile or something that can be turned into one
#' @param n     how many lines to read? (default is -1, i.e., all of them)
#' @param g     optional genome indicator (default is NULL) 
#'
#' @return      a MethBlockExperiment
#'
#' @details     realistically, this could be done with scanTabix, I just don't. 
#'
#' @import      Rsamtools
#'
#' @export 
#'
betaBlocksToMBE <- function(tbx, n=-1, g="hg19") {
  if (!is(tbx, "TabixFile")) tbx <- TabixFile(tbx)
  x <- readBetaBlocks(tbx, n=n, g=g)
  gr <- attr(x, "GRanges")
  attr(x, "GRanges") <- NULL # no need to duplicate it
  start(gr) <- start(gr) + 1 # quirk from beta_to_table
  MBE <- MethylationExperiment(assay=list(Beta=x), rowRanges=gr, 
                               colData=DataFrame(ID=colnames(x)),
                               metadata=list(sourceFile=tbx$path))
  class(MBE) <- "MethBlockExperiment"
  mainExpName(MBE) <- "CpG"
  return(MBE)
}
