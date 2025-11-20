#' extract the positional weights for an NMF factor
#' 
#' @param x   a MethylationBlockExperiment, usually 
#' @param y   which factor (1)
#' @param g   which genome (whichever one x claims to be using)
#' 
#' @return    a GRanges where the $score column is that factor's weights
#' 
#' @details   like duh, if there's no mcols(x)[["NMF"]], this won't work.
#'            You will want to use rtracklayer or igvR to view the result;
#'            seqinfo(gr) <- SeqinfoForUCSCGenome(g)[seqlevels(gr)] may help
#' 
#' @export
#'
factorWeight <- function(x, y=1, g=NULL) { 

  stopifnot("NMF" %in% names(mcols(x)))
  stopifnot(y <= ncol(mcols(x)[["NMF"]]))
  if (is.null(g)) g <- unique(genome(x))
  stopifnot(g %in% names(mcols(x)))

  gr <- as(mcols(x)[[g]], "GRanges")
  gr$score <- mcols(x)[["NMF"]][, y]
  genome(gr) <- g
  return(gr)

}
