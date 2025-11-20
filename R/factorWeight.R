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
#'            (see example for a demonstration of how this approach can work)
#' 
#' @examples
#' \dontrun{
#'   nmf5 <- factorWeight(MBE, 5)
#'   nmf7 <- factorWeight(MBE, 7) 
#'   library(rtracklayer) 
#'   for (i in c("nmf5","nmf7")) {
#'     j <- get(i)
#'     seqinfo(j) <- SeqinfoForUCSCGenome(unique(genome(j)))[seqlevels(j)]
#'     export(j, paste(i, unique(genome(j)), "bw", sep="."))
#'   }
#' }
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
