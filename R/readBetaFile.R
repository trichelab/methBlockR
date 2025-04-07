#' read wgbstools .beta file into a matrix (CpG bed into attr(, "GRanges"))
#' 
#' @param x     path to the file
#' @param CpG   optional but HIGHLY recommended CpG dictionary for x (NULL)
#' @param g     optional but HIGHLY recommended genome indicator (NULL)
#'
#' @return      a matrix of methylated (column 1) and total (column 2) reads
#'
#' @details     the output is nearly useless without the CpG dictionary. 
#' 
#' @import      GenomicRanges
#'
#' @export 
#'
readBetaFile <- function(x, CpG=NULL, g=NULL) { 
  N <- file.info(x)$size
  x <- matrix(readBin(x, "integer", N, size=1, signed=0), N/2, 2, byrow=TRUE)
  if (!is.null(CpG)) { 
    gr <- import(CpG)
    rownames(x) <- as.character(gr)
    if (!is.null(g)) genome(gr) <- g
    attr(x, "GRanges") <- gr
  }
  return(x) 
}
