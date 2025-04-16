#' read wgbstools .beta file into a matrix (CpG bed into attr(, "GRanges"))
#' 
#' @param x     path to the file
#' @param CpG   optional but HIGHLY recommended CpG dictionary for x (NULL)
#' @param g     optional but HIGHLY recommended genome indicator (NULL)
#'
#' @return      a matrix of methylated (column 1) and total (column 2) reads
#'
#' @details     the output is nearly useless without the CpG dictionary. If the
#'              size of the file is 56434896, we know g == "hg19", but the CpG
#'              dictionary for hg19 is 125MB, so we don't have it in data(). 
#' 
#' @import      GenomicRanges
#' @import      rtracklayer
#'
#' @export 
#'
readBetaFile <- function(x, CpG=NULL, g=NULL) { 
  N <- file.info(x)$size
  x <- matrix(readBin(x, "integer", N, size=1, signed=0), N/2, 2, byrow=TRUE)
  if (nrow(x) == 28217448) {
    message("This is an hg19 .beta file: 28217448 CpGs")
    g <- "hg19"
  }
  if (!is.null(CpG)) { 
    gr <- import(CpG)
    rownames(x) <- as.character(gr)
    if (!is.null(g)) genome(gr) <- g
    attr(x, "GRanges") <- gr
  }
  return(x) 
}
