#' read a tabix'ed beta block file from wgbstools, where columns 1-5 are coords
#'
#' @param tbx   either a TabixFile or something that can be turned into one
#' @param n     how many lines to read? (default is -1, i.e., all of them)
#' @param g     optional genome indicator (default is NULL) 
#'
#' @return      a matrix, possibly with a GRanges in attr(, "GRanges")
#'
#' @details     realistically, this could be done with scanTabix, I just don't. 
#'
#' @import      GenomicRanges
#' @import      Rsamtools
#'
#' @export 
#'
readBetaBlocks <- function(tbx, n=-1, g=NULL) {
  
  # first row is column names, usually either skipped or commented 
  cols <- sub("^(#)+", "", strsplit(readLines(tbx$path, n=1), "\t")[[1]])

  # the following rows are matrix columns, possibly preceded by coordinates
  df <- data.frame(do.call(rbind, strsplit(readLines(tbx$path, n=n+1)[-1], "\t")))
  names(df) <- cols
 
  # this is where, if coordinates are present, we extract them 
  cc <- grep("^(#)?(chr|start|end)(CpG)?$", names(df), value=TRUE)
  if (length(cc) > 0) { 
    message("Coordinates detected; placing a GRanges in attr(output, 'GRanges')")
    gr <- as(df[, cc], "GRanges")
    names(gr) <- as.character(gr) 
    for (i in grep("^(start|end)CpG$", names(mcols(gr)), value=TRUE)) {
      mcols(gr)[[i]] <- as.integer(mcols(gr)[[i]])
    }
    df <- df[, setdiff(names(df), cc)]
    rownames(df) <- names(gr) # to match coords
    if (!is.null(g)) genome(gr) <- g
  } else gr <- NULL 
  
  x <- data.matrix(df)/100 # quirk
  if (!is.null(gr)) attr(x, "GRanges") <- gr
  if (!(max(x, na.rm=TRUE) <= 1 & min(x, na.rm=TRUE) >= 0)) {
    warning("Values outside of 0-1 detected; squashing") 
    x[which(x > 1)] <- 1
    x[which(x < 0)] <- 0
  }
  message("Recovered ", nrow(x), " rows of ", ncol(x), " columns each.")
  fn <- fivenum(x)
  message("Median: ",fn[3],"; IQR ",fn[2],"-",fn[4],"; range ",fn[1],"-",fn[5])
  return(x)
    
}
