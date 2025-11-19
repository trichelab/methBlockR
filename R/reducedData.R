#' grab a data.frame of a reducedDim representation, name columns, bind colData
#' 
#' @param SCE     a SingleCellExperiment or something that acts like one 
#' @param name    the name of a reducedDim in SCE (will be used for names)
#' @param withcd  bind colData to the result? (TRUE) 
#' 
#' @return        a data.frame for feeding to ggplot2 or similar
#'
#' @export 
#'
reducedData <- function(SCE, name, withcd=TRUE) { 
 
  stopifnot(name %in% reducedDimNames(SCE))
  rd <- as.data.frame(reducedDim(SCE, name))
  names(rd) <- paste0(name, seq_len(ncol(rd))) 
  if (withcd) rd <- cbind(rd, as(colData(SCE), "data.frame"))
  return(rd)

}
