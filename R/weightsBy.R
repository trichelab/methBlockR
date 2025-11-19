#' plot the relative representation of some grouping by reducedData factor
#'
#' @param SCE   something that quacks like a SingleCellExperiment
#' @param by    name of a variable (usually a factor) in colData 
#' @param name  name of a reducedDim in SCE ("NMF")
#' @param add   optional confounder to include for e.g. faceting (NULL)
#' 
#' @return      a ggplot object
#'
#' @import reshape2
#' @import ggplot2 
#' 
#' @export
#'
weightsBy <- function(SCE, by, name="NMF", add=NULL) { 
 
  stopifnot(by %in% names(colData(SCE)))
  if (!is.null(add)) {
    id <- c("by", "add")
    stopifnot(add %in% names(colData(SCE)))
    wts <- cbind(reducedData(SCE, name, withcd=FALSE), 
                 by=colData(SCE)[, by], add=colData(SCE)[, add])
  } else { 
    id <- c("by")
    wts <- cbind(reducedData(SCE, name, withcd=FALSE), by=colData(SCE)[, by])
  }
  toPlot <- melt(wts, id.vars=id, variable.name="factor", value.name="weight")
  ggplot(toPlot) + 
    aes(x=weight, y=factor, fill=by) + 
    geom_bar(position="fill", stat="identity") + 
    theme_classic() 

}
