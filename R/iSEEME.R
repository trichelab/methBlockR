#' convenience wrapper for iSEE on a MethylationExperiment
#'
#' Note that this is a quick and dirty affair so don't expect much. 
#'
#' @param x             a MethylationExperiment with NMF and UMAP reducedDims
#' @param colorColumn   name of the colData column to color plots ("cell.type")
#'
#' @details presumes NMF and UMAP. you will almost certainly want to tweak this.
#'
#' @import iSEE
#'
#' @export
#'
iSEEME <- function(x, colorColumn = "cell.type") { 

  # reason for this will become obvious 
  stopifnot(all(c("UMAP","NMF") %in% reducedDimNames(x)))
  if (!colorColumn %in% names(colData(x))) colorColumn <- "subtype"

  iSEE(x,
     initial = 
       list(UMAP = new("ReducedDimensionPlot",
                       FontSize = 1.5,
                       PointSize = 4,
                       VisualBoxOpen = TRUE,
                       ColorBy = "Column data",
                       ColorByColumnData = colorColumn,
                       Type = "UMAP"
                       ),
            NMF  = new("ColumnDataPlot",
                       YAxis = "nmf1",
                       XAxis = "Column data",
                       XAxisColumnData = colorColumn,
                       FontSize = 1.5,
                       PointSize = 4,
                       DataBoxOpen = TRUE,
                       SelectionBoxOpen = TRUE,
                       ColorBy = "Column data",
                       ColorByColumnData = colorColumn,
                       ColumnSelectionDynamicSource = TRUE,
                       ColumnSelectionSource = "ReducedDimensionPlot1"
                       ),
            Item = new("ColumnDataTable",
                       SelectionBoxOpen = TRUE,
                       ColumnSelectionDynamicSource = TRUE,
                       ColumnSelectionSource = "ReducedDimensionPlot1"
                       )
            )
     )
}
