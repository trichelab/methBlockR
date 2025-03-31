#' convenience wrapper for iSEE on a MethylationExperiment
#'
#' Note that this is a quick and dirty affair so don't expect much. 
#'
#' @param   x   a MethylationExperiment with NMF and UMAP reducedDims
#'
#' @details presumes NMF and UMAP. you will almost certainly want to tweak this.
#'
#' @import iSEE
#'
#' @export
#'
iSEEME <- function(x) { 

  # reason for this will become obvious 
  stopifnot(all(c("UMAP","NMF") %in% reducedDimNames(x)))

  iSEE(x,
     initial = 
       list(UMAP = new("ReducedDimensionPlot",
                       FontSize = 1.5,
                       PointSize = 4,
                       VisualBoxOpen = TRUE,
                       ColorBy = "Column data",
                       ColorByColumnData = "subtype", 
                       TooltipColumnData = c("Pheno"),
                       Type = "UMAP"
                       ),
            NMF  = new("ColumnDataPlot",
                       YAxis = "nmf1",
                       XAxis = "Column data",
                       XAxisColumnData = "subtype",
                       FontSize = 1.5,
                       PointSize = 4,
                       DataBoxOpen = TRUE,
                       SelectionBoxOpen = TRUE,
                       ColorBy = "Column data",
                       ColorByColumnData = "subtype",
                       ColumnSelectionDynamicSource = TRUE,
                       ColumnSelectionSource = "ReducedDimensionPlot1",
                       TooltipColumnData = c("Pheno")
                       ),
            Item = new("ColumnDataTable",
                       SelectionBoxOpen = TRUE,
                       ColumnSelectionDynamicSource = TRUE,
                       ColumnSelectionSource = "ReducedDimensionPlot1"
                       )
            )
     )
}
