#' plotME
#'
#' plot a MethylationExperiment or MethBlockExperiment features + SNPs
#'
#' @param x a MethylationExperiment or MethBlockExperiment
#' @param k how many features to use (500)
#'
#' @details byExtremality() is called to determine features. A joint plot
#'          of highly variable features on the left and SNPs on the right
#'          is returned. This function is purely for a convenient first pass.
#'
#' @import ComplexHeatmap
#' @import circlize
#'
#' @export
#' 
plotME <- function(x, k=500, ...) {

    chr <- intersect(seqlevels(x), paste0("chr", 1:22))
    b <- getBeta(keepSeqlevels(x, chr, pruning.mode="coarse"))
  
    message("Finding highly extremal features...")
    N <- ncol(x)
    keepFeats <- rownames(x)[which(rowSums(is.na(getBeta(x)))/N < 0.5)]
    toPlot <- byExtremality(b[keepFeats,], k)
    message("Plotting DNA methylation...")
    jet <- colorRamp2(seq(0, 1, 0.125),
                      c("#00007F", "blue", "#007FFF", "cyan",
                        "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
    H1 <- Heatmap(as(t(toPlot), "matrix"), col=jet, name="Methylation", 
                  clustering_distance_columns="manhattan",
                  clustering_method_columns="ward.D2",
                  clustering_distance_rows="manhattan",
                  clustering_method_rows="ward.D2",
                  ...)
    if ("SNPs" %in% names(metadata(x))) { 
      message("Plotting SNPs...")
      H2 <- plotSNPcalls(x, rotate=TRUE)
      H1 + H2
    } else { 
      H1
    } 

}
