#' plotME
#'
#' plot a MethylationExperiment or MethBlockExperiment features + SNPs
#'
#' @param x         a MethylationExperiment or MethBlockExperiment
#' @param k         how many features to use (100)
#' @param minmean   minimum mean for byExtremality (0.2)
#' @param maxmean   maximum mean for byExtremality (0.8)
#' @param maxNA     maximum fraction NA for a feature (0.2) 
#' @param ...       parameters to pass to Heatmap
#'
#' @details byExtremality() is called to determine features. A joint plot
#'          of highly variable features on the left and SNPs on the right
#'          is returned. This function is purely for a convenient first pass.
#'          Recent experiences with low-coverage WGBS have convinced us that
#'          bounding minimum and maximum mean is important for byExtremality.
#'
#' @import ComplexHeatmap
#' @import circlize
#'
#' @export
#' 
plotME <- function(x, k=100, minmean=0.2, maxmean=0.8, maxNA=0.2, ...) {

    chr <- intersect(seqlevels(x), paste0("chr", 1:22))
    b <- assay(keepSeqlevels(x, chr, pruning.mode="coarse"), "Beta")
  
    message("Finding highly extremal features...")
    N <- ncol(x)
    keepFeats <- rownames(b)[which(rowSums(is.na(b)) / N <= maxNA)]
    toPlot <- byExtremality(b[keepFeats,], k, minmean=minmean, maxmean=maxmean)
    message("Plotting DNA methylation...")
    jet <- colorRamp2(seq(0, 1, 0.125),
                      c("#00007F", "blue", "#007FFF", "cyan",
                        "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
    H1 <- Heatmap(as(t(toPlot), "matrix"), col=jet, name="Methylation", 
                  clustering_distance_columns="manhattan",
                  clustering_method_columns="ward.D2",
                  clustering_distance_rows="manhattan",
                  clustering_method_rows="ward.D2",
                  column_names_gp = gpar(fontsize=9),
                  row_names_side="left",
                  ...)
    if ("SNPs" %in% names(metadata(x))) { 
      message("Plotting SNPs...")
      H2 <- suppressMessages(try(plotSNPcalls(x, rotate=TRUE, 
                                              column_names_gp=gpar(fontsize=9)),
                                 silent=TRUE))
      if (!inherits(H2, "try-error")) { 
        return(H1 + H2)
      } else { 
        message("SNP calling failed, omitting SNP plot")
      }
    }
    H1

}
