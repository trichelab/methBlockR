#' plotME
#'
#' plot a MethylationExperiment or MethBlockExperiment features + SNPs
#'
#' @param x         a MethylationExperiment or MethBlockExperiment
#' @param k         how many features to use (100)
#' @param minmean   minimum mean for byExtremality (0.2)
#' @param maxmean   maximum mean for byExtremality (0.8)
#' @param maxNA     maximum fraction NA for a feature (0.2) 
#' @param BPPARAM   BiocParallelParam() to pass to plotSNPcalls (SerialParam())
#' @param ...       parameters to pass to Heatmap
#'
#' @details byExtremality() is called to determine features. A joint plot
#'          of highly variable features on the left and SNPs on the right
#'          is returned. This function is purely for a convenient first pass.
#'          Recent experiences with low-coverage WGBS have convinced us that
#'          bounding minimum and maximum mean is important for byExtremality.
#'          If SNPs cannot be found, only the main methylation plot is returned.
#'          For reasons not yet clear, it's usually best to leave BPPARAM alone.
#'          For reasons quite clear, caching SNP calls is a very good idea, via
#'          metadata(x)$SNPcalls <- SNPcalls(x) # will speed up things A LOT
#'
#' @import ComplexHeatmap
#' @import BiocParallel
#' @import circlize
#'
#' @export
#' 
plotME <- function(x, k=100, minmean=0.2, maxmean=0.8, maxNA=0.2, BPPARAM=NULL, ...) {

    chr <- intersect(seqlevels(x), paste0("chr", 1:22))
    b <- assay(keepSeqlevels(x, chr, pruning.mode="coarse"), "Beta")
    message("Finding extremal features...")
    N <- ncol(x)
    keepFeats <- rownames(b)[which(rowSums(is.na(b)) / N <= maxNA)]
    toPlot <- byExtremality(b[keepFeats,], k, minmean=minmean, maxmean=maxmean)
    message("Plotting DNA methylation...")
    row_lab_size <- 80 / log2(nrow(x))
    col_lab_size <- 60 / log2(ncol(x))
    jet <- colorRamp2(seq(0, 1, 0.125),
                      c("#00007F", "blue", "#007FFF", "cyan",
                        "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
    H1 <- Heatmap(as(t(toPlot), "matrix"), col=jet, name="Methylation", 
                  clustering_distance_columns="manhattan",
                  clustering_method_columns="ward.D2",
                  clustering_distance_rows="manhattan",
                  clustering_method_rows="ward.D2",
                  show_parent_dend_line = FALSE,
                  column_names_gp = gpar(fontsize = col_lab_size),
                  row_names_gp = gpar(fontsize = row_lab_size),
                  row_names_side="left",
                  row_title_rot = 0,
                  row_gap = unit(1.5, "mm"),
                  ...)
    if ("SNPs" %in% names(metadata(x))) { 
      if (is.null(BPPARAM)) BPPARAM <- SerialParam(progressbar = TRUE)
      message("Calling and plotting SNPs...")
      H2 <- suppressMessages(try(plotSNPcalls(x, 
                                              qc=TRUE,
                                              BPPARAM=BPPARAM, 
                                              show_row_dend=FALSE),
                                 silent=TRUE))
      if (!inherits(H2, "try-error")) { 
        qcResult <- attr(H2@matrix, "qcResult")
        lst <- H1 + H2 
        draw(lst, row_split = qcResult)
      } else { 
        message("SNP calling failed, omitting SNP plot")
        draw(H1)
      }
    } else {
      draw(H1)
    }

}
