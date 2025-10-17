#' plotME
#'
#' plot a MethylationExperiment or MethBlockExperiment features + SNPs
#'
#' @param x         a MethylationExperiment or MethBlockExperiment
#' @param k         maximum number of features to use (100)
#' @param minmean   minimum mean for a row to be included (0.2)
#' @param maxmean   maximum mean for a row to be included (0.8)
#' @param maxNA     maximum fraction NA for a feature (0.2) 
#' @param splitBy   a column name on which to split rows (NULL)
#' @param rankBy    'extremality' or 'sd' (extremality)
#' @param pal       palette for SNP plot ("prp" or "jet") ("prp") 
#' @param ...       parameters to pass to Heatmap
#' @param BPPARAM   BiocParallelParam() to pass to plotSNPcalls (SerialParam())
#'
#' @details If SNPs are not found, only the main methylation plot is created.
#'          For reasons not yet clear, it's usually best to leave BPPARAM alone.
#'          For reasons quite clear, caching SNP calls is a very good idea, via
#'          metadata(x)$SNPcalls <- SNPcalls(x), and will speed up things A LOT
#'
#' @import ComplexHeatmap
#' @import BiocParallel
#' @import circlize
#'
#' @export
#' 
plotME <- function(x, k=100, minmean=0.2, maxmean=0.8, maxNA=0.2, splitBy=NULL, rankBy=c("extremality", "sd"), pal=c("prp", "jet"), ..., BPPARAM=NULL) {

    N <- ncol(x)
    pal <- match.arg(pal)
    rankBy <- match.arg(rankBy)
    chr <- intersect(seqlevels(x), paste0("chr", 1:22))
    b <- assay(keepSeqlevels(x, chr, pruning.mode="coarse"), "Beta")
    b <- b[which(rowSums(is.na(b)) / N <= maxNA), ]
    description <- c(extremality="extremal", sd="variable")

    message("Finding top ", k, " most ", description[rankBy], " features...")
    toPlot <-switch(rankBy, 
                    sd=bySD(b, k=k, min=minmean, max=maxmean),
                    extremality=byExtremality(b, k=k, min=minmean, max=maxmean))

    message("Plotting DNA methylation...")
    if (!is.null(splitBy)) { 
      H1 <- .plotMethSplit(toPlot, row_split = colData(x)[, splitBy], ...)
    } else {
      H1 <- .plotMethSplit(toPlot, ...)
    }
    if ("SNPs" %in% names(metadata(x))) { 
      if (is.null(BPPARAM)) BPPARAM <- SerialParam(progressbar = TRUE)
      message("plotting SNPs...")
      H2 <- suppressMessages(try(plotSNPcalls(x, 
                                              qc=TRUE,
                                              pal=pal,
                                              BPPARAM=BPPARAM,
                                              show_row_dend=FALSE),
                                 silent=TRUE))
      if (!inherits(H2, "try-error")) { 
        qcResult <- attr(H2@matrix, "qcResult")
        # bugfix for a rather bizarre situation
        H1 <- .plotMethSplit(toPlot[, rownames(H2@matrix)],
                             cluster_rows = FALSE) 
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


# helper fn
.plotMethSplit <- function(toPlot, ...) { 
    
  jet <- colorRamp2(seq(0, 1, 0.125),
                    c("#00007F", "blue", "#007FFF", 
                      "cyan", "#7FFF7F", "yellow", 
                      "#FF7F00", "red", "#7F0000"))
  row_size <- min(9, (60 / log2(nrow(toPlot))))
  col_size <- min(9, (60 / log2(ncol(toPlot))))
  Heatmap(as(t(toPlot), "matrix"), col=jet, name="Methylation", 
          column_names_gp = gpar(fontsize = col_size),
          row_names_gp = gpar(fontsize = row_size),
          clustering_distance_columns="manhattan",
          clustering_method_columns="ward.D2",
          clustering_distance_rows="manhattan",
          clustering_method_rows="ward.D2",
          show_parent_dend_line = FALSE,
          row_gap = unit(1, "mm"),
          row_names_side="left",
          row_title_rot = 0,
          ...)
}
