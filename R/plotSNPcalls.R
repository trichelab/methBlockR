#' call and plot SNPs from a matrix of SNP betas 
#' 
#' plot SNPs living in metadata(x)$SNPs, provided they match sampleNames(x)
#' 
#' @param x       an object with SNPs in its metadata(), or a SNP matrix 
#' @param rotate  rotate the subjects onto the side? (FALSE)
#' @param ...     other arguments passed on to Heatmap
#' 
#' @details The plotting itself is done by .plotSNPcallmat()
#' 
#' @import ComplexHeatmap
#' @import circlize
#' 
#' @export
plotSNPcalls <- function(x, rotate=FALSE, ...) { 

  if (is(x, "SummarizedExperiment")) {
    if ("SNPs" %in% names(metadata(x))) {
      SNPs <- metadata(x)$SNPs[, colnames(x)]
    } else {
      SNPs <- getBeta(x)[grep("^rs", rownames(x)), colnames(x)]
    }
    stopifnot(ncol(SNPs) == ncol(x))
  } else { 
    SNPs <- as.matrix(x)
  }
  SNPcalls <- SNPcalls(SNPs)
  if (rotate) SNPcalls <- t(SNPcalls)
  .plotSNPcallmat(SNPcalls, ...) 

}


# helper fn
.plotSNPcallmat <- function(SNPs_called, ...) {

  SNP_colors <- colorRamp2(seq(0, 2), c("#00007F", "yellow", "#7F0000"))
  toUse <- which(missingness(SNPs_called, 1) < .5)
  
  Heatmap(SNPs_called[toUse, ],
          col=SNP_colors, 
          name="Alleles",
          clustering_method_columns='ward.D2', 
          clustering_distance_columns='manhattan', 
          clustering_method_rows='ward.D2', 
          clustering_distance_rows='manhattan',
          ...)

}
