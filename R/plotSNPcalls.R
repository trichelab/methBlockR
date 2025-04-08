#' call and plot SNPs from a matrix of SNP betas 
#' 
#' plot SNPs living in metadata(x)$SNPs, provided they match sampleNames(x)
#' 
#' @param x       an object with SNPs in its metadata(), or a SNP matrix 
#' @param rotate  rotate the subjects onto the side? (FALSE)
#' @param qc      QC on the SNPs? (implies rotate) (FALSE)
#' @param ...     other arguments passed on to Heatmap
#' 
#' @details Plotting is done by .plotSNPcallmat() and .plotSNPheat()
#' 
#' @import ComplexHeatmap
#' @import circlize
#' 
#' @export
plotSNPcalls <- function(x, rotate=FALSE, qc=FALSE, ...) { 

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
  .plotSNPcallmat(SNPcalls, qc=qc...) 

}


# helper fn
.plotSNPcallmat <- function(SNPs_called, qc=FALSE, ...) {

  toUse <- which(missingness(SNPs_called, 1) < .5)
  if (qc) { 
    callRanges <- apply(colRanges(SNPs_called, na.rm=TRUE), 1, diff)
    qcColor <- c(fail="darkred", note="orange", pass="green")
    qcName <- c("fail", "note", "pass")
    qcRes <- factor(qcName[callRanges + 1])
    qcAnno <- rowAnnotation(qc = qcRes, col = list(qc = qcColor))
    .plotSNPheat(SNPs_called[toUse, ], right_annotation = qcAnno, ...)
  } else { 
    .plotSNPheat(SNPs_called[toUse, ], ...)
  }

}


# helper fn (DRY)
.plotSNPheat <- function(SNPs_called, ...) { 
  
  SNP_colors <- colorRamp2(seq(0, 2), c("#00007F", "yellow", "#7F0000"))
  Heatmap(SNPs_called,
          col=SNP_colors, 
          name="Alleles",
          clustering_method_columns='ward.D2', 
          clustering_distance_columns='manhattan', 
          clustering_method_rows='ward.D2', 
          clustering_distance_rows='manhattan',
          ...)

}
