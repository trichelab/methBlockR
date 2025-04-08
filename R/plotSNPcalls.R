#' call and plot SNPs from a matrix of SNP betas 
#' 
#' plot SNPs living in metadata(x)$SNPs, provided they match sampleNames(x)
#' 
#' @param x       an object with SNPs in its metadata(), or a SNP matrix 
#' @param rotate  rotate the subjects onto the side? (FALSE)
#' @param qc      QC on the SNPs? (implies rotate) (FALSE)
#' @param BPPARAM a BiocParallelParam() object, or (default) SerialParam()
#' @param ...     other arguments passed on to Heatmap
#' 
#' @details Plotting is done by .plotSNPcallmat() and .plotSNPheat()
#'          If SNP calls are found in metadata(x)$SNPcalls, they will be used,
#'          otherwise SNPs will be called on the fly (which can be sloooow). 
#'
#' @import ComplexHeatmap
#' @import circlize
#' 
#' @export
plotSNPcalls <- function(x, rotate=FALSE, qc=FALSE, BPPARAM=NULL, ...) { 

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
  if (is.null(BPPARAM)) BPPARAM <- SerialParam(progressbar = TRUE)
  if (is(x, "SummarizedExperiment") & "SNPcalls" %in% names(metadata(x))) { 
    message("Cached SNP calls found in metadata(x)$SNPcalls (yay)")
    calls <- metadata(x)$SNPcalls
    qcCol <- attr(calls, "qcColors")
    qcRes <- attr(calls, "qcResult")[colnames(x)]
    calls <- calls[, colnames(x)]
  } else {
    message("You can cache SNP calls: metadata(x)$SNPcalls <- SNPcalls(x)")
    calls <- SNPcalls(SNPs, BPPARAM=BPPARAM) # qc results added here
    qcRes <- attr(calls, "qcResult")
    qcCol <- attr(calls, "qcColors")
  }
  calls <- calls[which(missingness(calls, 1) < .5), ]
  if (rotate | qc) .plotSNPcallmat(t(calls), qc=qc,qcRes=qcRes,qcCol=qcCol,...) 
  else .plotSNPcallmat(calls, ...) 

}


# helper fn
.plotSNPcallmat <- function(calls, qc=FALSE, qcRes=NULL, qcCol=NULL, ...) {

  if (qc) {
  
    # avoid issues due to t() 
    attr(calls, "qcResult") <- qcRes
    attr(calls, "qcColors") <- qcCol
    if (any(qcRes != "pass")) {
      flag <- which(qcRes != "pass")
      flagged_gp <- gpar(fontsize = min(12, 40 / log2(length(flag))))
      flagged <- paste(names(qcRes)[flag], qcRes[flag], sep=": ")
      qcAnno <- rowAnnotation(qc = qcRes, 
                              link = anno_mark(at = flag, 
                                               labels = flagged,
                                               labels_gp = flagged_gp),
                              col = list(qc = qcCol, 
                                         link = qcCol))
      .plotSNPheat(calls, 
                   row_split = qcRes, 
                   cluster_row_slices = FALSE,
                   show_parent_dend_line = FALSE,
                   right_annotation = qcAnno,
                   show_row_names = FALSE, 
                   row_title_rot = 0,
                   ...)
    } else { 
      qcAnno <- rowAnnotation(qc = qcRes, col = list(qc = qcCol))
      .plotSNPheat(calls, right_annotation = qcAnno, ...)
    }
  } else { 
    .plotSNPheat(calls, ...)
  }

}


# helper fn (DRY)
.plotSNPheat <- function(calls, ...) { 
  
  Heatmap(calls,
          name="Alleles",
          clustering_method_columns='ward.D2', 
          clustering_distance_columns='manhattan', 
          column_names_gp = gpar(fontsize = min(9, 40 / log2(ncol(calls)))),
          clustering_method_rows='ward.D2', 
          clustering_distance_rows='manhattan',
          row_names_gp = gpar(fontsize = min(9, 60 / log2(nrow(calls)))),
          col=colorRamp2(seq(0, 2), c("#00007F", "yellow", "#7F0000")),
          ...)

}
