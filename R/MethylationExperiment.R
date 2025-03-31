#' MethylationExperiment() constructor function 
#'
#' @param ...   arguments passed to the SingleCellExperiment() constructor
#'
#' @return  a MethylationExperiment
#'
#' @details If openSesameToME is used, a matrix of total intensities can 
#'          be added to assays(object, "CN") with per-probe intensities.
#'          SingleCellExperiment::splitAltExps separates CpG, CpH & SNP probes. 
#'
#' @seealso     openSesameToME
#'
#' @importFrom  GenomicRanges         GRangesList
#' @importFrom  RaggedExperiment      RaggedExperiment
#' @importFrom  SingleCellExperiment  SingleCellExperiment
#'
#' @import      methods
#'
#' @rdname      MethylationExperiment
#'
#' @export
#'
MethylationExperiment <- function(...) {
  
  y <- SingleCellExperiment(...)
  class(y) <- "MethylationExperiment"
  y

}


#' Construct a MethylationExperiment from IDATs.
#'
#' @param IDATs       a vector of IDAT stubs (Basenames) with paths to them
#' @param BPPARAM     a BiocParallel::MulticoreParam() or similar for processing
#' @param intensities process total intensities into assay(x, "CN")? (TRUE)
#' @param cd          optional colData; rownames must match colnames of betas
#' @param ...         other arguments passed to sesame::openSesame()
#'
#' @details BPPARAM=MulticoreParam(progressbar=TRUE) is usually a good idea. 
#'          In the event that multiple array formats are detected, the function
#'          will attempt to rationalize the output by injecting NAs into the 
#'          probe-level output from the smaller arrays, and raising a warning.
#'
#' @seealso BiocParallel::MulticoreParam
#' @seealso MethylationExperiment
#' @seealso sesame::openSesame
#'
#' @importFrom  sesame openSesame
#' @import BiocParallel
#'
#' @export
#'
openSesameToME <- function(IDATs, BPPARAM=NULL, intensities=TRUE, cd=NULL, ...){

  # BPPARAM can default to a SerialParam, although that's dumb
  if (is.null(BPPARAM)) BPPARAM <- SerialParam(progressbar=TRUE)

  # process the IDAT files 
  message("Reading IDATs...")
  sdfs <- bplapply(IDATs, BPPARAM=BPPARAM, readIDATpair)
  rows <- vapply(sdfs, nrow, 1L) 
  if (length(unique(rows)) > 1) warning("Multiple array types detected.")

  # if that succeeds, read beta values
  message("Processing beta values...")
  Beta <- openSesame(sdfs, ..., BPPARAM=BPPARAM, func=getBetas)
  probes <- names(Beta[[which.max(rows)]])

  # wrap this better (DRY)
  if (length(unique(rows)) > 1) {
    message("Merging platforms...")
    Beta <- do.call(cbind, 
                    lapply(Beta, 
                           function(b) {
                             bb <- rep(NA_real_, length(probes))
                             names(bb) <- probes
                             pb <- intersect(probes, names(b))
                             bb[pb] <- b[pb]
                             return(bb)
                           }))
  }
  asys <- list(Beta=Beta)
  rd <- data.frame(row.names=rownames(asys[["Beta"]]))
  rd$type <- c(cg="CpG", ch="CpH", rs="SNP")[substr(rownames(rd), 1, 2)]

  # if requested, CN
  if (intensities) {
    message("Processing intensities for copy number analysis...")
    CN <- openSesame(sdfs, func=totalIntensities, BPPARAM=BPPARAM)
    if (length(unique(rows)) > 1) {
      message("Merging platforms...")
      CN <- do.call(cbind,
                    lapply(CN, 
                           function(i) {
                             ii <- rep(NA_real_, length(probes))
                             names(ii) <- probes
                             p <- intersect(probes, names(i))
                             ii[p] <- i[p]
                             return(ii)
                           }))
    }
    asys[["CN"]] <- CN
  }

  # construct the object & partition with splitAltExps
  ME <- MethylationExperiment(assays=asys, rowData=rd)
  ME <- splitAltExps(ME, rowData(ME)$type, ref="CpG")
  colnames(ME) <- basename(IDATs) # yes, it's a bit janky
  metadata(ME)[["SNPs"]] <- getBeta(altExp(ME, "SNP")) # for MBE()
  return(ME)

}


# convenience methods
setMethod("plot", c("MethylationExperiment", "missing"), 
  function(x, y, ...) plot(x, 500L, ...))
  
setMethod("plot", c("MethylationExperiment", "numeric"), 
  function(x, y, ...) {

    k <- intersect(seqlevels(x), paste0("chr", 1:22))
    b <- getBeta(keepSeqlevels(x, k, pruning.mode="coarse"))
    toPlot <- byExtremality(b, y)
    
    jet <- colorRamp2(seq(0, 1, 0.125),
                      c("#00007F", "blue", "#007FFF", "cyan",
                        "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
    H1 <- Heatmap(as(toPlot, "matrix"), col=jet, name="Methylation", 
                  clustering_distance_columns="manhattan",
                  clustering_method_columns="ward.D2",
                  clustering_distance_rows="manhattan",
                  clustering_method_rows="ward.D2",
                  ...)
    H2 <- plotSNPcalls(x, rotate=TRUE)
    H1 + H2

})
