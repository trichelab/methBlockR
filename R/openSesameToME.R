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
