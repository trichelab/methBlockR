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
#' @details SNPs will be read in first, so that any failures happen quickly. 
#'          BPPARAM=MulticoreParam(progressbar=TRUE) is a good starting point.
#'          We hope to eventually support FileSetExperiment as a backend,
#'          or some sort of materialized slice-able on-disk format.
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
  
  # if that succeeds, read beta values
  message("Processing beta values...")
  asys <- list(Beta=openSesame(sdfs, ..., BPPARAM=BPPARAM, func=getBetas))
  rd <- data.frame(row.names=rownames(asys[["Beta"]]))
  rd$type <- c(cg="CpG", ch="CpH", rs="SNP")[substr(rownames(rd), 1, 2)]

  # if requested, CN
  if (intensities) {
    message("Processing intensities for copy number analysis...")
    asys[["CN"]] <- openSesame(sdfs, func=totalIntensities, BPPARAM=BPPARAM)
  } 

  # construct the object & partition with splitAltExps
  ME <- MethylationExperiment(assays=asys, rowData=rd)
  ME <- splitAltExps(ME, rowData(ME)$type, ref="CpG")
  colnames(ME) <- basename(IDATs) # yes, it's a bit janky
  metadata(ME)[["SNPs"]] <- altExp(ME, "SNP") # for MBE()
  return(ME)

}
