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
