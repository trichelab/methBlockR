#' The MethBlockExperiment class (inherits from ME, which inherits from SCE). 
#'
#' The difference between a MethBlockExperiment and a MethylationExperiment is
#' that the former (this class) represents mappings(blocks) of correlated 
#' DNA methylation regions for at least one genome assembly (usually 2 or 3). 
#' Since methylation rates tend to be more stable across regions than at 
#' individual loci, this can be useful for low-input or cross-platform 
#' comparison. On the other hand, if one wants to model 'epi-clones' or
#' similar stochastic process based representations, a MethylationExperiment
#' or perhaps bsseq or epiallele representation may be preferable. Use the 
#' right tool for the experimental question at hand. We find that methylation 
#' blocks are often the right tool for modeling cell states and fates. 
#'
#' @param x       a MethylationExperiment
#' @param ...     additional arguments to asMethBlocks
#'
#' @details       Can go straight from IDATs to blocks by chaining this. 
#'
#' @return  a MethBlockExperiment 
#'
#' @import methods
#'
#' @export
#'
MethBlockExperiment <- function(x, ...) {

  if ("SNP" %in% altExpNames(x)) { 
    metadata(x)$SNPs <- getBeta(altExp(x, "SNP")) # with colnames
  }
  
  # bit of a kludge, but oh well
  y <- asMethBlocks(as(removeAltExps(x), "SingleCellExperiment"), ...)
  class(y) <- "MethBlockExperiment"
  y

}
