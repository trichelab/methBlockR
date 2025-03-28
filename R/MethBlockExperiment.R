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
#' Note that genome(MBE) <- 'hg38' or similar will attempt remapping of rows.
#' If this fails, the original (un-remapped) object will be returned instead.
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

  # drop non-CpG probes; keep SNP barcodes
  x <- unsplitAltExps(x, prefix.rows=FALSE)
  y <- asMethBlocks(x[grep("^cg", rownames(x)), ], ...)
  class(y) <- "MethBlockExperiment"
  y

}


#' @rdname  MethBlockExperiment
#'
#' @importMethodsFrom GenomeInfoDb "genome<-"
#'
#' @export 
#'
setReplaceMethod("genome", "MethBlockExperiment", 
  function(x, value) {

    y <- try(switchMethBlocksGenome(x, value))
    if (!inherits(y, "try-error")) return(y)
    else x

  }
)
