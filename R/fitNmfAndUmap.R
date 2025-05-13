#' performs NMF on the beta values, stuffs the factors into colData, runs UMAP
#' 
#' @param x   a MethBlockExperiment or similar thing-with-Beta-values
#' @param k   number of factors to find (30) (FIXME: use ARD instead)
#' @param ... other arguments to pass to RcppML::nmf
#'
#' @return    SCE, with new reducedDim()s, colData, nmf_fit & umap_fit metadata
#'
#' @details   the RSpectra import is to avoid a weird umap bug.
#'
#' @import    SingleCellExperiment
#' @import    RSpectra
#' @import    RcppML
#' @import    Matrix
#' @import    uwot
#'
#' @export
#'
fitNmfAndUmap <- function(SCE, k=30, ...) {

  if (!is(SCE, "SingleCellExperiment")) { 
    SCE <- as(SCE, "SingleCellExperiment") # if a SummarizedExperiment
  }
  k <- min(k, ncol(SCE) - 1)
  # mask <- as(is.na(getBeta(SCE)), "dMatrix")
  message("Fitting NMF with ", k, " components... ", appendLF = FALSE)
  nmf_fit <- nmf(as(.NA2ZERO(squeeze(getBeta(SCE))), "dMatrix"), mask="zeros", 
                 k=k, L1=0.01, ...)
  message("done.")
  reducedDim(SCE, "NMF") <- as(t(nmf_fit@h), "matrix") # so UMAP doesn't croak
  cols <- colnames(reducedDim(SCE, "NMF"))
  message("Adding ", cols[1], "-", cols[length(cols)], " to colData()... ", 
          appendLF = FALSE)
  colData(SCE)[, cols] <- reducedDim(SCE, "NMF")
  message("done.")
  metadata(SCE)$nmf_fit <- nmf_fit
  message("Fitting UMAP on NMF hat matrix... ", appendLF = FALSE)
  umap_fit <- umap(reducedDim(SCE, "NMF"), 
                   n_neighbors=min(round(ncol(SCE)/2), 15),
                   metric="cosine", ret_model=TRUE)
  reducedDim(SCE, "UMAP") <- umap_fit$embedding
  metadata(SCE)$umap_fit <- umap_fit
  message("done.")
  return(SCE)

}


# helper fn
.NA2ZERO <- function(x) {
  ifelse(is.na(x), 0, x)
}
