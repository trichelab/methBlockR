#' performs NMF on the beta values, stuffs the factors into colData, runs UMAP
#' 
#' @param SCE a SingleCellExperiment (or SummarizedExperiment to coerce to one)
#' @param k   number of factors to find (20)
#' @param ... other arguments to pass to RcppML::nmf_fit
#'
#' @return    SCE, with new reducedDim()s, colData, nmf_fit & umap_fit metadata
#'
#' @import    SingleCellExperiment
#' @import    methBlockR 
#' @import    RSpectra
#' @import    RcppML
#' @import    uwot
#'
#' @export
#'
fitNmfAndUmap <- function(SCE, k=20, ...) {

  .NA2ZERO <- function(x) ifelse(is.na(x), 0, x) # mask=NA no worky
  if (!is(SCE, "SingleCellExperiment")) { 
    SCE <- as(SCE, "SingleCellExperiment") # if a SummarizedExperiment
  }
  message("Fitting NMF with ", k, " components... ", appendLF = FALSE)
  nmf_fit <- nmf(.NA2ZERO(getBeta(SCE)), mask = "zeros", k = k, ...,
                 tol = 1e-3, L1 = 0.01, maxit = 1000, seed = 1234)
  message("done.")
  reducedDim(SCE, "NMF") <- as(t(nmf_fit@h), "matrix") # so UMAP doesn't croak
  cols <- colnames(reducedDim(SCE, "NMF"))
  message("Adding ", cols[1], "-", cols[length(cols)], " to colData()... ", 
          appendLF = FALSE)
  colData(SCE)[, cols] <- reducedDim(SCE, "NMF")
  message("done.")
  metadata(SCE)$nmf_fit <- nmf_fit
  message("Fitting UMAP on NMF hat matrix... ", appendLF = FALSE)
  umap_fit <- umap(reducedDim(SCE, "NMF"), metric="cosine", ret_model=TRUE)
  reducedDim(SCE, "UMAP") <- umap_fit$embedding
  metadata(SCE)$umap_fit <- umap_fit
  message("done.")
  return(SCE)

}
