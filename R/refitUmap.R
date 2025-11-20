#' refit a UMAP on a subset of NMF columns and stuff it in new reducedDim
#' 
#' @param x     a MethBlockExperiment or similar thing-with-Beta-values
#' @param k_sub which NMF factors to use for recomputing the UMAP
#' @param name  add the weights to rowData(x)? (TRUE) 
#' @param ...   other arguments to pass to umap
#'
#' @return    SCE, with new reducedDim()s, and `name` metadata (model)
#'
#' @details   the RSpectra import is to avoid a weird umap bug.
#'
#' @import    SingleCellExperiment
#' @import    RSpectra
#' @import    Matrix
#' @import    uwot
#'
#' @export
#'
refitUmap <- function(SCE, k_sub, name) { 

  if (!is(SCE, "SingleCellExperiment")) { 
    SCE <- as(SCE, "SingleCellExperiment") # if a SummarizedExperiment
  }

  stopifnot(all(k_sub %in% seq_len(ncol(reducedDim(SCE, "NMF")))))
  if (identical(sort(k_sub), seq_len(ncol(reducedDim(SCE, "NMF"))))) { 
    warning("You selected all of the existing NMF dimensions!")
  }

  message("Fitting UMAP on subsetted NMF hat matrix... ", appendLF = FALSE)
  umap_fit <- umap(reducedDim(SCE, "NMF")[, k_sub], 
                   n_neighbors=min(round(ncol(SCE)/2), 15),
                   metric="cosine", ret_model=TRUE)
  reducedDim(SCE, name) <- umap_fit$embedding
  metadata(SCE)[[name]] <- umap_fit
  message("done.")
  return(SCE)

}
