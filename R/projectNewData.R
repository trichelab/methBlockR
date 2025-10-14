#' project new data into an existing NMF model using RcppML::project
#' 
#' @param x   a MethBlockExperiment, SingleCellExperiment, or nmf object (query)
#' @param y   a MethBlockExperiment, SingleCellExperiment, or nmf object (model)
#' @param how how to project new data into the existing space? (both, nmf, umap)
#' @param ... other arguments to pass to RcppML::project and/or uwot::umap
#'
#' @return    a ReducedExperiment or FactorisedExperiment (see Details)
#'
#' @details   the type of returned object depends on "how": if "both" (the 
#'            default value), a ReducedExperiment with UMAP reduction will be 
#'            returned with a FactorisedExperiment in metadata()$nmf_projection.
#'            This provides both a UMAP projection and the NMF W and H matrices.
#'            If how == "nmf", only the FactorisedExperiment will be returned.
#'            If how == "umap", only the ReducedExperiment will be returned. 
#'            If we start supporting more methods, we'll figure something out.
#'
#' @import    SingleCellExperiment
#' @import    ReducedExperiment
#' @import    RcppML
#' @import    uwot
#'
#' @export
#'
projectNewData <- function(x, y, how=c("both", "nmf", "umap"), ...) {

  how <- match.arg(how)
  if (how %in% c("nmf", "both")) stopifnot("nmf_fit" %in% names(metadata(y)))
  if (how %in% c("umap", "both")) stopifnot("umap_fit" %in% names(metadata(y)))
  stopifnot(is(x, "SingleCellExperiment") | is(x, "matrix") | is(x, "Matrix"))
  stopifnot(.checkDims(x, y, how=how))
  stop("Need to finish projectNewData!")

}


# helper fn
.checkDims <- function(x, y, how = c("both", "nmf", "umap")) { 

  # determine if x is compatible with y based on NMF and/or UMAP feature vectors
  how <- match.arg(how)
  switch(how, 
         both=all(rownames(x) %in% rownames(metadata(y)$nmf_fit)),
         umap=(nrow(x) == metadata(y)$umap_fit$norig_col),
         nmf=all(rownames(x) %in% rownames(metadata(y)$nmf_fit)))

}
