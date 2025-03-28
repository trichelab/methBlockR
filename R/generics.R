if (!isGeneric("getBeta")) 
  setGeneric("getBeta", function(object, ...) standardGeneric("getBeta"))

if (!isGeneric("getM")) 
  setGeneric("getM", function(object, ...) standardGeneric("getM"))

if (!isGeneric("getCN")) 
  setGeneric("getCN", function(object, ...) standardGeneric("getCN"))

if (!isGeneric("getSNPs")) 
  setGeneric("getSNPs", function(object, ...) standardGeneric("getSNPs"))


# helper function
.B2M <- function(p) log2(p / (1 - p))


#' @export
#'
setMethod("getBeta", "MethylationExperiment", 
          function(object, ...) assay(object, "Beta"))

#' @export
#'
setMethod("getM", "MethylationExperiment", 
          function(object, ...) .B2M(assay(object, "Beta")))

#' @export
#'
setMethod("getCN", "SummarizedExperiment", 
          function(object, ...) assay(object, "CN"))

#' @export
#'
setMethod("getSNPs", "MethylationExperiment", 
          function(object, ...) assay(altExp(object, "SNPs"), "Beta"))

#' @export
#'
setMethod("getSNPs", "MethBlockExperiment", 
          function(object, ...) metadata(object)$SNPs)


# used in missingness() below
setClassUnion("num_or_char", c("numeric", "character"))

if (!isGeneric("missingness")) {
  setGeneric("missingness", 
             function(X, MARGIN, i, ...) standardGeneric("missingness"))
}

#' @export
#'
setMethod("missingness", c(X="ANY", MARGIN="missing", i="missing"),
          function(X, MARGIN, i, ...) missingness(X, MARGIN=2, i=1))

#' @export
#'
setMethod("missingness", c(X="ANY", MARGIN="numeric", i="missing"),
          function(X, MARGIN, i, ...) missingness(X, MARGIN=MARGIN, i=1))

#' @export
#'
setMethod("missingness", c(X="ANY", MARGIN="missing", i="num_or_char"),
          function(X, MARGIN, i, ...) missingness(X, MARGIN=2, i=i))

#' @export
#'
setMethod("missingness", c("SummarizedExperiment", "numeric", "num_or_char"),
  function(X, MARGIN, i, ...) {

    if (0L == dim(X)[MARGIN])
      stop("Nothing to evaluate: dim(<",class(X),">)[",MARGIN,"] == 0")
    if (0L == length(assays(X)))
      stop("Nothing to evaluate: length(assays(<",class(X),">)) == 0")
    if (is.numeric(i) & length(assays(X)) < i)
      stop("Nothing to evaluate: length(assays(<",class(X),">)) < ", i)
    if (is.character(i) & !(i %in% assayNames(X)))
      stop("Nothing to evaluate: ",i," not in assayNames(<",class(X),">)")
    if (!is(MARGIN, "integer") | !(MARGIN %in% c(1L, 2L)))
      stop("MARGIN must be an integer (1 or 2).")

    switch(MARGIN,
           rowSums(is.na(assay(X, i))) / ncol(X),
           colSums(is.na(assay(X, i))) / nrow(X))
  
  })

#' @export
#'
setMethod("missingness", c("matrix", "numeric", "num_or_char"),
  function(X, MARGIN, i, ...) {
    
    if (0L == dim(X)[MARGIN]) 
      stop("Nothing to evaluate: dim(<", class(X), ">)[", MARGIN, "] == 0")
    if (!is(MARGIN, "integer") | !(MARGIN %in% c(1L, 2L)))
      stop("MARGIN must be an integer (1 or 2).")
    
    switch(MARGIN,
           rowSums(is.na(X)) / ncol(X),
           colSums(is.na(X)) / nrow(X))

  })
