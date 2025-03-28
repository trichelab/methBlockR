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
          function(object, ...) assay(altExp(object, "SNP"), "Beta"))

#' @export
#'
setMethod("getSNPs", "MethBlockExperiment", 
          function(object, ...) metadata(object)$SNPs)


if (!isGeneric("missingness")) {
  setGeneric("missingness", 
             function(X, MARGIN, i, ...) standardGeneric("missingness"))
}

# used in missingness() below
setClassUnion("num_or_char", c("numeric", "character"))

# used in missingness() below
.missingnessParamsOk <- function(X, MARGIN, i, ...) { 

  msg <- NULL
  if (!(MARGIN %in% 1:2)) 
    msg <- c(msg, "MARGIN must be either 1 or 2 for missingness()")
  if (0L == dim(X)[MARGIN])
    msg <- c(msg, "Nothing to evaluate: dim(<",class(X),">)[",MARGIN,"] == 0")
  if (is(X, "SummarizedExperiment") & (0L == length(assays(X))))
    msg <- c(msg, "Nothing to evaluate: length(assays(<",class(X),">)) == 0")
  if (is.numeric(i) & is(X, "SummarizedExperiment") & (length(assays(X)) < i))
    msg <- c(msg, "Nothing to evaluate: length(assays(<",class(X),">)) < ", i)
  if (is.character(i) & is(X, "SummarizedExperiment") & !(i %in% assayNames(X)))
    msg <- c(msg,"Nothing to evaluate: ",i," not in assayNames(<",class(X),">)")
  if (!is.null(msg)) 
    stop(msg)
  if (is.null(msg))
    TRUE

}

#' @export
#'
setMethod("missingness", c(X="ANY", MARGIN="missing", i="missing"),
          function(X, MARGIN, i, ...) missingness(X, MARGIN=2L, i=1L))

#' @export
#'
setMethod("missingness", c(X="ANY", MARGIN="numeric", i="missing"),
          function(X, MARGIN, i, ...) missingness(X, MARGIN=MARGIN, i=1L))

#' @export
#'
setMethod("missingness", c(X="ANY", MARGIN="missing", i="num_or_char"),
          function(X, MARGIN, i, ...) missingness(X, MARGIN=2L, i=i))

#' @export
#'
setMethod("missingness", c("SummarizedExperiment", "numeric", "num_or_char"),
  function(X, MARGIN, i, ...) {

    if (.missingnessParamsOk(X, MARGIN, i, ...)) { 
      switch(MARGIN,
             rowSums(is.na(assay(X, i))) / ncol(X),
             colSums(is.na(assay(X, i))) / nrow(X))
    }

  })

#' @export
#'
setMethod("missingness", c("matrix", "numeric", "num_or_char"),
  function(X, MARGIN, i, ...) {
    
    if (.missingnessParamsOk(X, MARGIN, i, ...)) { 
      switch(MARGIN,
             rowSums(is.na(assay(X, i))) / ncol(X),
             colSums(is.na(assay(X, i))) / nrow(X))
    }

  })


#' @export
#'
setMethod("missingness", c("MethylationExperiment", "numeric", "num_or_char"),
  function(X, MARGIN, i, ...) {

    if (.missingnessParamsOk(X, MARGIN, i, ...)) { 
      applySCE(X, missingness, MARGIN=MARGIN, i=i, ...)
    }

  })

#' @export
#'
setMethod("missingness", c("MethBlockExperiment", "numeric", "num_or_char"),
  function(X, MARGIN, i, ...) {
    
    # override MethylationExperiment behavior
    if (.missingnessParamsOk(X, MARGIN, i, ...)) { 
      switch(MARGIN,
             rowSums(is.na(assay(X, i))) / ncol(X),
             colSums(is.na(assay(X, i))) / nrow(X))
    }

  })
