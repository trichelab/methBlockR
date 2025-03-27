# getGreen 
if (!isGeneric("getGreen")) {
  setGeneric("getGreen", 
             function(object) standardGeneric("getGreen"))
}

#' @export
#'
setMethod("getGreen", "SummarizedExperiment", 
          function(object) assay(object, "Green"))


# getRed
if (!isGeneric("getRed")) {
  setGeneric("getRed", 
             function(object) standardGeneric("getRed"))
}

#' @export
#'
setMethod("getRed", "SummarizedExperiment", 
          function(object) assay(object, "Red"))


# annotation
if (!isGeneric("annotation")) {
  setGeneric("annotation", 
             function(object, ...) standardGeneric("annotation"))
}

#' @export
#'
setMethod("annotation", "SummarizedExperiment", 
          function(object, ...) metadata(object)$annotation)


# annotation<-
if (!isGeneric("annotation<-")) {
  setGeneric("annotation<-", 
             function(object, ..., value) standardGeneric("annotation<-"))
}

#' @export
#'
setReplaceMethod("annotation", c(object="SummarizedExperiment", value="ANY"),
                 function(object, ..., value) {
                   metadata(object)$annotation <- value
                   return(object)
                 })


# preprocessMethod
if (!isGeneric("preprocessMethod")) {
  setGeneric("preprocessMethod", 
             function(object) standardGeneric("preprocessMethod"))
}

#' @export
#
setMethod("preprocessMethod", "SummarizedExperiment", 
          function(object) metadata(object)$preprocessMethod)



# helper fns
.M2B <- function(x) (2 ** x) / (1 + (2 ** x))
.B2M <- function(p) log2(p / (1 - p))

if (!isGeneric("getBeta")) {
  setGeneric("getBeta", 
             function(object, ...) standardGeneric("getBeta"))
}

#' @export
#
setMethod("getBeta", "SummarizedExperiment", 
          function(object, ...) {
            if ("Beta" %in% assayNames(object)) assay(object, "Beta")
            else if ("M" %in% assayNames(object)) .M2B(assay(object, "M"))
            else stop("No Beta or M assay found")
          })

if (!isGeneric("getM")) {
  setGeneric("getM", 
             function(object, ...) standardGeneric("getM"))
}

#' @export
#
setMethod("getM", "SummarizedExperiment", 
          function(object, ...) {
            if ("M" %in% assayNames(object)) assay(object, "M")
            else if ("Beta" %in% assayNames(object)) .B2M(assay(object, "Beta"))
            else stop("No Beta or M assay found")
          })


if (!isGeneric("getCN")) {
  setGeneric("getCN", 
             function(object, ...) standardGeneric("getCN"))
}

#' @export
#
setMethod("getCN", "SummarizedExperiment", 
          function(object, ...) assay(object, "CN"))

#' @export
#
setMethod("getCN", "MethylationExperiment", 
          function(object, ...) object@CN)


if (!isGeneric("getMeth")) {
  setGeneric("getMeth", 
             function(object) standardGeneric("getMeth"))
}

#' @export
#
setMethod("getMeth", "SummarizedExperiment", 
          function(object) getBeta(object) * (2 ** getCN(object)))


if (!isGeneric("getUnmeth")) {
  setGeneric("getUnmeth",
             function(object) standardGeneric("getUnmeth"))
}

#' @export
#
setMethod("getUnmeth", "SummarizedExperiment", 
          function(object) (1 - getBeta(object)) * (2 ** getCN(object)))


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


if (!isGeneric("getBlocks")) {
  setGeneric("getBlocks", 
             function(object, ...) standardGeneric("getBlocks"))
}

#' @export
#'
setMethod("getBlocks", "MethylationExperiment", 
          function(object, ...) metadata(object)$blocks)

#' @export
#'
setMethod("getBlocks", "MethBlockExperiment", 
          function(object, ...) object@blocks)


if (!isGeneric("getSNPs")) {
  setGeneric("getSNPs", 
             function(object, ...) standardGeneric("getSNPs"))
}

#' @export
#'
setMethod("getSNPs", "MethylationExperiment", 
          function(object, ...) object@SNPs)

#' @export
#'
setMethod("getBlocks", "MethBlockExperiment", 
          function(object, ...) object@blocks)
