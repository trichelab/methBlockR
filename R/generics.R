# getGreen 
if (!isGeneric("getGreen")) {
  setGeneric("getGreen", 
             function(object) standardGeneric("getGreen"))
}
setMethod("getGreen", "SummarizedExperiment", 
          function(object) assay(object, "Green"))


# getRed
if (!isGeneric("getRed")) {
  setGeneric("getRed", 
             function(object) standardGeneric("getRed"))
}
setMethod("getRed", "SummarizedExperiment", 
          function(object) assay(object, "Red"))


# annotation
if (!isGeneric("annotation")) {
  setGeneric("annotation", 
             function(object, ...) standardGeneric("annotation"))
}
setMethod("annotation", "SummarizedExperiment", 
          function(object, ...) metadata(object)$annotation)


# annotation<-
if (!isGeneric("annotation<-")) {
  setGeneric("annotation<-", 
             function(object, ..., value) standardGeneric("annotation<-"))
}
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
setMethod("preprocessMethod", "SummarizedExperiment", 
          function(object) metadata(object)$preprocessMethod)


# preprocessMethod<-
if (!isGeneric("preprocessMethod<-")) {
  setGeneric("preprocessMethod<-", 
             function(object, ..., value) standardGeneric("preprocessMethod<-"))
}
setReplaceMethod("preprocessMethod", 
                 c(object="SummarizedExperiment", value="ANY"),
                 function(object, ..., value) {
                   metadata(object)$preprocessMethod <- value
                   return(object)
                 })


# helper fns
.M2B <- function(x) (2 ** x) / (1 + (2 ** x))
.B2M <- function(p) log2(p / (1 - p))


# getBeta
if (!isGeneric("getBeta")) {
  setGeneric("getBeta", 
             function(object, ...) standardGeneric("getBeta"))
}
setMethod("getBeta", "SummarizedExperiment", 
          function(object, ...) {
            if ("Beta" %in% assayNames(object)) assay(object, "Beta")
            else if ("M" %in% assayNames(object)) .M2B(assay(object, "M"))
            else stop("No Beta or M assay found")
          })


# getM
if (!isGeneric("getM")) {
  setGeneric("getM", 
             function(object, ...) standardGeneric("getM"))
}
setMethod("getM", "SummarizedExperiment", 
          function(object, ...) {
            if ("M" %in% assayNames(object)) assay(object, "M")
            else if ("Beta" %in% assayNames(object)) .B2M(assay(object, "Beta"))
            else stop("No Beta or M assay found")
          })


# getCN
if (!isGeneric("getCN")) {
  setGeneric("getCN", 
             function(object, ...) standardGeneric("getCN"))
}
setMethod("getCN", "SummarizedExperiment", 
          function(object, ...) assay(object, "CN"))


# getMeth
if (!isGeneric("getMeth")) {
  setGeneric("getMeth", 
             function(object) standardGeneric("getMeth"))
}
setMethod("getMeth", "SummarizedExperiment", 
          function(object) getBeta(object) * (2 ** getCN(object)))


# getUnmeth
if (!isGeneric("getUnmeth")) {
  setGeneric("getUnmeth",
             function(object) standardGeneric("getUnmeth"))
}
setMethod("getUnmeth", "SummarizedExperiment", 
          function(object) (1 - getBeta(object)) * (2 ** getCN(object)))


# missingness: 0-1 value for how many NAs are present in each column (default) or row of X
if (!isGeneric("missingness")) {
  setGeneric("missingness", 
             function(X, MARGIN, i, ...) standardGeneric("missingness"))
}

# somewhat ironic default signature for 'missingness(X="ANY", MARGIN="missing", i="missing")
setMethod("missingness", c(X="ANY", MARGIN="missing", i="missing"),
          function(X, MARGIN, i, ...) missingness(X, MARGIN=2, i=1))

# more specific signature: 'missingness(X="ANY, MARGIN=int, i="missing")
setMethod("missingness", c(X="ANY", MARGIN="numeric", i="missing"),
          function(X, MARGIN, i, ...) missingness(X, MARGIN=MARGIN, i=1))

# avoid writing multiple methods for assay indexing 
setClassUnion("num_or_char", c("numeric", "character"))

# more specific signature: 'missingness(X="ANY", MARGIN="missing", i="num_or_char")
setMethod("missingness", c(X="ANY", MARGIN="missing", i="num_or_char"),
          function(X, MARGIN, i, ...) missingness(X, MARGIN=2, i=i))

# fully specified method for a SummarizedExperiment object
setMethod("missingness", c("SummarizedExperiment", "numeric", "num_or_char"),
          function(X, MARGIN, i, ...) {
            if (0L == dim(X)[MARGIN]) {
              stop(paste0("Nothing to evaluate: dim(<", class(X), ">)[", MARGIN, "] == 0"))
            }
            if (0L == length(assays(X))) {
              stop(paste0("Nothing to evaluate: length(assays(<", class(X), ">)) == 0"))
            }
            if (is.numeric(i) & length(assays(X)) < i) {
              stop(paste0("Nothing to evaluate: length(assays(<", class(X), ">)) < ", i))
            }
            if (is.character(i) & !(i %in% assayNames(X))) {
              stop(paste0("Nothing to evaluate: ", i, " not in assayNames(<", class(X), ">)"))
            }
            switch(MARGIN,
                   rowSums(is.na(assay(X, i))) / ncol(X),
                   colSums(is.na(assay(X, i))) / nrow(X))
          })

# fully specified method for any other type of matrix-like object
setMethod("missingness", c("matrix", "numeric", "num_or_char"),
          function(X, MARGIN, i, ...) {
            if (0L == dim(X)[MARGIN]) {
              stop(paste0("Nothing to evaluate: dim(<", class(X), ">)[", MARGIN, "] == 0"))
            }
            switch(MARGIN,
                   rowSums(is.na(X)) / ncol(X),
                   colSums(is.na(X)) / nrow(X))
          })
