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
setClassUnion("mat_or_df", c("matrix", "data.frame"))

# used in missingness() below
.missingnessParamsOk <- function(X, MARGIN, i=1L, ...) { 

  msg <- NULL
  if (!(MARGIN %in% 1:2)) {
    msg <- c(msg, "MARGIN must be 1 or 2 for missingness(X, MARGIN)")
  } else if (0L == dim(X)[MARGIN]) {
    msg <- c(msg, "Nothing to do: dim(<",class(X),">)[",MARGIN,"] == 0")
  }

  if (is(X, "SummarizedExperiment")) {
    if (0L == length(assays(X)))
      msg <- c(msg, "Nothing to do: length(assays(<", class(X), ">)) == 0")
    if (is.numeric(i) & length(assays(X)) < i)
      msg <- c(msg, "Nothing to do: length(assays(<", class(X), ">)) < ", i)
    if (is.character(i) & !(i %in% assayNames(X)))
      msg <- c(msg,"Nothing to do: ", i, " %in% assayNames(<", class(X), ">)")
  }

  if (!is.null(msg)) {
    message("Error: ", msg) 
    return(FALSE)
  } else {
    return(TRUE)
  }

}

#' missingness() 
#'
#' What fraction of a rectangular data structure is NA, by row or by column?
#'
#' @param X       a rectangular data structure
#' @param MARGIN  which margin to tally on (1=rows, 2=columns) (2) 
#' @param i       which assay to tally, if X is a SummarizedExperiment
#' @param ...     other parameters, currently ignored
#'
#' @examples  
#'
#' vals <- runif(n=100)
#' is.na(vals) <- as.logical(rbinom(length(vals), size=1, prob=0.1))
#' mat <- matrix(vals, ncol=20, nrow=5)
#' colnames(mat) <- paste0("column", seq_len(ncol(mat)))
#' rownames(mat) <- paste0("row", seq_len(nrow(mat)))
#' missingness(mat)
#' missingness(mat, 1)
#'
#' @seealso applySCE
#' @seealso rowSums
#' @seealso colSums
#' @seealso is.na
#'
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
setMethod("missingness", c("mat_or_df", "numeric", "num_or_char"),
  function(X, MARGIN, i, ...) {
    
    if (.missingnessParamsOk(X, MARGIN=MARGIN)) { 
      switch(MARGIN, rowSums(is.na(X))/ncol(X), colSums(is.na(X))/nrow(X))
    }

  })

#' @export
#'
setMethod("missingness", c("SummarizedExperiment", "numeric", "num_or_char"),
  function(X, MARGIN, i, ...) {

    if (.missingnessParamsOk(X, MARGIN, i)) {
      missingness(assay(X, i), MARGIN=MARGIN)
    }

  })


#' @export
#' @importMethodsFrom SummarizedExperiment show
#' 
setMethod("show", "MethBlockExperiment", 
  function(object) {
    callNextMethod()
    cat("genome:", unique(genome(object)), "\n")
  })
