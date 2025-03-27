#' MethylationExperiment() constructor function 
#'
#' A MethylationExperiment object must have its ncols() the same across 
#' assays, CN, and SNPs (if the latter are not NULL). The rownames are 
#' usually different (as are the number of rows) for these slots. 
#'
#' @param ...   arguments passed to the SingleCellExperiment() constructor
#' @param CN    a copynumber representation (see Details) or NULL (the default).
#' @param SNPs  a SNP matrix or NULL (the default).
#'
#' @return  a MethylationExperiment
#'
#' @details the @CN slot can be a GRangesList, a RaggedExperiment, or NULL. 
#'          If openSesameToME is used, a matrix of total intensities can 
#'          be added to assays(object, "CN") with per-probe intensities, but
#'          the @CN slot exists to hold CN estimates across _multiple_ probes,
#'          and is meant to be populated via conumee2::CNV.fit or similar.
#'
#' @seealso openSesameToME
#'
#' @import methods
#'
#' @importFrom  GenomicRanges         GRangesList
#' @importFrom  RaggedExperiment      RaggedExperiment
#' @importFrom  SingleCellExperiment  SingleCellExperiment
#'
#' @rdname  MethylationExperiment
#'
#' @export
#'
MethylationExperiment <- function(..., CN=NULL, SNPs=NULL) {
  
  y <- SingleCellExperiment(...)
  class(y) <- "MethylationExperiment"
  y@SNPs <- SNPs
  y@CN <- CN
  y

}


#' Construct a MethylationExperiment from IDATs.
#'
#' @param IDATs       a vector of IDAT stubs (Basenames) with paths to them
#' @param BPPARAM     a BiocParallel::MulticoreParam() or similar for processing
#' @param intensities process CN (total intensity)? (TRUE) 
#' @param cd          optional colData; rownames must match colnames of betas
#' @param ...         other arguments passed to sesame::openSesame()
#'
#' @details SNPs will be read in first, so that any failures happen quickly. 
#'          BPPARAM=MulticoreParam(progressbar=TRUE) is a good starting point.
#'          We hope to eventually support FileSetExperiment as a backend.
#'
#' @seealso BiocParallel::MulticoreParam
#' @seealso MethylationExperiment
#' @seealso sesame::openSesame
#'
#' @importFrom  sesame openSesame
#' @import BiocParallel
#'
#' @export
#'
openSesameToME <- function(IDATs, BPPARAM=NULL, intensities=TRUE, cd=NULL, ...){

  # BPPARAM can default to a SerialParam, although that's dumb
  if (is.null(BPPARAM)) BPPARAM <- SerialParam(progressbar=TRUE)

  # try to read SNPs first; will fail faster than betas or CN 
  message("Processing SNPs for sample barcoding...")
  SNPs <- openSesame(IDATs, BPPARAM=BPPARAM, func=getAFs) 
   
  # if that succeeds, read beta values
  message("Processing beta values...")
  betas <- openSesame(IDATs, ..., BPPARAM=BPPARAM, func=getBetas)

  # skeletal colData and a backup for rowname mismatches: 
  cd0 <- DataFrame(path=IDATs, row.names=colnames(betas))
  if (!is.null(cd)) {
    if (!identical(colnames(betas), rownames(cd))) { 
      message("rownames(cd) != colnames(betas). Discarding.")
      cd <- cd0
    }
  } else cd <- cd0

  # if requested, CN
  if (intensities) {
    message("Processing intensities for copy number analysis...")
    CN <- openSesame(IDATs, func=totalIntensities, BPPARAM=BPPARAM)
    MethylationExperiment(assays=list(Beta=betas, CN=CN), SNPs=SNPs, colData=cd)
  } else { 
    MethylationExperiment(assays=list(Beta=betas), SNPs=SNPs, colData=cd)
  }

}


setValidity2("MethylationExperiment", 
  function(object) {
 
    NC <- NCOL(object)
    msg <- NULL
    
    if (!is.null(getSNPs(object))) {
      if (NCOL(getSNPs(object)) != NC) {
        msg <- c(msg, "ncol(getSNPs(object)) != ncol(object)")
      }
    }

    if (!is.null(getCN(object))) { 

      if (is(getCN(object), "RaggedExperiment")) {
        if (NCOL(getCN(object)) != NC) {
          msg <- c(msg, "ncol(getCN(object)) != ncol(object)")
        }
      }

      if (is(getCN(object), "GRangesList") | is(getCN(object), "CNV.analysis")){
        if (length(getCN(object)) != NC) {
          msg <- c(msg, "length(getCN(object)) != ncol(object)")
        }
      }

    }

    if (length(msg)) {
      msg
    } else TRUE

})


#' @rdname  MethylationExperiment
#'
#' @export 
#'
setReplaceMethod("colnames", "MethylationExperiment",
  function(x, value) {

    y <- callNextMethod() # SingleCellExperiment parent method
    if (!is.null(getSNPs(x))) colnames(y@SNPs) <- value
    if (!is.null(getCN(x))) colnames(y@CN) <- value
    return(y)

  }
)


# add generic CN replacer
if (!isGeneric("CN<-")) {
  setGeneric("CN<-", function(x, value) standardGeneric("CN<-"))
}


#' @rdname  MethylationExperiment
#'
#' @export 
#'
setReplaceMethod("CN", "MethylationExperiment",
  function(x, value) {

    if (!is.null(value)) {
      if (is(value, "RaggedExperiment")) {
        stopifnot(identical(colnames(value), colnames(x)))
      } else {
        stopifnot(identical(names(value), colnames(x)))
      }
    }

    x@CN <- value
    validObject(x)
    x

  }
)


# add generic SNPs replacer
if (!isGeneric("SNPs<-")) {
  setGeneric("SNPs<-", function(x, value) standardGeneric("SNPs<-"))
}


#' @rdname  MethylationExperiment
#'
#' @export 
#'
setReplaceMethod("SNPs", "MethylationExperiment", 
  function(x, value) {

    if (!is.null(value)) stopifnot(identical(colnames(value), colnames(x)))
    x@SNPs <- value
    validObject(x)
    x

  }
)


#' @rdname  MethylationExperiment
#'
#' @export 
#'
setMethod("[", "MethylationExperiment", 
  function(x, i, j, drop=TRUE) {

    if (!missing(j)) {
      if (is.character(j)) {
        fmt <- paste0("<", class(x), ">[,j] index out of bounds: %s")
        j <- SummarizedExperiment:::.SummarizedExperiment.charbound(
          j, colnames(x), fmt
        )
      }
      j <- as.vector(j)
     
      SNPs <- getSNPs(x)
      if (!is.null(SNPs)) SNPs <- SNPs[, j, drop=FALSE]

      CN <- getCN(x)
      if (!is.null(CN)) {
        if (is(CN, "RaggedExperiment")) CN <- CN[, j, drop=FALSE]
        else CN <- CN[j] # GRL, CNV.analysis, etc.
      }
    }

    out <- callNextMethod()
    BiocGenerics:::replaceSlots(out, SNPs=SNPs, CN=CN, check=FALSE)

  }
)
