#' Given an object with methylation block coordinates, switch the genome
#'
#' Now keeping track of genome-to-genome mappings in rowData() of the result
#'
#' @param x         a GenomicRatioSet or something like it with methBlock coords
#' @param to        the new genome to use (hg19 or hg38) 
#'
#' @return          an object very much like x, but mapped onto `to`
#' 
#' @seealso         asMethBlocks
#'
#' @import          GenomicRanges 
#'
#' @export
switchMethBlocksGenome <- function(x, to=c("hg38", "hg19")) {

  stopifnot(is(x, "SummarizedExperiment"))
  
  to <- match.arg(to)
  from <- unique(genome(x))
  if (is.na(from)) stop("The genome of the object cannot be NA!") 
  if (from == to) {
    message("Genome is already ", to, ". No need to switch.")
    return(x) 
  } else {
    message("Switching genome from ", from, " to ", to, "...") 
    data("mbHg19Hg38")
    mb <- subset(mbHg19Hg38, mbHg19Hg38[[from]] %in% rownames(x))
    mappable <- nrow(mb)
    mb <- subset(mb, !is.na(mb[[to]]))
    mapped <- nrow(mb)
    message("Found ", mapped, " / ", mappable, 
            " blocks mapping from ", from, " to ", to, ".")
    x <- x[mb[[from]], ]
    rowRanges(x) <- as(mb[[to]], "GRanges")
    rowData(x)[[from]] <- mb[[from]]
    rowData(x)[[to]] <- mb[[to]]
    rownames(x) <- mb[[to]]
    genome(rowRanges(x)) <- to
    return(x)
  } 

}


# helper used to construct mbHg19Hg38
.mbMapByChr <- function(x, chr, from=c("hg19","hg38"), BPPARAM=NULL) {
  from <- match.arg(from)
  to <- setdiff(c("hg19","hg38"), from)
  x <- subset(x, !is.na(x[[to]]) & grepl(paste0(chr, ":"), x[[to]]))
  if (is.null(BPPARAM)) BPPARAM <- MulticoreParam(progressbar=TRUE) 
  for (v in c(from, to)) x[[v]] <- as.character(x[[v]])
  splt <- split(x[, to, drop=FALSE], x[[from]])
  uniqMap <- names(which(sapply(splt, nrow) == 1))
  message("Found ", length(uniqMap), " unique mappings ",
          "from ", from, " ", chr, " to ", to, ".")
  uniqueMapped <- sapply(splt[uniqMap], function(x) x[[to]])
  multiMap <- setdiff(names(splt), uniqMap)
  message("Found ", length(multiMap), " split mappings ",
          "from ", from, " ", chr, " to ", to, ".")
  message("Breaking ties by choosing the most common mapping.")
  multiMapped <- do.call(c, 
                         bplapply(splt[multiMap], 
                                  function(x) {
                                    tbl <- table(x[[to]])
                                    names(tbl)[which.max(tbl)]
                                  },
                                  BPPARAM=BPPARAM))
  res <- c(uniqueMapped, multiMapped)[x[[from]]]
  df <- data.frame(from=names(res), to=res)
  names(df) <- c(from, to)
  return(df)
} 
# system.time(chr1 <- .mbMapByChr(methBlox, "chr1", "hg19")) # about 10 minutes


# helper function
sharedMethBlockPairs <- function(mb, g=c("hg19","hg38")) {

  g <- match.arg(g)
  g0 <- "hg38"
  mb <- subset(mb, !is.na(mb[[g]]) & !is.na(mb[[g0]]))
  mb$inBlocks <- mb[["inHg19mb"]] & mb[["inHg38mb"]]
  mb$pairing <- paste0(as.character(mb[[g]]), "|", as.character(mb[[g0]]))
  dupes <- mb$pairing[which(duplicated(mb$pairing))]
  mb$dupe <- mb$pairing %in% dupes
  mb$hits <- 1
  mb$hits <- rowsum(mb$hits, mb$pairing)[mb$pairing, 1]
  mb <- subset(mb, !duplicated(pairing))[, c(g, g0, "inBlocks", "hits")]
  rownames(mb) <- NULL

  message("Pairing ", g, " and ", g0, " methylation block ranges...")
  paired <- Pairs(as(mb[[g]], "GRanges"), as(mb[[g0]], "GRanges"))
  names(paired) <- as.character(mb[[g]])
  mcols(paired)$inBlocks <- mb$inBlocks
  dupes <- names(paired)[duplicated(names(paired))]
  mcols(paired)$splitMap <- names(paired) %in% dupes
  ndups <- sum(names(paired) %in% dupes)
  mcols(paired)$hits <- mb$hits
  nuniq <- length(paired) - ndups
  message("Found ", ndups, " split mappings and ", 
          nuniq, " unique mappings from ", g, " to ", g0, ".")
  return(paired)

}

# helper function
splitSharedMethBlockPairs <- function(mb, g=c("hg19","hg38")) {

  paired <- sharedMethBlockPairs(mb, g)
  splt <- subset(paired, mcols(paired)$splitMap)
  ssplt <- split(as.character(second(splt)), as.character(first(splt)))

}
