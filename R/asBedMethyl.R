#' Dump DNA methylation values into a .bedmethyl file for e.g. MARLIN or Lamprey
#'
#' @param x         a SummarizedExperiment-like object (DO NOT dump huge ones)
#' @param BPPARAM   a BiocParallelParam (default is SerialParam(progressbar=T))
#' 
#' @return list     a list of filenames dumped
#'
#' @details         This function will squawk if x has more than 50 columns. 
#'                  If x is a methBlockExperiment, fromMethBlocks is used to
#'                  ensure that locus-level representations exist in output.
#'                  The output is neither compressed nor tabixed at present.
#'                  bedMethyl column spec (with notes on R implementation):
#'                  | ---------- | --------- | ------------ |
#'                  | *name*     | *type*    | *notes*      |
#'                  | ---------- | --------- | ------------ |
#'                  | chrom      | character |              |
#'                  | ---------- | --------- | ------------ |
#'                  | start      | integer   |              |
#'                  | ---------- | --------- | ------------ |
#'                  | end        | integer   |              |
#'                  | ---------- | --------- | ------------ |
#'                  | name       | character | m or h       |
#'                  | ---------- | --------- | ------------ |
#'                  | score      | integer   | 0-1000, pval |
#'                  | ---------- | --------- | ------------ |
#'                  | strand     | character | +/-/.        |
#'                  | ---------- | --------- | ------------ |
#'                  | thickStart | integer   |              |
#'                  | ---------- | --------- | ------------ |
#'                  | thickEnd   | integer   |              |
#'                  | ---------- | --------- | ------------ |
#'                  | color      | character | '255,0,0'    |
#'                  | ---------- | --------- | ------------ |
#'                  | nValidCov  | character | set to 100   |
#'                  | ---------- | --------- | ------------ |
#'                  | percMod    | numeric   | 0-100%       |
#'                  | ---------- | --------- | ------------ |
#'                  | nMod       | integer   | methylated   |
#'                  | ---------- | --------- | ------------ |
#'                  | nCanon     | integer   | unmethlyated |
#'                  | ---------- | --------- | ------------ |
#'                  | nOther     | integer   | other calls  |
#'                  | ---------- | --------- | ------------ |
#'                  | nDelete    | integer   | reads w/dels |
#'                  | ---------- | --------- | ------------ |
#'                  | nFail      | integer   | failed calls |
#'                  | ---------- | --------- | ------------ |
#'                  | nDiff      | integer   | variant base |
#'                  | ---------- | --------- | ------------ |
#'                  | nNoCall    | integer   | ref nocalls  |
#'                  | ---------- | --------- | ------------ |
#'                  For more, consult the [bedMethyl description](https://genome.ucsc.edu/goldenpath/help/bedMethyl.html).
#'
#' @import BiocParallel
#' @import Seqinfo
#'
#' @export
asBedMethyl <- function(x, BPPARAM=SerialParam(progressbar=TRUE)) {

  x <- sort(sortSeqlevels(x)) ## must be sorted in order to tabix later
  g <- unique(genome(x))
  stopifnot(!is.null(g))
  if (class(BPPARAM) == "SerialParam" & ncol(x) > 50) { 
    warning("You are dumping more than 50 samples at a time, serially.") 
    warning("Either use a MulticoreParam object or dump fewer samples.") 
  }
  
  # if x is a methBlockExperiment, expand using fromMethBlocks()
  if (is(x, "methBlockExperiment")) {
    stop("Need to automate fromMethBlocks() using data(probes)")
    message("Expanding x using fromMethBlocks()...")
    data("probes", package="methBlockR")
    loci <- as(mcols(probes)[[g]], "GRanges") 
    x <- fromMethBlocks(x, loci)
  }

  # tidy this up ONCE
  gr <- granges(x)
  gr$score <- 0L
  gr$name <- "m" # always CpG 
  mcols(gr) <- mcols(gr)[, c("score", "name")]
  gr$thickStart <- start(gr)
  gr$thickEnd <- end(gr) 
  gr$color <- "255,0,0"    # per UCSC demo 
  gr$nValidCov <- 100L     # beads == reads for us
  gr$percMod <- NA_real_   # will swap beta values
  gr$nMod <- NA_integer_   # round(nValidCov * percMod)
  gr$nCanon <- NA_integer_ # nValidCov - nMod
  gr$nOther <- 0L          # can't detect on array
  gr$nDelete <- 0L         # can't detect on array
  gr$nFail <- 0L           # can't detect on array
  gr$nDiff <- 0L           # can't detect on array
  gr$nNoCall <- 0L         # can't detect on array

  targets <- colnames(x)
  names(targets) <- colnames(x)
  output <- bplapply(targets, .dumpBedMethyl, x=x, gr=gr, BPPARAM=BPPARAM)
  invisible(output)

}


# add coordinates for features 
.dumpBedMethyl <- function(y, x, gr) {

  what <- "CpG"
  fmt <- "bedMethyl"
  yy <- gsub(" ", "", y)
  g <- unique(genome(gr))
  fname <- paste(yy, what, g, fmt, sep=".")
  header <- paste0('track type="', fmt, '" ', 
                   'name="', y, '" ', 
                   'description="', g, " ", what, " methylation for ", y, '" ',
                   'visibility="pack"')
  columns <- c("seqnames", "start", "end", "name", "score", "strand", # BED6
               "thickStart", "thickEnd", "color",                     # BED9
               "nValidCov",                                           # preset
               "percMod", "nMod", "nCanon",                           # varies
               "nOther", "nDelete", "nFail", "nDiff", "nNoCall")      # preset
  stopifnot(all(columns %in% names(as(gr[1], "data.frame"))))

  gr$percMod <- round(getBeta(x[, y])[names(gr),] * gr$nValidCov)
  gr$nMod <- as.integer(gr$percMod)
  gr$nCanon <- as.integer(gr$nValidCov - gr$nMod)
  contents <- as(subset(gr, !is.na(gr$percMod)), "data.frame")[, columns]

  message("Dumping ", outfile, "...")
  cat(header, "\n", file=fname)
  write.table(contents,
              file=fname, append=T, row.names=F, col.names=F, quote=F, sep="\t")
  return(fname)
}
