#' @importClassesFrom RaggedExperiment RaggedExperiment
setClassUnion("CNorNULL", c("GRangesList", "RaggedExperiment", "NULL"))

setClassUnion("matrixorNULL", c("matrix", "NULL"))

setClass("MethylationExperiment", 
                contains="SingleCellExperiment",
                slots=c(CN="CNorNULL", SNPs="matrixorNULL"))

setClass("MethBlockExperiment", contains="MethylationExperiment")
