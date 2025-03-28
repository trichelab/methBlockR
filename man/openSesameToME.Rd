% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MethylationExperiment.R
\name{openSesameToME}
\alias{openSesameToME}
\title{Construct a MethylationExperiment from IDATs.}
\usage{
openSesameToME(IDATs, BPPARAM = NULL, intensities = TRUE, cd = NULL, ...)
}
\arguments{
\item{IDATs}{a vector of IDAT stubs (Basenames) with paths to them}

\item{BPPARAM}{a BiocParallel::MulticoreParam() or similar for processing}

\item{intensities}{process total intensities into assay(x, "CN")? (TRUE)}

\item{cd}{optional colData; rownames must match colnames of betas}

\item{...}{other arguments passed to sesame::openSesame()}
}
\description{
Construct a MethylationExperiment from IDATs.
}
\details{
SNPs will be read in first, so that any failures happen quickly.
BPPARAM=MulticoreParam(progressbar=TRUE) is a good starting point.
We hope to eventually support FileSetExperiment as a backend,
or some sort of materialized slice-able on-disk format.
}
\seealso{
BiocParallel::MulticoreParam

MethylationExperiment

sesame::openSesame
}
