% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotME.R
\name{plotME}
\alias{plotME}
\title{plotME}
\usage{
plotME(
  x,
  k = 100,
  minmean = 0.2,
  maxmean = 0.8,
  maxNA = 0.2,
  splitBy = NULL,
  rankBy = c("extremality", "sd"),
  pal = c("prp", "jet"),
  ...,
  BPPARAM = NULL
)
}
\arguments{
\item{x}{a MethylationExperiment or MethBlockExperiment}

\item{k}{maximum number of features to use (100)}

\item{minmean}{minimum mean for a row to be included (0.2)}

\item{maxmean}{maximum mean for a row to be included (0.8)}

\item{maxNA}{maximum fraction NA for a feature (0.2)}

\item{splitBy}{a column name on which to split rows (NULL)}

\item{rankBy}{'extremality' or 'sd' (extremality)}

\item{pal}{palette for SNP plot ("prp" or "jet") ("prp")}

\item{...}{parameters to pass to Heatmap}

\item{BPPARAM}{BiocParallelParam() to pass to plotSNPcalls (SerialParam())}
}
\description{
plot a MethylationExperiment or MethBlockExperiment features + SNPs
}
\details{
If SNPs cannot be found, only the main methylation plot is returned.
For reasons not yet clear, it's usually best to leave BPPARAM alone.
For reasons quite clear, caching SNP calls is a very good idea, via
metadata(x)$SNPcalls <- SNPcalls(x), and will speed up things A LOT
}
