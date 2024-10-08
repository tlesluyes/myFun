% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BED_metrics.R
\name{BED_metrics}
\alias{BED_metrics}
\title{BED_metrics}
\usage{
BED_metrics(BED, verbose = TRUE)
}
\arguments{
\item{BED}{a data.frame, a GRanges object or a path to a BED file}

\item{verbose}{a boolean, whether to print metrics (default: TRUE)}
}
\value{
A named list of metrics (number of chromosomes, number of regions and total size of the regions) before and after removing overlaps (\code{GenomicRanges::reduce()}). Strand information is not considered
}
\description{
Give BED metrics
}
\details{
This function provides several metrics of interest from a BED file/object.
}
\examples{
# BED format is 0-based for starts!
BED=data.frame(chr=c(c(1,1:3)), start=c(0, 1e3, 0, 1e3), end=c(2e3, 5e3, 2e3, 4e3))
BED_metrics(BED)
}
\author{
tlesluyes
}
