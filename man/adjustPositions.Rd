% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adjustPositions.R
\name{adjustPositions}
\alias{adjustPositions}
\title{adjustPositions}
\usage{
adjustPositions(
  DF,
  CHRsize,
  chr_column = "chr",
  start_column = "start",
  end_column = "end",
  suffix = "_adj"
)
}
\arguments{
\item{DF}{a data.frame}

\item{CHRsize}{a data.frame from the \code{load_CHRsize} function}

\item{chr_column}{a column name with chromosome information (default: "chr")}

\item{start_column}{a column name with start position (default: "start")}

\item{end_column}{a column name with end position (default: "end")}

\item{suffix}{a suffix for the adjusted positions (default: "_adj")}
}
\value{
A data.frame with adjusted genomic positions
}
\description{
Adjust genomic positions
}
\details{
This function adjusts genomic positions according to the chromosome sizes. The first nucleotide of chromosome 2 corresponds to the size of the chromosome 1 + 1bp and so on.
}
\examples{
DF=data.frame(chr=c(1:3), start=rep(1e6, 3), end=rep(125e6, 3))
load_CHRsize("hg19")
adjustPositions(DF, CHRsize)
}
\author{
tlesluyes
}
