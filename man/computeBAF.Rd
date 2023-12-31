% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/logR_BAF_purity_ploidy.R
\name{computeBAF}
\alias{computeBAF}
\title{computeBAF}
\usage{
computeBAF(nMajor, nMinor, purity, digits = 4)
}
\arguments{
\item{nMajor}{the number of copies of the major allele}

\item{nMinor}{the number of copies of the minor allele}

\item{purity}{the purity estimate of the tumour}

\item{digits}{a numeric, the number of digits to round to (default: 4)}
}
\value{
A vector of two numbers representing the BAF values
}
\description{
Compute the theoretical BAF values for a given segment
}
\details{
This function computes the theoretical BAF values for a given segment (from nMajor, nMinor and purity values).
}
\examples{
# A 2+1 state in a tumour with 90\% purity
computeBAF(2, 1, 0.9)
# A 1+0 state in a tumour with 60\% purity
computeBAF(1, 0, 0.6)
}
\seealso{
https://doi.org/10.1038/s41592-020-01013-2
}
\author{
tlesluyes
}
