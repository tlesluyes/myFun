% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/load_cytoband.R
\name{load_cytoband}
\alias{load_cytoband}
\title{load_cytoband}
\usage{
load_cytoband(assembly)
}
\arguments{
\item{assembly}{an assembly (hg19 or hg38)}
}
\value{
A data.frame with the \code{cytoband} information
}
\description{
Load \code{cytoband} information
}
\details{
This function loads \code{cytoband} information for a given assembly. It is then available as a data.frame called \code{cytoband} in the environment.
}
\examples{
load_cytoband("hg38"); head(cytoband)
}
\author{
tlesluyes
}
