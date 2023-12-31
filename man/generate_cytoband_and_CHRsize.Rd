% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_cytoband_and_CHRsize.R
\name{generate_cytoband_and_CHRsize}
\alias{generate_cytoband_and_CHRsize}
\title{generate_cytoband_and_CHRsize}
\usage{
generate_cytoband_and_CHRsize(cytoband_file)
}
\arguments{
\item{cytoband_file}{a cytoband file}
}
\value{
A list with both the \code{cytoband} and \code{CHRsize} information
}
\description{
Generate \code{cytoband} and \code{CHRsize} information
}
\details{
This function generates \code{cytoband} and \code{CHRsize} information from a cytoband file. This can be obtained from the UCSC table browser -> select a genome/assembly -> "Mapping and Sequencing" -> "Chromosome Band" (not the ideogram version!) -> "get output" -> Remove the first "#" character (keep the header!).
}
\seealso{
\verb{load_CHRsize("hg38"); load_cytoband("hg38")}
}
\author{
tlesluyes
}
