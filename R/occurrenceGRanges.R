#' @title occurrenceGRanges
#' @description Get the occurrence of events
#' @details This function gets the occurrence of events in a list of GRanges objects. All objects must have the same metadata columns and metadata must be TRUE/FALSE.
#' @param myGRList a list of GRanges objects, each object should correspond to one CNA profile
#' @param myMetadata a vector of metadata to consider
#' @return A GRanges object with nSamples as the total number of samples and metadata columns with the occurrence of events
#' @examples
#' require("GenomicRanges")
#' GR1=GRanges(seqnames="1", ranges=IRanges(start=1, end=1000), Gain=TRUE, Loss=FALSE)
#' GR2=GRanges(seqnames="1", ranges=IRanges(start=10, end=2000), Gain=FALSE, Loss=TRUE)
#' occurrenceGRanges(list(GR1, GR2),c("Gain", "Loss"))
#' @author tlesluyes
#' @export
occurrenceGRanges=function(myGRList, myMetadata) {
  checkGRlist(myGRList)
  stopifnot(all(sapply(myGRList, function(x) all(myMetadata %in% names(GenomicRanges::mcols(x)))))) # Make sure metadata exists
  OUT=GenomicRanges::disjoin(Reduce(c, myGRList)) # Create a list of all regions
  OUT$nSamples=GenomicRanges::countOverlaps(OUT, GenomicRanges::GRangesList(myGRList)) # Get the number of samples for each region
  for (i in myMetadata) { # For each metadata
    stopifnot(all(sapply(myGRList, function(x) all(GenomicRanges::mcols(x)[, i] %in% c(TRUE, FALSE))))) # Make sure metadata is TRUE/FALSE
    GenomicRanges::mcols(OUT)[, i]=GenomicRanges::countOverlaps(OUT, GenomicRanges::GRangesList(lapply(myGRList, function(x) x[GenomicRanges::mcols(x)[, i], ]))) # Get the number of samples with specific metadata for each region
  }; rm(i)
  return(OUT)
}