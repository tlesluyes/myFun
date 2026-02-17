#' @title occurrenceGRanges
#' @description Get the occurrence of events
#' @details This function gets the occurrence of events in a list of GRanges objects. All objects must have the same metadata columns and metadata must be TRUE/FALSE.
#' @param myGRList a list of GRanges objects, each object should correspond to one CNA profile
#' @param myMetadata a vector of metadata to consider
#' @param computeP whether to compute the percentage of samples with specific metadata for each region (default: TRUE)
#' @return A GRanges object with nSamples as the total number of samples and metadata columns with the occurrence of events
#' @examples
#' GR1 <- GenomicRanges::GRanges(seqnames="1",
#'                               ranges=IRanges::IRanges(start=1, end=1000),
#'                               Gain=TRUE, Loss=FALSE)
#' GR2 <- GenomicRanges::GRanges(seqnames="1",
#'                               ranges=IRanges::IRanges(start=10, end=2000),
#'                               Gain=FALSE, Loss=TRUE)
#' occurrenceGRanges(list(GR1, GR2), c("Gain", "Loss"))
#' @author tlesluyes
#' @export
occurrenceGRanges <- function(myGRList, myMetadata, computeP=TRUE) {
  checkGRlist(myGRList)
  stopifnot(typeof(myMetadata)=="character")
  stopifnot(length(computeP) == 1 && is.logical(computeP))
  stopifnot(all(sapply(myGRList, function(x) all(myMetadata %in% names(mcols(x)))))) # Make sure metadata exists
  OUT <- disjoin(Reduce(c, myGRList)) # Create a list of all regions
  OUT$nSamples <- countOverlaps(OUT, GRangesList(myGRList)) # Get the number of samples for each region
  for (i in myMetadata) { # For each metadata
    stopifnot(all(sapply(myGRList, function(x) all(mcols(x)[, i] %in% c(TRUE, FALSE))))) # Make sure metadata is TRUE/FALSE
    mcols(OUT)[, i] <- countOverlaps(OUT, GRangesList(lapply(myGRList, function(x) x[mcols(x)[, i], ]))) # Get the number of samples with specific metadata for each region
    if (isTRUE(computeP)) {
      stopifnot(! paste0(i, "_p") %in% names(mcols(OUT))) # Make sure percentage column does not already exist
      mcols(OUT)[, paste0(i, "_p")] <- round(mcols(OUT)[, i] / OUT$nSamples * 100, 4) # Compute percentage
    }
  }; rm(i)
  return(OUT)
}
