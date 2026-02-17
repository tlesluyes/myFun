#' @title checkGR
#' @description Check that the given object is a valid GRanges object
#' @details This function checks that the given object is a valid GRanges object.
#' @param myGR a GRanges object
#' @param checkOverlaps whether to check for overlapping ranges in the GRanges object (default: TRUE)
#' @return TRUE if the input is a valid GRanges object
#' @examples
#' GR1 <- GenomicRanges::GRanges(seqnames="1",
#'                               ranges=IRanges::IRanges(start=c(1, 1001), end=c(1000, 2000)))
#' checkGR(GR1)
#' GR2 <- GenomicRanges::GRanges(seqnames="1",
#'                               ranges=IRanges::IRanges(start=c(1, 500), end=c(1000, 1000)))
#' checkGR(GR2, checkOverlaps=FALSE)
#' @author tlesluyes
#' @export
checkGR <- function(myGR, checkOverlaps=TRUE) {
  stopifnot(inherits(myGR, "GRanges")) # Entry must be a GRanges object
  if(checkOverlaps) {
    stopifnot(all(countOverlaps(myGR, myGR)==1)) # GRanges object must not have overlapping ranges
  }
  invisible(TRUE)
}

#' @title checkGRlist
#' @description Check that the given object is a list of valid GRanges objects
#' @details This function checks that the given object is a list of valid GRanges objects.
#' @param myGRList a list of GRanges objects
#' @param checkOverlaps whether to check for overlapping ranges in the GRanges objects (default: TRUE)
#' @return TRUE if the input is a list of valid GRanges objects
#' @examples
#' GR1 <- GenomicRanges::GRanges(seqnames="1", ranges=IRanges::IRanges(start=1, end=1000))
#' GR2 <- GenomicRanges::GRanges(seqnames="1", ranges=IRanges::IRanges(start=10, end=2000))
#' checkGRlist(list(GR1, GR2))
#' @author tlesluyes
#' @export
checkGRlist <- function(myGRList, checkOverlaps=TRUE) {
  stopifnot(typeof(myGRList)=="list") # myGRList must be a list
  stopifnot(length(myGRList)>1) # myGRList must have several entries
  stopifnot(all(sapply(myGRList, checkGR, checkOverlaps=checkOverlaps))) # All entries must be valid GRanges objects
  invisible(TRUE)
}
