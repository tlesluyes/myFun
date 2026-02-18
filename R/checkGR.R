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

#' @title cleanGRlistMetadata
#' @description Remove regions with missing metadata in at least one sample from a list of GRanges objects
#' @details This function removes regions that have missing metadata in at least one sample from a list of GRanges objects.
#' @param myGRList a list of harmonized GRanges objects (all granges must be identical across samples)
#' @param metadataCols a character vector of metadata column names to check for missing values
#' @param checkOverlaps whether to check for overlapping ranges in the GRanges objects (default: TRUE)
#' @return A list of cleaned GRanges objects with the same regions but only those with complete metadata across all samples
#' @examples
#' GR1 <- GenomicRanges::GRanges(seqnames=rep("1", 3),
#'                               ranges=IRanges::IRanges(start=c(1, 1001, 10001),
#'                                                       end=c(1000, 10000, 20000)),
#'                               CNstatus=c(NA, "2+1", "1+1"))
#' GR2 <- GenomicRanges::GRanges(seqnames=rep("1", 2),
#'                               ranges=IRanges::IRanges(start=c(500, 10001),
#'                                                       end=c(10000, 25000)),
#'                               CNstatus=c("2+1", "1+1"))
#' myGRList <- harmonizeGRanges(list(GR1=GR1, GR2=GR2), keepHoles=TRUE)
#' cleanGRlistMetadata(myGRList, metadataCols="CNstatus")
#' @author tlesluyes
#' @export
cleanGRlistMetadata <- function(myGRList, metadataCols, checkOverlaps=TRUE) {
  stopifnot(checkGRlist(myGRList, checkOverlaps=checkOverlaps))
  stopifnot(length(metadataCols)>0 && is.character(metadataCols))
  stopifnot(all(sapply(2:length(myGRList), function(x) {identical(granges(myGRList[[1]]), granges(myGRList[[x]]))}))) # Make sure regions are strictly identical
  stopifnot(all(sapply(myGRList, function(x) all(metadataCols %in% names(mcols(x))))))
  KEEP <- apply(do.call(cbind, lapply(myGRList, function(x) mcols(x)[, metadataCols])), 1, function(x) all(! is.na(x)))
  if (! all(KEEP)) {
    print(paste0("The following regions have been removed because of missing metadata in at least one sample: ", paste0(granges(myGRList[[1]])[! KEEP], collapse=", ")))
  }
  myGRList <- lapply(myGRList, function(x) x[KEEP, ])
  return(myGRList)
}
