#' @title checkGRlist
#' @description Check that the given object is a list of GRanges objects
#' @details This function checks that the given object is a list of GRanges objects.
#' @param myGRList a list of GRanges objects
#' @return TRUE if the input is a list of GRanges objects
#' @examples
#' require("GenomicRanges")
#' GR1=GRanges(seqnames="1", ranges=IRanges(start=1, end=1000))
#' GR2=GRanges(seqnames="1", ranges=IRanges(start=10, end=2000))
#' checkGRlist(list(GR1, GR2))
#' @author tlesluyes
#' @export
checkGRlist=function(myGRList) {
  stopifnot(typeof(myGRList)=="list") # myGRList must be a list
  stopifnot(length(myGRList)>1) # myGRList must have several entries
  stopifnot(all(sapply(myGRList, inherits, "GRanges"))) # All entries must be GRanges objects
  invisible(TRUE)
}