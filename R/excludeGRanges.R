#' @title excludeGRanges
#' @description Exclude GRanges regions
#' @details This function excludes GRanges regions from a reference GRanges object.
#' @param GR.ref a GRanges object to be filtered
#' @param GR.toremove a GRanges object containing regions to exclude
#' @return A filtered GRanges object
#' @examples
#' GR1 <- GenomicRanges::GRanges(seqnames = rep(1, 2),
#'                               ranges=IRanges::IRanges(start=c(1, 5e5+1), end=c(5e5, 1e6)),
#'                               score=c(3, 2))
#' GR2 <- GenomicRanges::GRanges(seqnames = 1,
#'                               ranges=IRanges::IRanges(start=4e5+1, end=6e5))
#' excludeGRanges(GR1, GR2)
#' @author tlesluyes
#' @export
excludeGRanges <- function(GR.ref, GR.toremove) {
  stopifnot(inherits(GR.ref, "GRanges") && inherits(GR.toremove, "GRanges"))
  GR.ref.tmp <- disjoin(c(GR.ref, GR.toremove))
  hits <- findOverlaps(GR.ref.tmp, GR.toremove)
  GR.ref.tmp <- GR.ref.tmp[-unique(queryHits(hits))]
  GR.ref <- pintersect(findOverlapPairs(GR.ref, GR.ref.tmp))
  GR.ref$hit <- NULL
  return(GR.ref)
}
