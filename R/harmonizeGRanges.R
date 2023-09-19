#' @title harmonizeGRanges
#' @description Harmonize GRanges objects
#' @details This function harmonizes GRanges objects by keeping only regions covered by all samples.
#' @param myGRList a list of GRanges objects, each object should correspond to one CNA profile
#' @return A list of harmonized GRanges objects
#' @examples
#' require("GenomicRanges")
#' GR1=GRanges(seqnames="1", ranges=IRanges(start=1, end=1000), nMajor=1, nMinor=1)
#' GR2=GRanges(seqnames="1", ranges=IRanges(start=10, end=2000), nMajor=2, nMinor=1)
#' harmonizeGRanges(list(GR1, GR2))
#' @author tlesluyes
#' @export
harmonizeGRanges=function(myGRList) {
  checkGRlist(myGRList)
  ALL_REGIONS=GenomicRanges::disjoin(Reduce(c, myGRList)) # Create a list of all regions
  ALL_REGIONS=ALL_REGIONS[apply(do.call(cbind, lapply(myGRList, function(x) GenomicRanges::countOverlaps(ALL_REGIONS, x)==1)), 1, all)] # Only keep regions covered by all samples (+remove intra-sample overlaps)
  myGRList=lapply(myGRList, function(x) GenomicRanges::pintersect(IRanges::findOverlapPairs(x, ALL_REGIONS))) # Get GRanges information for those regions
  myGRList=lapply(myGRList, function(x) {names(x)=1:length(x); return(x)}) # Reset names as some regions could have been discarded
  stopifnot(sapply(2:length(myGRList), function(x) identical(GenomicRanges::granges(myGRList[[1]]), GenomicRanges::granges(myGRList[[x]])))) # Make sure regions are strictly identical
  return(myGRList)
}