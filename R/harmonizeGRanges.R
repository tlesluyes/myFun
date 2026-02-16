#' @title harmonizeGRanges
#' @description Harmonize GRanges objects
#' @details This function harmonizes GRanges objects by keeping only regions covered by all samples.
#' @param myGRList a list of GRanges objects, each object should correspond to one CNA profile
#' @param cores a numeric, the number of cores to use (default: 1)
#' @param keepHoles a logical, whether to keep holes (regions not covered by all samples) in the output (default: FALSE)
#' @return A list of harmonized GRanges objects
#' @examples
#' GR1 <- GenomicRanges::GRanges(seqnames="1",
#'                               ranges=IRanges::IRanges(start=1, end=1000),
#'                               nMajor=1, nMinor=1)
#' GR2 <- GenomicRanges::GRanges(seqnames="1",
#'                               ranges=IRanges::IRanges(start=10, end=2000),
#'                               nMajor=2, nMinor=1)
#' harmonizeGRanges(list(GR1, GR2))
#' @author tlesluyes
#' @export
harmonizeGRanges <- function(myGRList, cores=1, keepHoles=FALSE) {
  stopifnot(length(cores)==1 && is.numeric(cores) && cores>=1)
  stopifnot(length(keepHoles)==1 && is.logical(keepHoles))
  if (cores>1) registerDoParallel(cores=cores)
  checkGRlist(myGRList)
  ALL_REGIONS <- disjoin(Reduce(c, myGRList)) # Create a list of all regions
  if (isTRUE(keepHoles)) {
    myGRList <- foreach(x=myGRList, .final=function(x) setNames(x, names(myGRList))) %dopar% {return(sort(pintersect(findOverlapPairs(c(x, setdiff(ALL_REGIONS, x)), ALL_REGIONS))))}
  } else {
    ALL_REGIONS <- ALL_REGIONS[apply(foreach(x=myGRList, .combine=cbind) %dopar% {countOverlaps(ALL_REGIONS, x)==1}, 1, all)] # Only keep regions covered by all samples (+remove intra-sample overlaps)
    myGRList <- foreach(x=myGRList, .final=function(x) setNames(x, names(myGRList))) %dopar% {pintersect(findOverlapPairs(x, ALL_REGIONS))} # Get GRanges information for those regions
    myGRList <- foreach(x=myGRList, .final=function(x) setNames(x, names(myGRList))) %dopar% {names(x) <- 1:length(x); return(x)} # Reset names as some regions could have been discarded
  }
  stopifnot(all(foreach(x=2:length(myGRList), .combine=c) %dopar% {identical(granges(myGRList[[1]]), granges(myGRList[[x]]))})) # Make sure regions are strictly identical
  return(myGRList)
}
