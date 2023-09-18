#' @title computeMD
#' @description Compute the Manhattan distance (MD)
#' @details This function computes the Manhattan distance (MD) between two profiles (as GRanges objects).
#' @param GR1 a GRanges object corresponding to a single CNA profile
#' @param GR2 a GRanges object corresponding to a single CNA profile
#' @param nMajor a metadata column name for the major allele (default: "nMajor")
#' @param nMinor a metadata column name for the minor allele (default: "nMinor")
#' @param convertMb a boolean, the MD will be converted to megabases if set to TRUE (default: FALSE)
#' @return A numeric value representing the MD
#' @examples
#' require(GenomicRanges)
#' GR1=GRanges(seqnames=rep("1", 3),
#'             ranges=IRanges(start=c(1, 1001, 10001),end=c(1000, 10000, 20000)),
#'             nMajor=c(1, 2, 1),
#'             nMinor=c(1, 1, 1))
#' GR2=GRanges(seqnames=rep("1", 2),
#'             ranges=IRanges(start=c(500, 10001),end=c(10000, 25000)),
#'             nMajor=c(2, 1),
#'             nMinor=c(1, 1))
#' # in this example:
#' #    Region 500-1000 (size=501) is 1+1 for GR1 and 2+1 for GR2
#' #    Region 1001-20000 (size=19000) is identical between GR1 and GR2 (both 2+1 and 1+1)
#' #    MD is: (abs(2-1)+abs(1-1))*501 = 501
#' computeMD(GR1, GR2)
#' @author tlesluyes
#' @export
computeMD=function(GR1, GR2, nMajor="nMajor", nMinor="nMinor", convertMb=FALSE) {
  require(GenomicRanges)
  stopifnot(length(convertMb)==1 && is.logical(convertMb))
  profiles=harmonizeGRanges(list(GR1, GR2))
  stopifnot(all(sapply(profiles, function(x) all(c(nMajor, nMinor) %in% names(mcols(x))))))
  stopifnot(all(sapply(profiles, function(x) is.numeric(mcols(x)[, nMajor]))))
  stopifnot(all(sapply(profiles, function(x) is.numeric(mcols(x)[, nMinor]))))
  MD=sum((abs(mcols(profiles[[1]])[, nMajor]-mcols(profiles[[2]])[, nMajor])+abs(mcols(profiles[[1]])[, nMinor]-mcols(profiles[[2]])[, nMinor]))*width(profiles[[1]]))
  if (convertMb) MD=MD/1e6
  return(MD)
}