#' @title computeMD_batch
#' @description Compute the Manhattan distance (MD) for a batch of samples
#' @details This function computes the Manhattan distance (MD) between multiple profiles (as a list of GRanges objects).
#' @param myGRList a list of GRanges objects, each object should correspond to one CNA profile
#' @param cores a numeric, the number of cores to use (default: 1)
#' @param min_seg_size a numeric, the minimum segment size (in bp) to consider (default: 0)
#' @param nMajor a metadata column name for the major allele (default: "nMajor")
#' @param nMinor a metadata column name for the minor allele (default: "nMinor")
#' @param convertMb a boolean, the MD will be converted to megabases if set to TRUE (default: FALSE)
#' @return A matrix of MD values
#' @examples
#' require("GenomicRanges")
#' GR1=GRanges(seqnames=rep("1", 3),
#'             ranges=IRanges(start=c(1, 1001, 10001), end=c(1000, 10000, 20000)),
#'             nMajor=c(1, 2, 1),
#'             nMinor=c(1, 1, 1))
#' GR2=GRanges(seqnames=rep("1", 2),
#'             ranges=IRanges(start=c(500, 10001), end=c(10000, 25000)),
#'             nMajor=c(2, 1),
#'             nMinor=c(1, 1))
#' GR3=GRanges(seqnames="1",
#'             ranges=IRanges(start=500, end=25000),
#'             nMajor=1,
#'             nMinor=1)
#' myGRList=list(GR1, GR2, GR3)
#' names(myGRList)=c("GR1", "GR2", "GR3")
#' computeMD_batch(myGRList)
#' @author tlesluyes
#' @export
computeMD_batch=function(myGRList, cores=1, min_seg_size=0, nMajor="nMajor", nMinor="nMinor", convertMb=FALSE) {
  if (cores>1) doParallel::registerDoParallel(cores=cores)
  checkGRlist(myGRList)
  stopifnot(length(convertMb)==1 && is.logical(convertMb))
  stopifnot(all(sapply(myGRList, function(x) all(c(nMajor, nMinor) %in% names(GenomicRanges::mcols(x))))))
  stopifnot(all(sapply(myGRList, function(x) is.numeric(GenomicRanges::mcols(x)[, nMajor]))))
  stopifnot(all(sapply(myGRList, function(x) is.numeric(GenomicRanges::mcols(x)[, nMinor]))))
  if (is.null(names(myGRList))) {
    NAMES=as.character(1:length(myGRList))
  } else {
    NAMES=names(myGRList)
  }
  myGRList=harmonizeGRanges(myGRList, cores=cores)
  if (min_seg_size>0) {
    TO_KEEP=IRanges::width(myGRList[[1]])>min_seg_size
    myGRList=lapply(myGRList, function(x) x[TO_KEEP, ])
    rm(TO_KEEP)
  }
  nMajor_matrix=do.call(cbind, lapply(myGRList, function(x) { return(GenomicRanges::mcols(x)[, nMajor]) }))
  nMinor_matrix=do.call(cbind, lapply(myGRList, function(x) { return(GenomicRanges::mcols(x)[, nMinor]) }))
  WIDTHS=IRanges::width(myGRList[[1]])
  MD=foreach::foreach(i=1:length(myGRList), .combine=cbind) %dopar% {
    if (i %% 50==0) print(paste0(i, "/", length(myGRList)))
    return(apply((abs(nMajor_matrix-nMajor_matrix[, i])+abs(nMinor_matrix-nMinor_matrix[, i]))*WIDTHS, 2, sum))
  }
  rownames(MD)=NAMES
  colnames(MD)=NAMES
  if (convertMb) MD=MD/1e6
  return(MD)
}