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
#' GR1 <- GenomicRanges::GRanges(seqnames=rep("1", 3),
#'                               ranges=IRanges::IRanges(start=c(1, 1001, 10001),
#'                                                       end=c(1000, 10000, 20000)),
#'                               nMajor=c(1, 2, 1),
#'                               nMinor=c(1, 1, 1))
#' GR2 <- GenomicRanges::GRanges(seqnames=rep("1", 2),
#'                               ranges=IRanges::IRanges(start=c(500, 10001),
#'                                                       end=c(10000, 25000)),
#'                               nMajor=c(2, 1),
#'                               nMinor=c(1, 1))
#' # in this example:
#' #    Region 500-1000 (size=501) is 1+1 for GR1 and 2+1 for GR2
#' #    Region 1001-20000 (size=19000) is identical between GR1 and GR2 (both 2+1 and 1+1)
#' #    MD is: (abs(2-1)+abs(1-1))*501 = 501
#' computeMD(GR1, GR2)
#' @author tlesluyes
#' @export
computeMD <- function(GR1, GR2, nMajor="nMajor", nMinor="nMinor", convertMb=FALSE) {
  checkGRlist(list(GR1, GR2))
  stopifnot(length(convertMb)==1 && is.logical(convertMb))
  stopifnot(all(sapply(list(GR1, GR2), function(x) all(c(nMajor, nMinor) %in% names(mcols(x))))))
  stopifnot(all(sapply(list(GR1, GR2), function(x) is.numeric(mcols(x)[, nMajor]))))
  stopifnot(all(sapply(list(GR1, GR2), function(x) is.numeric(mcols(x)[, nMinor]))))
  profiles <- harmonizeGRanges(list(GR1, GR2))
  MD <- sum((abs(mcols(profiles[[1]])[, nMajor]-mcols(profiles[[2]])[, nMajor])+abs(mcols(profiles[[1]])[, nMinor]-mcols(profiles[[2]])[, nMinor]))*width(profiles[[1]]))
  if (convertMb) MD <- MD/1e6
  return(MD)
}

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
#' GR1 <- GenomicRanges::GRanges(seqnames=rep("1", 3),
#'                               ranges=IRanges::IRanges(start=c(1, 1001, 10001),
#'                                                       end=c(1000, 10000, 20000)),
#'                               nMajor=c(1, 2, 1),
#'                               nMinor=c(1, 1, 1))
#' GR2 <- GenomicRanges::GRanges(seqnames=rep("1", 2),
#'                               ranges=IRanges::IRanges(start=c(500, 10001),
#'                                                       end=c(10000, 25000)),
#'                               nMajor=c(2, 1),
#'                               nMinor=c(1, 1))
#' GR3 <- GenomicRanges::GRanges(seqnames="1",
#'                               ranges=IRanges::IRanges(start=500,
#'                                                       end=25000),
#'                               nMajor=1,
#'                               nMinor=1)
#' myGRList <- list(GR1, GR2, GR3)
#' names(myGRList) <- c("GR1", "GR2", "GR3")
#' computeMD_batch(myGRList)
#' @author tlesluyes
#' @export
computeMD_batch <- function(myGRList, cores=1, min_seg_size=0, nMajor="nMajor", nMinor="nMinor", convertMb=FALSE) {
  if (cores>1) registerDoParallel(cores=cores)
  checkGRlist(myGRList)
  stopifnot(length(convertMb)==1 && is.logical(convertMb))
  stopifnot(all(sapply(myGRList, function(x) all(c(nMajor, nMinor) %in% names(mcols(x))))))
  stopifnot(all(sapply(myGRList, function(x) is.numeric(mcols(x)[, nMajor]))))
  stopifnot(all(sapply(myGRList, function(x) is.numeric(mcols(x)[, nMinor]))))
  if (is.null(names(myGRList))) {
    NAMES <- as.character(1:length(myGRList))
  } else {
    NAMES <- names(myGRList)
  }
  myGRList <- harmonizeGRanges(myGRList, cores=cores)
  if (min_seg_size>0) {
    TO_KEEP <- width(myGRList[[1]])>min_seg_size
    myGRList <- lapply(myGRList, function(x) x[TO_KEEP, ])
    rm(TO_KEEP)
  }
  nMajor_matrix <- do.call(cbind, lapply(myGRList, function(x) { return(mcols(x)[, nMajor]) }))
  nMinor_matrix <- do.call(cbind, lapply(myGRList, function(x) { return(mcols(x)[, nMinor]) }))
  WIDTHS <- width(myGRList[[1]])
  MD <- foreach(i=1:length(myGRList), .combine=cbind) %dopar% {
    i <- get("i")
    if (i %% 50==0) print(paste0(i, "/", length(myGRList)))
    return(apply((abs(nMajor_matrix-nMajor_matrix[, i])+abs(nMinor_matrix-nMinor_matrix[, i]))*WIDTHS, 2, sum))
  }
  rownames(MD) <- NAMES
  colnames(MD) <- NAMES
  if (convertMb) MD <- MD/1e6
  return(MD)
}
