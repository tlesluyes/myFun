#' @title computeISA
#' @description Compute the inter-sample agreement (ISA)
#' @details This function computes the inter-sample agreement (ISA) between two profiles (as GRanges objects). This corresponds to the fraction of the genome (%) with the same CN status.
#' @param GR1 a GRanges object corresponding to a single CNA profile
#' @param GR2 a GRanges object corresponding to a single CNA profile
#' @param CNstatus a metadata column name for the copy-number status (default: "CNstatus"). Can be total (e.g. "3") or allele-specific (e.g. "2+1")
#' @return A percentage representing the ISA
#' @examples
#' GR1 <- GenomicRanges::GRanges(seqnames=rep("1", 3),
#'                               ranges=IRanges::IRanges(start=c(1, 1001, 10001),
#'                                                       end=c(1000, 10000, 20000)),
#'                               CNstatus=c("1+1", "2+1", "1+1"))
#' GR2 <- GenomicRanges::GRanges(seqnames=rep("1", 2),
#'                               ranges=IRanges::IRanges(start=c(500, 10001),
#'                                                       end=c(10000, 25000)),
#'                               CNstatus=c("2+1", "1+1"))
#' # in this example:
#' #    Region 500-1000 (size=501) is 1+1 for GR1 and 2+1 for GR2
#' #    Region 1001-20000 (size=19000) is identical between GR1 and GR2 (both 2+1 and 1+1)
#' #    ISA is: 19000/19501 = 97.43%
#' computeISA(GR1, GR2)
#' @author tlesluyes
#' @export
computeISA <- function(GR1, GR2, CNstatus="CNstatus") {
  checkGRlist(list(GR1, GR2))
  stopifnot(all(sapply(list(GR1, GR2), function(x) CNstatus %in% names(mcols(x)))))
  profiles <- harmonizeGRanges(list(GR1, GR2))
  profiles <- cleanGRlistMetadata(profiles, metadataCols=CNstatus)
  sameCN <- which(mcols(profiles[[1]])[, CNstatus]==mcols(profiles[[2]])[, CNstatus])
  return(sum(width(profiles[[1]][sameCN, ]))/sum(width(profiles[[1]]))*100)
}

#' @title computeISA_batch
#' @description Compute the inter-sample agreement (ISA) for a batch of samples
#' @details This function computes the inter-sample agreement (ISA) between multiple profiles (as a list of GRanges objects).
#' @param myGRList a list of GRanges objects, each object should correspond to one CNA profile
#' @param cores a numeric, the number of cores to use (default: 1)
#' @param min_seg_size a numeric, the minimum segment size (in bp) to consider (default: 0)
#' @param CNstatus a metadata column name for the copy-number status (default: "CNstatus"). Can be total (e.g. "3") or allele-specific (e.g. "2+1")
#' @return A matrix of ISA values
#' @examples
#' GR1 <- GenomicRanges::GRanges(seqnames=rep("1", 3),
#'                               ranges=IRanges::IRanges(start=c(1, 1001, 10001),
#'                                                       end=c(1000, 10000, 20000)),
#'                               CNstatus=c("1+1", "2+1", "1+1"))
#' GR2 <- GenomicRanges::GRanges(seqnames=rep("1", 2),
#'                               ranges=IRanges::IRanges(start=c(500, 10001),
#'                                                       end=c(10000, 25000)),
#'                               CNstatus=c("2+1", "1+1"))
#' GR3 <- GenomicRanges::GRanges(seqnames="1",
#'                               ranges=IRanges::IRanges(start=500,
#'                                                       end=25000),
#'                               CNstatus="1+1")
#' myGRList <- list(GR1, GR2, GR3)
#' names(myGRList) <- c("GR1", "GR2", "GR3")
#' computeISA_batch(myGRList)
#' @author tlesluyes
#' @export
computeISA_batch <- function(myGRList, cores=1, min_seg_size=0, CNstatus="CNstatus") {
  if (cores>1) registerDoParallel(cores=cores)
  checkGRlist(myGRList)
  stopifnot(all(sapply(myGRList, function(x) all("CNstatus" %in% names(mcols(x))))))
  if (is.null(names(myGRList))) {
    NAMES <- as.character(1:length(myGRList))
  } else {
    NAMES <- names(myGRList)
  }
  myGRList <- harmonizeGRanges(myGRList, cores=cores)
  myGRList <- cleanGRlistMetadata(myGRList, metadataCols=CNstatus)
  if (min_seg_size>0) {
    TO_KEEP <- width(myGRList[[1]])>min_seg_size
    myGRList <- lapply(myGRList, function(x) x[TO_KEEP, ])
    rm(TO_KEEP)
  }
  CNstatus_matrix <- do.call(cbind, lapply(myGRList, function(x) { return(mcols(x)[, CNstatus]) }))
  WIDTHS <- width(myGRList[[1]])
  SUM_WIDTH <- sum(WIDTHS)
  ISA <- foreach(i=1:length(myGRList), .combine=cbind) %dopar% {
    i <- get("i")
    if (i %% 50==0) print(paste0(i, "/", length(myGRList)))
    SAME <- CNstatus_matrix==CNstatus_matrix[, i]
    return(sapply(1:length(myGRList), function(x) sum(WIDTHS[SAME[, x]])/SUM_WIDTH))
  }
  ISA <- ISA*100
  rownames(ISA) <- NAMES
  colnames(ISA) <- NAMES
  return(ISA)
}
