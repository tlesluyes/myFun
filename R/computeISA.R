#' @title computeISA
#' @description Compute the inter-sample agreement (ISA)
#' @details This function computes the inter-sample agreement (ISA) between two profiles (as GRanges objects). This corresponds to the fraction of the genome (%) with the same CN status.
#' @param GR1 a GRanges object corresponding to a single CNA profile
#' @param GR2 a GRanges object corresponding to a single CNA profile
#' @param CNstatus a metadata column name for the copy-number status (default: "CNstatus"). Can be total (e.g. "3") or allele-specific (e.g. "2+1")
#' @return A percentage representing the ISA
#' @examples
#' require("GenomicRanges")
#' GR1=GRanges(seqnames=rep("1", 3),
#'             ranges=IRanges(start=c(1, 1001, 10001),end=c(1000, 10000, 20000)),
#'             CNstatus=c("1+1", "2+1", "1+1"))
#' GR2=GRanges(seqnames=rep("1", 2),
#'             ranges=IRanges(start=c(500, 10001),end=c(10000, 25000)),
#'             CNstatus=c("2+1", "1+1"))
#' # in this example:
#' #    Region 500-1000 (size=501) is 1+1 for GR1 and 2+1 for GR2
#' #    Region 1001-20000 (size=19000) is identical between GR1 and GR2 (both 2+1 and 1+1)
#' #    ISA is: 19000/19501 = 97.43%
#' computeISA(GR1, GR2)
#' @author tlesluyes
#' @export
computeISA=function(GR1, GR2, CNstatus="CNstatus") {
  profiles=harmonizeGRanges(list(GR1, GR2))
  stopifnot(all(sapply(profiles, function(x) CNstatus %in% names(GenomicRanges::mcols(x)))))
  sameCN=which(GenomicRanges::mcols(profiles[[1]])[, CNstatus]==GenomicRanges::mcols(profiles[[2]])[, CNstatus])
  return(sum(IRanges::width(profiles[[1]][sameCN, ]))/sum(IRanges::width(profiles[[1]]))*100)
}