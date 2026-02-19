#' @title example_CNAs
#' @description Generate random CNA profiles
#' @details This function generates random CNA profiles for a given assembly.
#' @param n an integer, the number of CNA profiles to generate (default: 1)
#' @param assembly a character, the assembly to use (default: "hg38")
#' @param seed an integer, the seed for reproducibility (default: 1234)
#' @return A list of GRanges objects, each corresponding to one CNA profile
#' @examples
#' example_CNAs(n=3, assembly="hg38")
#' @author tlesluyes
#' @export
example_CNAs <- function(n=1, assembly="hg38", seed=1234) {
  load_CHRsize(assembly)
  CHRsize <- CHRsize[1:22, ]
  mydefaultGR <- GRanges(seqnames=CHRsize$chr,
                         ranges=IRanges(start=rep(1, nrow(CHRsize)), end=CHRsize$size),
                         nMajor=1,
                         nMinor=1,
                         gain=FALSE,
                         loss=FALSE,
                         seqinfo = get_Seqinfo(assembly))
  chr_arms <- get_arms(assembly)
  set.seed(seed)
  myGRlist <- lapply(1:n, function(i) {
    myGR <- mydefaultGR
    nEvents <- sample(1:10, 1)
    parameters <- data.frame(chr=sample(as.character(seqnames(mydefaultGR)), nEvents, replace=FALSE),
                             category=sample(c("whole", "p", "q"), prob = c(0.5, 0.25, 0.25), nEvents, replace = TRUE),
                             gain=sample(c(TRUE, FALSE), nEvents, replace = TRUE))
    parameters$loss <- ! parameters$gain
    CNAs <- Reduce(c, lapply(1:nrow(parameters), function(i) {
      GRanges(seqnames = parameters$chr[i],
              ranges = IRanges(start=ifelse(parameters$category[i] %in% c("whole", "p"), 1, chr_arms$start[chr_arms$arm==paste0(parameters$chr[i], parameters$category[i])]),
                               end=ifelse(parameters$category[i] %in% c("whole", "q"), CHRsize$size[which(CHRsize$chr==parameters$chr[i])], chr_arms$end[chr_arms$arm==paste0(parameters$chr[i], parameters$category[i])])),
              nMajor=1+as.numeric(parameters$gain[i]),
              nMinor=1-as.numeric(parameters$loss[i]),
              gain=parameters$gain[i],
              loss=parameters$loss[i],
              seqinfo = seqinfo(mydefaultGR))
    }))
    CNAs <- sort(CNAs)
    myGR <- sort(c(excludeGRanges(myGR, CNAs), CNAs))
    myGR$CNstatus <- paste0(myGR$nMajor, "+", myGR$nMinor)
    stopifnot(checkGR(myGR))
    return(myGR)
  })
  names(myGRlist) <- paste0("G", 1:n)
  return(myGRlist)
}
