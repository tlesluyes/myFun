#' @title computeLogR
#' @description Compute the theoretical logR value for a given segment
#' @details This function computes the theoretical logR value for a given segment (from nMajor, nMinor, purity and ploidy values). Since logR isn't allele-specific, ntot can be used instead of nMajor (and nMinor should set to 0).
#' @param nMajor the number of copies of the major allele
#' @param nMinor the number of copies of the minor allele
#' @param purity the purity estimate of the tumour
#' @param ploidy the ploidy estimate of the tumour
#' @param digits a numeric, the number of digits to round to (default: 4)
#' @return A number representing the logR value
#' @seealso https://doi.org/10.1038/s41592-020-01013-2
#' @examples
#' # A 2+1 state in a diploid tumour with 90% purity
#' computeLogR(2, 1, 0.9, 2)
#' # A loss of 1 copy (2+1) in a pseudo-tetraploid tumour with 60% purity
#' computeLogR(2, 1, 0.6, 3.5)
#' @author tlesluyes
#' @export
computeLogR <- function(nMajor, nMinor, purity, ploidy, digits=4) {
  stopifnot(all(is.numeric(c(nMajor, nMinor, purity, ploidy, digits))))
  round(log2((purity*(nMajor+nMinor)+2*(1-purity))/(purity*ploidy+2*(1-purity))), digits)
}

#' @title computeBAF
#' @description Compute the theoretical BAF values for a given segment
#' @details This function computes the theoretical BAF values for a given segment (from nMajor, nMinor and purity values).
#' @param nMajor the number of copies of the major allele
#' @param nMinor the number of copies of the minor allele
#' @param purity the purity estimate of the tumour
#' @param digits a numeric, the number of digits to round to (default: 4)
#' @return A vector of two numbers representing the BAF values
#' @seealso https://doi.org/10.1038/s41592-020-01013-2
#' @examples
#' # A 2+1 state in a tumour with 90% purity
#' computeBAF(2, 1, 0.9)
#' # A 1+0 state in a tumour with 60% purity
#' computeBAF(1, 0, 0.6)
#' @author tlesluyes
#' @export
computeBAF <- function(nMajor, nMinor, purity, digits=4) {
  stopifnot(all(is.numeric(c(nMajor, nMinor, purity, digits))))
  BAF1 <- (purity*nMinor+(1-purity))/(purity*(nMajor+nMinor)+2*(1-purity))
  BAF2 <- (purity*nMajor+(1-purity))/(purity*(nMajor+nMinor)+2*(1-purity))
  return(round(c(BAF1, BAF2), digits))
}

#' @title computeFit
#' @description Compute the purity/ploidy fit for a given segment
#' @details This function computes the purity/ploidy fit (rho, psi and psit) for a given segment (from logR, BAF, proposed nMajor and proposed nMinor).
#' @param logR the logR value of the segment
#' @param BAF the BAF value of the segment (upper band only so the value should be in the 0.5-1 space)
#' @param nMajor the number of copies of the major allele
#' @param nMinor the number of copies of the minor allele
#' @param gamma the gamma parameter is platform-dependent and represents the expected logR decrease in a diploid sample where one copy is lost (should be 1 for HTS data and 0.55 for SNP arrays)
#' @param digits a numeric, the number of digits to round to (default: 4)
#' @return A list with the rho (=purity), psi (=total ploidy) and psit (=tumour ploidy) values
#' @examples
#' # A segment has logR=0.5361 and BAF=0.3448/0.6552
#' # What is the purity/ploidy fit if I believe that the segment is 2+1?
#' computeFit(0.5361, 0.6552, 2, 1, 1) # purity=90%; ploidy=2
#' @author tlesluyes
#' @export
#' @seealso https://doi.org/10.1038/s41592-020-01013-2
computeFit <- function(logR, BAF, nMajor, nMinor, gamma, digits=4) {
  stopifnot(all(is.numeric(c(logR, BAF, nMajor, nMinor, gamma, digits))))
  rho <- (2*BAF-1)/(2*BAF-BAF*(nMajor+nMinor)-1+nMajor)
  psi <- (rho*(nMajor+nMinor)+2-2*rho)/(2^(logR/gamma))
  psit <- (psi-2*(1-rho))/rho
  return(list(rho=round(rho, digits), psi=round(psi, digits), psit=round(psit, digits)))
}

#' @title reestimate_ploidy
#' @description Compute the re-estimated ploidy for a given sample
#' @details This function computes the re-estimated ploidy for a given sample (from its old purity/ploidy fit and the re-estimated purity).
#' @param rho.old old purity estimate
#' @param psit.old old ploidy estimate
#' @param rho.new new purity estimate
#' @param WGD number of WGD events (0 if there is no WGD)
#' @param digits a numeric, the number of digits to round to (default: 4)
#' @return A number representing the re-estimated ploidy
#' @examples
#' # A pseudo-diploid sample has purity=74% and ploidy=2.4
#' # What is the re-estimated ploidy if I believe that the sample has purity=61%?
#' reestimate_ploidy(0.74, 2.4, 0.61, 0)
#' @author tlesluyes
#' @export
reestimate_ploidy <- function(rho.old, psit.old, rho.new, WGD, digits=4) {
  stopifnot(all(is.numeric(c(rho.old, psit.old, rho.new, digits))))
  COEF <- 2*(WGD+1)
  return(round(((rho.old*psit.old)+COEF*(rho.new-rho.old))/rho.new, digits))
}

#' @title reestimate_purity
#' @description Compute the re-estimated purity for a given sample
#' @details This function computes the re-estimated purity for a given sample in the context of a jump in ploidy (so the matched ploidy needs to be doubled or halved).
#' @param rho.old old purity estimate
#' @param psit.old old ploidy estimate
#' @param switch a character ("double" or "halve") indicating whether the ploidy should be doubled or halved
#' @param digits a numeric, the number of digits to round to (default: 4)
#' @return A number representing the re-estimated purity
#' @examples
#' # A sample has purity=74% and ploidy=2.4 but the CNA profile needs to be doubled
#' # What is the re-estimated purity?
#' reestimate_purity(0.74, 2.4, "double")
#' @author tlesluyes
#' @export
reestimate_purity <- function(rho.old, psit.old, switch, digits=4) {
  stopifnot(all(is.numeric(c(rho.old, psit.old, digits))))
  stopifnot(switch %in% c("double", "halve"))
  if (switch=="double") {
    return(round(rho.old/(2-rho.old), digits))
  } else {
    return(round(2*rho.old/(rho.old+1), digits))
  }
}
