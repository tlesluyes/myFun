#' @title load_cytoband
#' @description Load `cytoband` information
#' @details This function loads `cytoband` information for a given assembly. It is then available as a data.frame called `cytoband` in the environment.
#' @param assembly an assembly (hg19 or hg38)
#' @return A data.frame with the `cytoband` information
#' @examples load_cytoband("hg38"); head(cytoband)
#' @author tlesluyes
#' @export
load_cytoband=function(assembly) {
  if (assembly=="hg19") {
    data("cytoband_hg19", package="myFun")
  } else if (assembly=="hg38") {
    data("cytoband_hg38", package="myFun")
  } else {
    stop("Unsupported assembly")
  }
}