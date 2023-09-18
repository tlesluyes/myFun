#' @title load_CHRsize
#' @description Load `CHRsize` information
#' @details This function loads `CHRsize` information for a given assembly. It is then available as a data.frame called `CHRsize` in the environment.
#' @param assembly an assembly (hg19 or hg38)
#' @return A data.frame with the `CHRsize` information
#' @examples load_CHRsize("hg38"); head(CHRsize)
#' @author tlesluyes
#' @export
load_CHRsize=function(assembly) {
  if (assembly=="hg19") {
    data("CHRsize_hg19", package="myFun")
  } else if (assembly=="hg38") {
    data("CHRsize_hg38", package="myFun")
  } else {
    stop("Unsupported assembly")
  }
}