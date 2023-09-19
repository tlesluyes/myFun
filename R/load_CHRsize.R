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
    message("Loading hg19 data")
    #utils::data("CHRsize_hg19", package="myFun")
    load(system.file("extdata", "CHRsize_hg19.rda", package="myFun"), envir=.GlobalEnv)
  } else if (assembly=="hg38") {
    message("Loading hg38 data")
    #utils::data("CHRsize_hg38", package="myFun")
    load(system.file("extdata", "CHRsize_hg38.rda", package="myFun"), envir=.GlobalEnv)
  } else {
    stop("Unsupported assembly")
  }
}