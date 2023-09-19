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
    message("Loading hg19 data")
    #utils::data("cytoband_hg19", package="myFun")
    load(system.file("extdata", "cytoband_hg19.rda", package="myFun"), envir=.GlobalEnv)
  } else if (assembly=="hg38") {
    message("Loading hg38 data")
    #utils::data("cytoband_hg38", package="myFun")
    load(system.file("extdata", "cytoband_hg38.rda", package="myFun"), envir=.GlobalEnv)
  } else {
    stop("Unsupported assembly")
  }
}