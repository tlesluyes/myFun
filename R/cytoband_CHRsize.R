#' @title generate_cytoband_and_CHRsize
#' @description Generate `cytoband` and `CHRsize` information
#' @details This function generates `cytoband` and `CHRsize` information from a cytoband file. This can be obtained from the UCSC table browser -> select a genome/assembly -> "Mapping and Sequencing" -> "Chromosome Band" (not the ideogram version!) -> "get output" -> Remove the first "#" character (keep the header!).
#' @param cytoband_file a cytoband file
#' @return A list with both the `cytoband` and `CHRsize` information
#' @seealso `load_CHRsize("hg38"); load_cytoband("hg38")`
#' @author tlesluyes
#' @export
generate_cytoband_and_CHRsize <- function(cytoband_file) {
  cytoband <- read.table(paste0(cytoband_file), header=TRUE, stringsAsFactors=FALSE, sep="\t")
  colnames(cytoband)[1:3] <- c("chr", "start", "end")
  cytoband <- cytoband[which(cytoband$chr %in% paste0(rep(c("", "chr"), each=24), rep(c(1:22, "X", "Y"), 2))), ]
  cytoband$start <- cytoband$start+1
  COLORS <- data.frame(gieStain=c("acen", "gneg", paste0("gpos", seq(100, 25, -25)), "gvar", "stalk"),
                       colors=c("darkred", "grey90", "black", "grey40", "grey60", "grey80", "black", "grey60"),
                       row.names=1)
  stopifnot(all(cytoband$gieStain %in% rownames(COLORS)))
  cytoband$color <- COLORS[cytoband$gieStain, "colors"]
  rm(COLORS)
  CHRsize <- data.frame(number=1:24, chr=c(1:22, "X", "Y"), size=sapply(c(1:22, "X", "Y"), function(x) max(cytoband$end[which(cytoband$chr==paste0("chr", x))])), row.names=1)
  CHRsize$middle <- 0
  for (i in 1:24) {
    if (i==1) CHRsize$middle[i] <- CHRsize$size[i]/2
    else CHRsize$middle[i] <- sum(as.numeric(CHRsize$size[1:(i-1)]))+CHRsize$size[i]/2
  }; rm(i)
  CHRsize$sum <- cumsum(as.numeric(CHRsize$size))
  CHRsize$add <- cumsum(as.numeric(c(0, CHRsize$size[1:23])))
  cytoband$chr <- gsub("^chr", "", cytoband$chr)
  cytoband <- adjustPositions(cytoband, CHRsize)
  cytoband <- cytoband[order(cytoband$start_adj), ]
  return(list(cytoband=cytoband, CHRsize=CHRsize))
}

#' @title load_CHRsize
#' @description Load `CHRsize` information
#' @details This function loads `CHRsize` information for a given assembly. It is then available as a data.frame called `CHRsize` in the environment.
#' @param assembly an assembly (hg19 or hg38)
#' @return A data.frame with the `CHRsize` information
#' @examples load_CHRsize("hg38"); head(CHRsize)
#' @author tlesluyes
#' @export
load_CHRsize <- function(assembly) {
  if (assembly=="hg19") {
    message("Loading hg19 data")
    load(system.file("extdata", "CHRsize_hg19.rda", package="myFun"), envir=.GlobalEnv)
  } else if (assembly=="hg38") {
    message("Loading hg38 data")
    load(system.file("extdata", "CHRsize_hg38.rda", package="myFun"), envir=.GlobalEnv)
  } else {
    stop("Unsupported assembly")
  }
}

#' @title load_cytoband
#' @description Load `cytoband` information
#' @details This function loads `cytoband` information for a given assembly. It is then available as a data.frame called `cytoband` in the environment.
#' @param assembly an assembly (hg19 or hg38)
#' @return A data.frame with the `cytoband` information
#' @examples load_cytoband("hg38"); head(cytoband)
#' @author tlesluyes
#' @export
load_cytoband <- function(assembly) {
  if (assembly=="hg19") {
    message("Loading hg19 data")
    load(system.file("extdata", "cytoband_hg19.rda", package="myFun"), envir=.GlobalEnv)
  } else if (assembly=="hg38") {
    message("Loading hg38 data")
    load(system.file("extdata", "cytoband_hg38.rda", package="myFun"), envir=.GlobalEnv)
  } else {
    stop("Unsupported assembly")
  }
}
