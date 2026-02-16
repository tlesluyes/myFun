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
  cytoband$chr <- gsub("chr", "", cytoband$chr)
  cytoband$start_adj <- cytoband$start
  cytoband$end_adj <- cytoband$end
  for (i in c(2:22, "X", "Y")) {
    INDEX <- which(cytoband$chr==i)
    cytoband$start_adj[INDEX] <- cytoband$start_adj[INDEX]+CHRsize$add[which(CHRsize$chr==i)]
    cytoband$end_adj[INDEX] <- cytoband$end_adj[INDEX]+CHRsize$add[which(CHRsize$chr==i)]
    rm(INDEX)
  }; rm(i)
  cytoband <- cytoband[order(cytoband$start_adj), ]
  return(list(cytoband=cytoband, CHRsize=CHRsize))
}
