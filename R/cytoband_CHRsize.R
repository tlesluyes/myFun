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
#' @param assembly an assembly (hg19, hg38 or CHM13)
#' @return A data.frame with the `CHRsize` information
#' @examples load_CHRsize("hg38"); head(CHRsize)
#' @author tlesluyes
#' @export
load_CHRsize <- function(assembly) {
  stopifnot(length(assembly)==1 && is.character(assembly))
  if (assembly=="hg19") {
    if (environmentName(parent.frame())=="R_GlobalEnv") message("Loading hg19 CHRsize data")
    load(system.file("extdata", "CHRsize_hg19.rda", package="myFun"), envir=parent.frame())
  } else if (assembly=="hg38") {
    if (environmentName(parent.frame())=="R_GlobalEnv") message("Loading hg38 CHRsize data")
    load(system.file("extdata", "CHRsize_hg38.rda", package="myFun"), envir=parent.frame())
  } else if (assembly=="CHM13") {
    if (environmentName(parent.frame())=="R_GlobalEnv") message("Loading CHM13 (v2.0) CHRsize data")
    load(system.file("extdata", "CHRsize_CHM13.rda", package="myFun"), envir=parent.frame())
  } else {
    stop("Unsupported assembly")
  }
}

#' @title load_cytoband
#' @description Load `cytoband` information
#' @details This function loads `cytoband` information for a given assembly. It is then available as a data.frame called `cytoband` in the environment.
#' @param assembly an assembly (hg19, hg38 or CHM13)
#' @return A data.frame with the `cytoband` information
#' @examples load_cytoband("hg38"); head(cytoband)
#' @author tlesluyes
#' @export
load_cytoband <- function(assembly) {
  stopifnot(length(assembly)==1 && is.character(assembly))
  if (assembly=="hg19") {
    if (environmentName(parent.frame())=="R_GlobalEnv") message("Loading hg19 cytoband data")
    load(system.file("extdata", "cytoband_hg19.rda", package="myFun"), envir=parent.frame())
  } else if (assembly=="hg38") {
    if (environmentName(parent.frame())=="R_GlobalEnv") message("Loading hg38 cytoband data")
    load(system.file("extdata", "cytoband_hg38.rda", package="myFun"), envir=parent.frame())
  } else if (assembly=="CHM13") {
    if (environmentName(parent.frame())=="R_GlobalEnv") message("Loading CHM13 (v2.0) cytoband data")
    load(system.file("extdata", "cytoband_CHM13.rda", package="myFun"), envir=parent.frame())
  } else {
    stop("Unsupported assembly")
  }
}

#' @title get_arms
#' @description Get chromosome arms information
#' @details This function gets chromosome arms information for a given assembly.
#' @param assembly an assembly (e.g. hg38) or a data.frame with expected cytoband information (chr, start, end, name, gieStain)
#' @param withCentromeres whether to include centromeric regions (default: TRUE)
#' @return A data.frame with genomic regions, can be converted to GRanges on the fly
#' @examples get_arms("hg38")
#' @author tlesluyes
#' @export
get_arms <- function(assembly, withCentromeres=TRUE) {
  if (length(assembly)==1 && is.character(assembly)) {
    load_cytoband(assembly)
  } else if (is.data.frame(assembly)) {
    cytoband <- assembly
  } else {
    stop("Unsupported input")
  }
  stopifnot(length(withCentromeres)==1 && is.logical(withCentromeres))
  stopifnot(all(c("chr", "start", "end", "name") %in% colnames(cytoband)))
  stopifnot(all(grepl("p|q", cytoband$name)))
  if (isFALSE(withCentromeres)) {
    stopifnot("gieStain" %in% colnames(cytoband))
    cytoband <- cytoband[which(cytoband$gieStain!="acen"), ]
  }
  cytoband$arm <- paste0(cytoband$chr, gsub("^(p|q).*$", "\\1", cytoband$name))
  cytoband <- split(cytoband, factor(cytoband$arm, levels=unique(cytoband$arm)))
  cytoband <- lapply(cytoband, function(x) {
    x <- x[order(x$start, x$end), ]
    return(data.frame(chr=unique(x$chr), start=min(x$start), end=max(x$end), arm=unique(x$arm)))
  })
  cytoband <- do.call(rbind, cytoband)
  return(cytoband)
}

#' @title get_centromeres
#' @description Get centromeres information
#' @details This function gets centromeres information for a given assembly.
#' @param assembly an assembly (e.g. hg38) or a data.frame with expected cytoband information (chr, start, end, gieStain)
#' @return A data.frame with genomic regions, can be converted to GRanges on the fly
#' @examples get_centromeres("hg38")
#' @author tlesluyes
#' @export
get_centromeres <- function(assembly) {
  if (length(assembly)==1 && is.character(assembly)) {
    load_cytoband(assembly)
  } else if (is.data.frame(assembly)) {
    cytoband <- assembly
  } else {
    stop("Unsupported input")
  }
  stopifnot(all(c("chr", "start", "end", "gieStain") %in% colnames(cytoband)))
  cytoband <- cytoband[which(cytoband$gieStain=="acen"), ]
  cytoband <- split(cytoband, factor(cytoband$chr, levels=unique(cytoband$chr)))
  cytoband <- lapply(cytoband, function(x) {
    return(data.frame(chr=unique(x$chr), start=min(x$start), end=max(x$end)))
  })
  cytoband <- do.call(rbind, cytoband)
  return(cytoband)
}

#' @title get_Seqinfo
#' @description Get Seqinfo object
#' @details This function gets Seqinfo object for a given assembly.
#' @param assembly an assembly (e.g. hg38) or a data.frame with expected cytoband information (chr, size)
#' @param genome genome name (e.g. hg38), required if `assembly` is a data.frame (default: NA)
#' @return A Seqinfo object
#' @examples get_Seqinfo("hg38")
#' @author tlesluyes
#' @export
get_Seqinfo <- function(assembly, genome=NA) {
  if (length(assembly)==1 && is.character(assembly)) {
    load_CHRsize(assembly)
    genome <- assembly
  } else if (is.data.frame(assembly)) {
    CHRsize <- assembly
    stopifnot(length(genome)==1 && is.character(genome))
  } else {
    stop("Unsupported input")
  }
  stopifnot(all(c("chr", "size") %in% colnames(CHRsize)))
  myseqinfo <- Seqinfo(seqnames=CHRsize$chr, seqlengths=CHRsize$size, isCircular=rep(FALSE, nrow(CHRsize)), genome=genome)
  return(myseqinfo)
}

#' @title seqinfo2GR
#' @description Convert Seqinfo to GRanges
#' @details This function generates a GRanges object from a Seqinfo object.
#' @param myseqinfo a Seqinfo object
#' @return A GRanges object
#' @examples seqinfo2GR(get_Seqinfo("hg38"))
#' @author tlesluyes
#' @export
seqinfo2GR <- function(myseqinfo) {
  stopifnot(is(myseqinfo, "Seqinfo"))
  return(GRanges(seqnames=seqnames(myseqinfo), ranges=IRanges(start=1, end=seqlengths(myseqinfo)), seqinfo=myseqinfo))
}

#' @title add_arm_info
#' @description Add arm information to a GRanges object
#' @details This function adds arm information (e.g. 1q, 17p) to a GRanges object.
#' @param myGR a GRanges object
#' @param assembly an assembly (e.g. hg38)
#' @return A GRanges object with arm information added
#' @examples add_arm_info(GenomicRanges::GRanges(seqnames=c("1", "17"),
#'                                               ranges=IRanges::IRanges(start=c(150e6, 1),
#'                                                                       end=c(200e6, 20e6))),
#'                        "hg38")
#' @author tlesluyes
#' @export
add_arm_info <- function(myGR, assembly) {
  stopifnot(checkGR(myGR))
  myseqinfo <- get_Seqinfo(assembly)
  if (any(is.na(genome(myGR)))) {
    seqlevels(myGR) <- seqlevels(myseqinfo)
    seqinfo(myGR) <- myseqinfo
  }
  chr_arms <- makeGRangesFromDataFrame(get_arms(assembly), seqinfo = myseqinfo, keep.extra.columns = TRUE)
  myGR <- harmonizeGRanges(list(myGR=myGR, chr_arms=chr_arms))
  myGR[["myGR"]]$arm <- myGR[["chr_arms"]]$arm
  myGR <- myGR[["myGR"]]
  return(myGR)
}
