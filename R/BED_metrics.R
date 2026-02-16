#' @title BED_metrics
#' @description Give BED metrics
#' @details This function provides several metrics of interest from a BED file/object.
#' @param BED a data.frame, a GRanges object or a path to a BED file
#' @param verbose a boolean, whether to print metrics (default: TRUE)
#' @return A named list of metrics (number of chromosomes, number of regions and total size of the regions) before and after removing overlaps (`GenomicRanges::reduce()`). Strand information is not considered
#' @examples
#' # BED format is 0-based for starts!
#' BED <- data.frame(chr=c(c(1,1:3)), start=c(0, 1e3, 0, 1e3), end=c(2e3, 5e3, 2e3, 4e3))
#' BED_metrics(BED)
#' @author tlesluyes
#' @export
BED_metrics <- function(BED, verbose=TRUE) {
  stopifnot(length(verbose)==1 && is.logical(verbose))
  stopifnot(is.data.frame(BED) || class(BED)=="GRanges" || file.exists(BED))
  if (is.character(BED)) {
    BED <- read.table(BED, sep="\t", header=FALSE, stringsAsFactors=FALSE)
  }
  if (is.data.frame(BED)) {
    stopifnot(nrow(BED)>=1)
    stopifnot(ncol(BED)>=3)
    BED <- BED[, 1:3]
    BED <- makeGRangesFromDataFrame(BED, seqnames.field=colnames(BED)[1], start.field=colnames(BED)[2], end.field=colnames(BED)[3], starts.in.df.are.0based = TRUE)
  }
  metrics <- list()
  metrics$chromosomes <- length(unique(seqnames(BED)))
  metrics$regions <- length(BED)
  metrics$size <- sum(width(BED))
  BED <- reduce(BED)
  metrics$chromosomes_reduced <- length(unique(seqnames(BED)))
  metrics$regions_reduced <- length(BED)
  metrics$size_reduced <- sum(width(BED))
  if (verbose) {
    print("BED metrics:")
    print(paste0("  Raw information: ", metrics$chromosomes, " chromosomes; ", metrics$regions, " regions; ", prettyNum(metrics$size, big.mark=",", scientific=FALSE), " bp"))
    print(paste0("       No overlap: ", metrics$chromosomes_reduced, " chromosomes; ", metrics$regions_reduced, " regions; ", prettyNum(metrics$size_reduced, big.mark=",", scientific=FALSE), " bp"))
  }
  invisible(metrics)
}
