#' @title adjustPositions
#' @description Adjust genomic positions
#' @details This function adjusts genomic positions according to the chromosome sizes. The first nucleotide of chromosome 2 corresponds to the size of the chromosome 1 + 1bp and so on.
#' @param DF a data.frame
#' @param CHRsize a data.frame from the `load_CHRsize` function
#' @param chr_column a column name with chromosome information (default: "chr")
#' @param start_column a column name with start position (default: "start")
#' @param end_column a column name with end position (default: "end"). If end does not exist (e.g. SNP positions), it should be set to the same column as start so it will be ignored
#' @param suffix a suffix for the adjusted positions (default: "_adj")
#' @return A data.frame with adjusted genomic positions
#' @examples
#' DF <- data.frame(chr=c(1:3), start=rep(1e6, 3), end=rep(125e6, 3))
#' load_CHRsize("hg19")
#' adjustPositions(DF, CHRsize)
#' @author tlesluyes
#' @export
adjustPositions <- function(DF, CHRsize, chr_column="chr", start_column="start", end_column="end", suffix="_adj") {
  stopifnot(is.data.frame(DF))
  stopifnot(is.data.frame(CHRsize))
  stopifnot(all(c(chr_column, start_column, end_column) %in% colnames(DF)))
  stopifnot(all(c("chr", "add") %in% colnames(CHRsize)))
  stopifnot(all(gsub("^chr", "", DF[, chr_column]) %in% CHRsize$chr))
  DF[, paste0(start_column, suffix)] <- DF[, start_column]
  if (end_column != start_column) DF[, paste0(end_column, suffix)] <- DF[, end_column]
  matches <- match(gsub("^chr", "", DF[, chr_column]), CHRsize$chr)
  stopifnot(! any(is.na(matches)))
  DF[, paste0(start_column, suffix)] <- DF[, paste0(start_column, suffix)] + CHRsize$add[matches]
  if (end_column != start_column) DF[, paste0(end_column, suffix)] <- DF[, paste0(end_column, suffix)] + CHRsize$add[matches]
  return(DF)
}
