#' @title adjustPositions
#' @description Adjust genomic positions
#' @details This function adjusts genomic positions according to the chromosome sizes. The first nucleotide of chromosome 2 corresponds to the size of the chromosome 1 + 1bp and so on.
#' @param input a GRanges object (its seqinfo must be set) or a data.frame
#' @param CHRsize a data.frame from the `load_CHRsize` function to be used if `input` is a data.frame (default: NULL)
#' @param chr_column a column name with chromosome information (default: "chr"; unused if `input` is a GRanges object)
#' @param start_column a column name with start position (default: "start"; unused if `input` is a GRanges object)
#' @param end_column a column name with end position (default: "end"; unused if `input` is a GRanges object). If end does not exist (e.g. SNP positions), it should be set to the same column as start so it will be ignored
#' @param suffix a suffix for the adjusted positions (default: "_adj")
#' @return An object (depending on the input type) with adjusted genomic positions
#' @examples
#' CNAprofile.GR <- example_CNAs(n=1, assembly="hg38")[[1]]
#' CNAprofile.DF <- data.frame(CNAprofile.GR)
#' load_CHRsize("hg38")
#' output.DF <- adjustPositions(CNAprofile.DF, CHRsize, chr_column="seqnames")
#' output.GR <- adjustPositions(CNAprofile.GR)
#' stopifnot(all(output.DF$start_adj == output.GR$start_adj))
#' stopifnot(all(output.DF$end_adj == output.GR$end_adj))
#' @author tlesluyes
#' @export
adjustPositions <- function(input, CHRsize=NULL, chr_column="chr", start_column="start", end_column="end", suffix="_adj") {
  if (is.data.frame(input)) {
    stopifnot(is.data.frame(CHRsize))
    stopifnot(all(c(chr_column, start_column, end_column) %in% colnames(input)))
    stopifnot(all(c("chr", "add") %in% colnames(CHRsize)))
    stopifnot(all(gsub("^chr", "", input[, chr_column]) %in% CHRsize$chr))
    input[, paste0(start_column, suffix)] <- input[, start_column]
    if (end_column != start_column) input[, paste0(end_column, suffix)] <- input[, end_column]
    matches <- match(gsub("^chr", "", input[, chr_column]), CHRsize$chr)
    stopifnot(! any(is.na(matches)))
    input[, paste0(start_column, suffix)] <- input[, paste0(start_column, suffix)] + CHRsize$add[matches]
    if (end_column != start_column) input[, paste0(end_column, suffix)] <- input[, paste0(end_column, suffix)] + CHRsize$add[matches]
  } else if (is(input, "GRanges")) {
    stopifnot(! any(is.na(seqlengths(input))))
    mcols(input)[, paste0("start", suffix)] <- start(input)
    if (any(start(input)!=end(input))) mcols(input)[, paste0("end", suffix)] <- end(input)
    CHRsize <- data.frame(chr=names(seqlengths(input)), size=seqlengths(input))
    CHRsize$add <- cumsum(as.numeric(c(0, CHRsize$size[1:(nrow(CHRsize)-1)])))
    matches <- match(as.character(seqnames(input)), CHRsize$chr)
    stopifnot(! any(is.na(matches)))
    mcols(input)[, paste0("start", suffix)] <- mcols(input)[, paste0("start", suffix)] + CHRsize$add[matches]
    if (any(start(input)!=end(input))) mcols(input)[, paste0("end", suffix)] <- mcols(input)[, paste0("end", suffix)] + CHRsize$add[matches]
  } else {
    stop("Unsupported input type")
  }
  return(input)
}
