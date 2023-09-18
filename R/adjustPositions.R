#' @title adjustPositions
#' @description Adjust genomic positions
#' @details This function adjusts genomic positions according to the chromosome sizes. The first nucleotide of chromosome 2 corresponds to the size of the chromosome 1 + 1bp and so on.
#' @param DF a data.frame
#' @param CHRsize a data.frame from the `load_CHRsize` function
#' @param chr_column a column name with chromosome information (default: "chr")
#' @param start_column a column name with start position (default: "start")
#' @param end_column a column name with end position (default: "end")
#' @param suffix a suffix for the adjusted positions (default: "_adj")
#' @return A data.frame with adjusted genomic positions
#' @examples
#' DF=data.frame(chr=c(1:3), start=rep(1e6, 3), end=rep(125e6, 3))
#' load_CHRsize("hg19")
#' adjustPositions(DF, CHRsize)
#' @author tlesluyes
#' @export
adjustPositions=function(DF, CHRsize, chr_column="chr", start_column="start", end_column="end", suffix="_adj") {
  stopifnot(all(c(chr_column, start_column, end_column) %in% colnames(DF)))
  DF[, chr_column]=gsub("^chr", "", DF[, chr_column])
  DF[, paste0(start_column, suffix)]=DF[, start_column]
  DF[, paste0(end_column, suffix)]=DF[, end_column]
  for (CHR in c(2:22, "X", "Y")) {
    INDEX=which(DF[, chr_column]==CHR)
    if (length(INDEX>0)) {
      DF[INDEX, paste0(start_column, suffix)]=DF[INDEX, paste0(start_column, suffix)]+CHRsize$add[which(CHRsize$chr==CHR)]
      DF[INDEX, paste0(end_column, suffix)]=DF[INDEX, paste0(end_column, suffix)]+CHRsize$add[which(CHRsize$chr==CHR)]
    }
    rm(INDEX)
  }; rm(CHR)
  return(DF)
}