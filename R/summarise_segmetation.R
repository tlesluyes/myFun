#' @title summarise_segmetation
#' @description Summarise segmentation data
#' @details This function summarises segmentation data, typically logR and/or BAF values for individual SNPs or loci.
#' @param DF a data.frame with segmentation data
#' @param col_chr a string, the name of the column containing the chromosome
#' @param col_start a string, the name of the column containing the start position
#' @param col_end a string, the name of the column containing the end position (can be the same as col_start for SNP-based segmentation where start=end)
#' @param col_values a vector of strings, the names of the columns containing the values of interest (logR, BAF, etc.)
#' @return A named list with segments being a data.frame with the summarised information and IDs being a list of SNPs/loci associated with the different segments
#' @examples
#' DF <- data.frame(chr=c(rep("chr1", 10),rep("chr2", 6)),
#'                  pos=c(1:10*1e3, 1:6*1e3),
#'                  logR=c(rep(0, 4), rep(0.54, 3), rep(0, 3), rep(-0.86, 3), rep(0, 3)),
#'                  BAF=c(rep(0.5, 4), rep(0.34, 3), rep(0.5, 3), rep(0.09, 3), rep(0.5, 3)),
#'                  row.names=paste0("SNP_", 1:16))
#' DF
#' summarise_segmetation(DF, "chr", "pos", "pos", c("logR", "BAF"))
#' @author tlesluyes
#' @export
summarise_segmetation <- function(DF, col_chr, col_start, col_end, col_values) {
  tmp <- data.frame(chr=DF[, col_chr], start=DF[, col_start], end=DF[, col_end], row.names=rownames(DF))
  tmp <- cbind(tmp, values=DF[, col_values])
  colnames(tmp)[4:ncol(tmp)] <- col_values
  myRLE <- rle(apply(tmp[, c(1, 4:ncol(tmp))], 1, paste0, collapse="_"))
  tmp$segment <- NA
  for (i in 1:length(myRLE$lengths)) {
    if (i==1) {
      tmp$segment[1:myRLE$lengths[1]] <- 1
    } else {
      tmp$segment[(sum(myRLE$lengths[1:(i-1)])+1):sum(myRLE$lengths[1:i])] <- i
    }
  }; rm(i)
  tmp <- split(tmp, tmp$segment)
  IDs <- lapply(tmp, rownames)
  tmp <- do.call(rbind, lapply(tmp, function(x) {
    tmp2 <- data.frame(chr=x$chr[1],
                       start=min(x$start),
                       end=max(x$end),
                       markers=nrow(x))
    for (col_value in col_values) {
      tmp2[, col_value] <- x[1, col_value]
    }; rm(col_value)
    tmp2$segment <- x$segment[1]
    rownames(tmp2) <- x$segment[1]
    return(tmp2)
  }))
  return(list(segments=tmp, IDs=IDs))
}
