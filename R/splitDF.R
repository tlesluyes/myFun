#' @title splitDF
#' @description Split a data.frame
#' @details This function splits a data.frame into a list of data.frames.
#' @param DF a data.frame to split
#' @param chunks a number of chunks to obtain
#' @param shuffle a boolean, whether to shuffle the data.frame before splitting (default: FALSE)
#' @param seed a number, the seed for the random number generator (default: 1234)
#' @return A list of data.frames
#' @examples
#' DF <- data.frame(a=1:26, b=letters)
#' splitDF(DF, 3)
#' @author tlesluyes
#' @export
splitDF <- function(DF, chunks, shuffle=FALSE, seed=1234) {
  stopifnot(length(chunks)==1 && all.equal(chunks, as.integer(chunks)))
  stopifnot(is.data.frame(DF) && nrow(DF)>chunks)
  stopifnot(length(shuffle)==1 && is.logical(shuffle))
  stopifnot(length(seed)==1 && all.equal(seed, as.integer(seed)))
  if (shuffle) {
    set.seed(seed)
    DF <- DF[sample(1:nrow(DF), nrow(DF)), ]
  }
  return(split(DF, factor(sort(1:nrow(DF) %% chunks))))
}
