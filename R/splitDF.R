#' @title splitDF
#' @description Split a data.frame
#' @details This function splits a data.frame into a list of data.frames.
#' @param DF a data.frame to split
#' @param chunks a number of chunks to obtain
#' @return A list of data.frames
#' @examples
#' DF=data.frame(a=1:26, b=letters)
#' splitDF(DF, 3)
#' @author tlesluyes
#' @export
splitDF=function(DF, chunks) {
  return(split(DF, factor(sort(1:nrow(DF) %% chunks))))
}