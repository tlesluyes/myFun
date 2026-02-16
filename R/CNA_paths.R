#' @title all_paths
#' @description Get all possibles paths between two copy-number states
#' @param start a vector of length 2 (representing a copy-number state; e.g. c(1, 1) represents a 1+1 state), defining where to start
#' @param end a vector of length 2 (representing a copy-number state; e.g. c(1, 1) represents a 1+1 state), defining where to end
#' @param WGD a boolean defining if WGD events are allowed
#' @param max_path_size an integer defining the maximum path size
#' @param simplify a boolean defining if consecutive and opposite alterations (e.g. +1/+0 and then -1/-0) are allowed
#' @param path a string defining the current path
#' @param path_size an integer defining the current path size
#' @param OUT a vector of possible paths
#' @return A vector of all possible paths given as characters (separator=";")
#' @author tlesluyes
#' @noRd
all_paths <- function(start, end, WGD, max_path_size, simplify, path, path_size, OUT) {
  # If the path is valid, store it
  if (identical(start, end)) OUT <- c(OUT, path)
  # Prevent infinite recursive calls
  if (path_size>max_path_size) return(OUT)
  # Get the alterations since the last WGD
  lastAlterations <- strsplit(path, ";")[[1]]
  if ("WGD" %in% lastAlterations) lastAlterations <- lastAlterations[(which.max(lastAlterations=="WGD")+1):length(lastAlterations)]
  # This is the major allele, allow alterations if there is no HD (>0)
  if (start[1]>0) {
    if (!simplify || ! "-1/-0" %in% lastAlterations) OUT <- all_paths(c(start[1]+1, start[2]), end, WGD, max_path_size, simplify, paste(path, "+1/+0", sep=";"), path_size+1, OUT)
    if (start[1]-1>=start[2] && (!simplify || ! "+1/+0" %in% lastAlterations)) OUT <- all_paths(c(start[1]-1, start[2]), end, WGD, max_path_size, simplify, paste(path, "-1/-0", sep=";"), path_size+1, OUT)
    # This is the minor allele, allow alterations if there is no LOH (>0)
    if (start[2]>0) {
      if (start[1]>=start[2]+1 && (!simplify || ! "-0/-1" %in% lastAlterations)) OUT <- all_paths(c(start[1], start[2]+1), end, WGD, max_path_size, simplify, paste(path, "+0/+1", sep=";"), path_size+1, OUT)
      if (!simplify || ! "+0/+1" %in% lastAlterations) OUT <- all_paths(c(start[1], start[2]-1), end, WGD, max_path_size, simplify, paste(path, "-0/-1", sep=";"), path_size+1, OUT)
    }
  }
  # If WGD is allowed, try it
  if (WGD) OUT <- all_paths(start*2, end, WGD, max_path_size, simplify, paste(path, "WGD", sep=";"), path_size+1, OUT)
  # Return the vector of valid paths
  return(OUT)
}

#' @title get_all_paths
#' @description Get all possibles paths between two copy-number states
#' @details This function returns all possible paths between two copy-number states. The expected input is allele-specific (with two values), but it can be used for total copy-number by setting c(ntot, 0). Possible events include: +1/+0 (gain of the major allele), -1/-0 (loss of the major allele), +0/+1 (gain of the minor allele), -0/-1 (loss of the minor allele) and WGD.
#' @param start a vector of length 2 (representing a copy-number state; e.g. c(1, 1) represents a 1+1 state), defining where to start
#' @param end a vector of length 2 (representing a copy-number state; e.g. c(1, 1) represents a 1+1 state), defining where to end
#' @param WGD a boolean defining if WGD events are allowed
#' @param max_path_size an integer defining the maximum path size
#' @param simplify a boolean defining if consecutive and opposite alterations (e.g. +1/+0 and then -1/-0) are allowed
#' @return A vector of all possible paths given as characters (separator=";")
#' @author tlesluyes
#' @examples
#' # Diploid baseline (1+1) turns into 2+1
#' print(get_all_paths(start=c(1, 1), end=c(2, 1), WGD=TRUE))
#' # Chromosome X in males (1+0) is gained (5 copies)
#' print(get_all_paths(start=c(1, 0), end=c(5, 0), WGD=TRUE))
#' @export
get_all_paths <- function(start, end, WGD, max_path_size=5, simplify=TRUE) {
  # Check input parameters
  stopifnot(length(start)==2 && all(is.numeric(start)))
  stopifnot(start[1]>=start[2])
  stopifnot(length(end)==2 && all(is.numeric(end)))
  stopifnot(end[1]>=end[2])
  stopifnot(length(max_path_size)==1 && is.numeric(max_path_size))
  stopifnot(length(WGD)==1 && is.logical(WGD))
  stopifnot(length(simplify)==1 && is.logical(simplify))
  stopifnot(start[1]>0 || (start[1]==0 && end[1]==0))
  stopifnot(start[2]>0 || (start[2]==0 && end[2]==0))
  # Run the recursive function
  paths <- all_paths(start=start,
                     end=end,
                     WGD=WGD,
                     max_path_size=max_path_size,
                     simplify=simplify,
                     path=c("Start"),
                     path_size=1,
                     OUT=c())
  # Throw an error if no path is found
  if (length(paths)==0) {
    stop("Error: No valid paths found, try increasing max_path_size")
  }
  # Reformat the output
  stopifnot(all(sapply(strsplit(paths, ";"), "[", 1)=="Start"))
  paths <- gsub("^Start;?", "", paths)
  return(paths)
}

#' @title get_shortest_path
#' @description Get the shortest path among several
#' @details This function returns the shortest possible path. It should be used after running the get_all_paths function or can be used as long as the input format is correct.
#' @param paths all possible paths to consider
#' @param wanted_WGD a numeric value defining the number of WGD events wanted (can be NA to allow for any possibility, including no event at all; default: NA)
#' @param count_WGD a boolean defining if the number of WGD events should be counted (default: FALSE)
#' @return A numeric value representing the minimal number of events, its name represents the full path
#' @author tlesluyes
#' @examples
#' # Diploid baseline (1+1) turns into 2+1
#' print(get_shortest_path(get_all_paths(start=c(1, 1), end=c(2, 1), WGD=TRUE)))
#' # Chromosome X in males (1+0) is gained (5 copies)
#' print(get_shortest_path(get_all_paths(start=c(1, 0), end=c(5, 0), WGD=TRUE)))
#' @export
get_shortest_path <- function(paths, wanted_WGD=NA, count_WGD=FALSE) {
  # Check input parameters
  stopifnot(length(paths)>0 && all(is.character(paths)))
  stopifnot(length(wanted_WGD)==1 && (is.na(wanted_WGD) || is.numeric(wanted_WGD)))
  stopifnot(length(count_WGD)==1 && is.logical(count_WGD))
  # Check paths
  paths <- strsplit(paths, ";")
  stopifnot(all(unlist(paths) %in% c("+1/+0", "+0/+1", "-1/-0", "-0/-1", "WGD")))
  # If a number of WGD is wanted, go for it
  if (!is.na(wanted_WGD)) {
    n_WGD <- sapply(paths, function(x) length(which(x=="WGD")))
    stopifnot(wanted_WGD %in% n_WGD)
    paths <- paths[which(n_WGD==wanted_WGD)]
  }
  path_size <- sapply(paths, length)
  if (!count_WGD) path_size <- path_size-sapply(paths, function(x) length(which(x=="WGD")))
  names(path_size) <- sapply(paths, paste, collapse=";")
  # Keep the shortest path
  path_size <- path_size[which.min(path_size)]
  return(path_size)
}
