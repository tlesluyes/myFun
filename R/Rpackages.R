#' @title Rpackages
#' @description List installed packages and determine their source
#' @details This function lists installed packages and determine whether they are base packages or come from CRAN/Bioconductor or if they are external (GitHub, SourceForge, etc.). This function requires an internet connection.
#' @param CRAN_URL the CRAN URL (default: "http://cran.us.r-project.org")
#' @param Bioconductor_URL the Bioconductor URL (default: "https://www.bioconductor.org/packages/release/bioc/")
#' @return A data.frame with the installed packages and an additional column: Source (possible values: Base, CRAN, Bioconductor, External)
#' @examples head(Rpackages())
#' @author tlesluyes
#' @export
Rpackages=function(CRAN_URL="http://cran.us.r-project.org",
                   Bioconductor_URL="https://www.bioconductor.org/packages/release/bioc/") {
  # Pick CRAN packages
  CRAN=as.data.frame(available.packages(repos=CRAN_URL))[, c("Package"), drop=FALSE]
  CRAN$Source="CRAN"
  # Pick Bioconductor packages
  url=url(Bioconductor_URL, "rb")
  Bioconductor=as.data.frame(rvest::html_table(rvest::read_html(url))[[1]])[, c("Package"), drop=FALSE]
  close(url)
  Bioconductor$Source="Bioconductor"
  # Pick installed packages
  myPackages=as.data.frame(installed.packages())
  # Determine the source of each installed package
  myPackages$Source=foreach::foreach(i=myPackages$Package, .combine=c) %do% {
    # Test whether it is from CRAN and/or Bioconductor
    FROM=c(ifelse(i %in% CRAN$Package, "CRAN", NA),
           ifelse(i %in% Bioconductor$Package, "Bioconductor", NA))
    FROM=FROM[!is.na(FROM)]
    # If it is not from CRAN or Bioconductor, test whether it is a base package
    if (length(FROM)==0 && i %in% c("base", "compiler", "datasets", "grDevices", "graphics", "grid", "methods", "parallel", "splines", "stats", "stats4", "tcltk", "tools", "translations", "utils")) FROM="Base"
    # Otherwise, it looks like an external package
    if (length(FROM)==0) FROM="External"
    # The collapse is needed if a package exists in both CRAN and Bioconductor (not sure this could happen though)
    return(paste0(FROM, collapse="/"))
  }
  return(myPackages)
}