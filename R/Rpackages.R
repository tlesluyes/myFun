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
  Bioconductor=as.data.frame(html_table(read_html(url))[[1]])[, c("Package"), drop=FALSE]
  close(url)
  Bioconductor$Source="Bioconductor"
  # Pick installed packages
  myPackages=as.data.frame(installed.packages())
  # Determine the source of each installed package
  myPackages$Source=foreach(i=myPackages$Package, .combine=c) %do% {
    i=get("i")
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

#' @title cleanStringDeps
#' @description This function cleans the string of dependencies
#' @param x a string to be cleaned
#' @return A cleaned string
#' @author tlesluyes
#' @noRd
cleanStringDeps=function(x) {
  x=gsub("\\\n", "", x)
  x=gsub("\\( ?>=? ?[0-9.\\-]+ ?\\)", "", x)
  x=gsub("^ ", "", x)
  x=gsub(" $", "", x)
  x=as.character(x)
  return(x)
}

#' @title myGroupsDeps
#' @description This function defines groups based on the number of dependencies.
#' @param x a numeric vector
#' @return A vector of groups. "No" means no dependencies, "Few" means less than 10 dependencies and "Lot" means at least 10 dependencies
#' @author tlesluyes
#' @noRd
myGroupsDeps=function(x) {
  out=rep("TBD", length(x))
  out[x==0]="No"
  out[x>0 & x<10]="Few"
  out[x>=10]="Lot"
  return(out)
}

#' @title RpackageDependencies
#' @description Show the package dependencies
#' @details Given a folder of R packages, this function reads the DESCRIPTION files of the installed packages and shows their dependencies.
#' @param customFolder a vector of folder names (default: NULL; .libPaths() is used)
#' @param customDependencyTypes a vector of dependency types, possible values are: "Depends", "Imports", "LinkingTo", "Suggests" and "Enhances" (default: c("Depends", "Imports", "LinkingTo"))
#' @param customColours a named vector of colours. Names must correspond to the dependency types and values must be valid colours (default: NULL; an internal colour scheme is used)
#' @param simplifyNetwork a boolean defining if the network should be simplified, i.e. the R base packages are removed (default: TRUE)
#' @param saveFile a string defining the name of the HTML file where the network should be saved (default: NULL; no file is saved)
#' @return A list with nodes (a data.frame of R packages), links (a data.frame of package dependencies) and plot (a network plot using networkD3)
#' @examples
#' myDep=RpackageDependencies()
#' print(head(myDep$nodes))
#' print(head(myDep$links))
#' @author tlesluyes
#' @export
RpackageDependencies=function(customFolder=NULL, customDependencyTypes=NULL, customColours=NULL, simplifyNetwork=TRUE, saveFile=NULL) {
  # Check input parameters
  if (is.null(customFolder)) {
    customFolder=.libPaths()
  } else {
    stopifnot(all(dir.exists(customFolder)))
  }
  customFolder=dir(customFolder, full.names=TRUE)
  customFolder=customFolder[which(!duplicated(basename(customFolder)))]
  names(customFolder)=basename(customFolder)
  customFolder=customFolder[file.exists(paste0(customFolder, "/DESCRIPTION"))]
  stopifnot(length(customFolder)>0)

  if (is.null(customDependencyTypes)) {
    customDependencyTypes=c("Depends", "Imports", "LinkingTo")
  } else {
    stopifnot(all(customDependencyTypes %in% c("Depends", "Imports", "LinkingTo", "Suggests", "Enhances")))
  }

  if (is.null(customColours)) {
    customColours=c("Depends"="#9e9e9e33",
                    "Imports"="#9e9e9e33",
                    "LinkingTo"="#df536b33",
                    "Suggests"="#f5c71033",
                    "Enhances"="#61d04f33")
  } else {
    stopifnot(all(names(customColours) %in% customDependencyTypes))
  }

  stopifnot(length(simplifyNetwork)==1 && is.logical(simplifyNetwork))

  PACKAGES=foreach(DIR=customFolder, .final=function(x) setNames(x, names(customFolder))) %do% {
    DIR=get("DIR")
    # Read the DESCRIPTION file
    DESC=read.dcf(paste0(DIR, "/DESCRIPTION"))
    OUT=list()
    # For each dependency type
    for (i in customDependencyTypes) {
      if (i %in% colnames(DESC)) {
        # Get the cleaned list of dependencies
        OUT[[i]]=setdiff(strsplit(cleanStringDeps(DESC[, i]), " ?, ?")[[1]], "R")
        if (length(OUT[[i]])==0) OUT[[i]]=NA
      } else {
        OUT[[i]]=NA
      }
    }; rm(i)
    return(OUT)
  }

  network=list()
  # Create a data.frame of links where source is the package and target is the dependency
  network$links=foreach(INDEX=1:length(PACKAGES), .combine=rbind) %do% {
    INDEX=get("INDEX")
    tmp=data.frame(source="TBD", target="TBD", type="TBD")
    for (i in names(PACKAGES[[INDEX]])) {
      if (length(PACKAGES[[INDEX]][[i]])==1 && is.na(PACKAGES[[INDEX]][[i]])) next
      tmp=rbind(tmp, data.frame(source=names(PACKAGES)[INDEX], target=PACKAGES[[INDEX]][[i]], type=i))
    }; rm(i)
    if (nrow(tmp)==1) return(NULL)
    tmp=tmp[-1, ]
    return(tmp)
  }
  # Create a data.frame of nodes
  network$nodes=data.frame(name=sort(union(names(customFolder), unlist(lapply(PACKAGES, unlist)))))
  if (simplifyNetwork) {
    BASE=c("base", "compiler", "datasets", "grDevices", "graphics", "grid", "methods", "parallel", "splines", "stats", "stats4", "tcltk", "tools", "translations", "utils")
    network$links=network$links[-which(network$links$source %in% BASE | network$links$target %in% BASE), ]
    network$nodes=network$nodes[-which(network$nodes$name %in% BASE), , drop=FALSE]
    rm(BASE)
  }

  # Get the number of targets and sources for all packages
  network$nodes$ntarget=sapply(network$nodes$name, function(x) length(which(network$links$target==x)))
  network$nodes$nsource=sapply(network$nodes$name, function(x) length(which(network$links$source==x)))
  network$nodes$color=myGroupsDeps(network$nodes$nsource)

  # Get the node indexes for the links
  CORR=data.frame(INDEX=1:nrow(network$nodes), row.names=network$nodes$name)
  network$links$source_index=CORR[network$links$source, "INDEX"]-1 # Indexes are 0-based
  network$links$target_index=CORR[network$links$target, "INDEX"]-1 # Indexes are 0-based
  rm(CORR)
  colourScale='d3.scaleOrdinal().domain(["No", "Few", "Lot"]).range(["#000000", "#0000ff", "#ff0000"]);'
  network$links$value=1
  # Create the plot
  network$plot=forceNetwork(Links=network$links,
                            Nodes=network$nodes,
                            Source="source_index",
                            Target="target_index",
                            NodeID="name",
                            Nodesize="ntarget",
                            Group="color",
                            zoom=TRUE,
                            arrows=TRUE,
                            Value="value",
                            linkColour=customColours[network$links$type],
                            opacity=1,
                            colourScale=colourScale,
                            radiusCalculation="Math.sqrt(d.nodesize*2)+3")
  # Save the plot if needed
  if (!is.null(saveFile)) {
    saveNetwork(network$plot, saveFile)
  }
  # Get the number of free nodes
  FREE=which(network$nodes$ntarget==0 & network$nodes$nsource==0)
  if (length(FREE)>0) {
    print(paste0(length(FREE), " free nodes: ", paste0(network$nodes$name[FREE], collapse=", ")))
  }
  # Return the final object
  invisible(network)
}