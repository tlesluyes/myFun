# myFun
myFun is a collection of my favorite R functions, packaged for simplicity.

![R package version](https://img.shields.io/github/r-package/v/tlesluyes/myFun?color=blue) ![License](https://img.shields.io/github/license/tlesluyes/myFun?label=License&color=blue) ![GitHub workflow status](https://img.shields.io/github/actions/workflow/status/tlesluyes/myFun/r.yml?logo=github&label=R%20CMD%20check)

## Installation

### Dependencies
```R
install.packages(c("doParallel", "foreach"))
BiocManager::install(c("GenomicRanges", "IRanges"))
```

### myFun
```R
devtools::install_github("tlesluyes/myFun")
```