# myFun
myFun is a collection of my favorite R functions, packaged for simplicity.

[![R package version](https://img.shields.io/github/r-package/v/tlesluyes/myFun?color=blue)](manual/myFun.pdf) [![License](https://img.shields.io/github/license/tlesluyes/myFun?label=License&color=blue)](LICENSE) [![GitHub workflow status](https://img.shields.io/github/actions/workflow/status/tlesluyes/myFun/r.yml?logo=github&label=R%20CMD%20check)](https://github.com/tlesluyes/myFun/actions/)

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

## Some examples
What are the logR and BAF values for a 2+1 segment in a diploid tumour with 70% purity?
```R
computeLogR(2, 1, 0.70, 2)
computeBAF(2, 1, 0.70, 2)
```

What are the chromosome sizes in hg38?
```R
load_CHRsize("hg38")
```

How to get the genomic locations covered in different CNA profiles?
```R
require("GenomicRanges")
GR1=GRanges(seqnames="1", ranges=IRanges(start=1, end=1000), nMajor=1, nMinor=1)
GR2=GRanges(seqnames="1", ranges=IRanges(start=10, end=2000), nMajor=2, nMinor=1)
GR3=GRanges(seqnames="1", ranges=IRanges(start=c(5, 21), end=c(20, 2000)), nMajor=c(1, 2), nMinor=c(1, 0))
harmonizeGRanges(list(GR1, GR2, GR3))
```

## Manual
All of the functions are described in [manual/myFun.pdf](manual/myFun.pdf).