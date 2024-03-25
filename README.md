# myFun
myFun is a collection of my favorite R functions, packaged for simplicity.

[![R package version](https://img.shields.io/github/r-package/v/tlesluyes/myFun?color=blue)](manual/myFun.pdf) [![License](https://img.shields.io/github/license/tlesluyes/myFun?label=License&color=blue)](LICENSE) [![GitHub workflow status](https://img.shields.io/github/actions/workflow/status/tlesluyes/myFun/r.yml?logo=github&label=R%20CMD%20check)](https://github.com/tlesluyes/myFun/actions/)

## Installation

### Dependencies
```R
# General dependencies
install.packages(c("devtools", "BiocManager"))
# Package dependencies
install.packages(c("doParallel", "foreach", "rvest", "networkD3"))
BiocManager::install(c("GenomicRanges", "IRanges"))
```

### myFun
```R
devtools::install_github("tlesluyes/myFun")
```
## Manual
All of the functions are described in [manual/myFun.pdf](manual/myFun.pdf).

## Some examples
### What are logR and BAF values for a given copy-number status?
```R
# A 2+1 segment in a diploid tumour with 70% purity

computeLogR(2, 1, 0.70, 2)
# 0.433

computeBAF(2, 1, 0.70, 2)
# 0.37 0.63
```

### What is the purity/ploidy fit if a segment has a specific copy-number status?
```R
# A segment has logR=0.5361 and BAF=0.3448/0.6552, we expect a 2+1 copy-number status
computeFit(0.5361, 0.6552, 2, 1, 1)
# $rho
#[1] 0.9002
#
#$psi
#[1] 2.0001
#
#$psit
#[1] 2.0001
```

### What are re-estimated purity/ploidy values with changes in purity/ploidy?
```R
# A pseudo-diploid sample has purity=74% and ploidy=2.4. What is the re-estimated ploidy if I believe that the sample has purity=61%?
reestimate_ploidy(0.74, 2.4, 0.61, 0)
# 2.4852

# A sample has purity=74% and ploidy=2.4 but the CNA profile needs to be doubled. What is the re-estimated purity?
reestimate_purity(0.74, 2.4, "double")
# 0.5873
```

### What are the chromosome sizes and bands in hg38?
```R
load_CHRsize("hg38")
head(CHRsize)
# chr      size     middle        sum        add
#   1 248956422  124478211  248956422          0
#   2 242193529  370053186  491149951  248956422
#   3 198295559  590297730  689445510  491149951
#   4 190214555  784552788  879660065  689445510
#   5 181538259  970429194 1061198324  879660065
#   6 170805979 1146601314 1232004303 1061198324

load_cytoband("hg38")
head(cytoband)
# chr    start      end   name gieStain  color start_adj  end_adj
#   1        1  2300000 p36.33     gneg grey90         1  2300000
#   1  2300001  5300000 p36.32   gpos25 grey80   2300001  5300000
#   1  5300001  7100000 p36.31     gneg grey90   5300001  7100000
#   1  7100001  9100000 p36.23   gpos25 grey80   7100001  9100000
#   1  9100001 12500000 p36.22     gneg grey90   9100001 12500000
#   1 12500001 15900000 p36.21   gpos50 grey60  12500001 15900000

# Empty plot, CGH-like
plot(NULL, ylim=c(0,1), xlim=c(1, CHRsize$sum[nrow(CHRsize)]), xaxt='n', yaxt="n", xlab="Chromosomes", ylab="Whatever", cex.lab=1.25, xaxs="i")
abline(v=CHRsize$sum[2:nrow(CHRsize)], col="grey")
axis(1, at=CHRsize$middle, labels = CHRsize$chr, las=2)
```

### 

### How to get the genomic locations covered in different CNA profiles?
```R
require("GenomicRanges")
GR1=GRanges(seqnames="1", ranges=IRanges(start=1, end=1000), nMajor=1, nMinor=1)
GR2=GRanges(seqnames="1", ranges=IRanges(start=10, end=2000), nMajor=2, nMinor=1)
GR3=GRanges(seqnames="1", ranges=IRanges(start=c(5, 21), end=c(20, 2000)), nMajor=c(1, 2), nMinor=c(1, 0))
harmonizeGRanges(list(GR1, GR2, GR3))
# [[1]]
# GRanges object with 2 ranges and 3 metadata columns:
#     seqnames    ranges strand |    nMajor    nMinor       hit
#        <Rle> <IRanges>  <Rle> | <numeric> <numeric> <logical>
#   1        1     10-20      * |         1         1      TRUE
#   2        1   21-1000      * |         1         1      TRUE
#   -------
#   seqinfo: 1 sequence from an unspecified genome; no seqlengths
# 
# [[2]]
# GRanges object with 2 ranges and 3 metadata columns:
#     seqnames    ranges strand |    nMajor    nMinor       hit
#        <Rle> <IRanges>  <Rle> | <numeric> <numeric> <logical>
#   1        1     10-20      * |         2         1      TRUE
#   2        1   21-1000      * |         2         1      TRUE
#   -------
#   seqinfo: 1 sequence from an unspecified genome; no seqlengths
# 
# [[3]]
# GRanges object with 2 ranges and 3 metadata columns:
#     seqnames    ranges strand |    nMajor    nMinor       hit
#        <Rle> <IRanges>  <Rle> | <numeric> <numeric> <logical>
#   1        1     10-20      * |         1         1      TRUE
#   2        1   21-1000      * |         2         0      TRUE
#   -------
#   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

### How to measure distances between CNA profiles?
```R
require("GenomicRanges")
GR1=GRanges(seqnames=rep("1", 3),
            ranges=IRanges(start=c(1, 1001, 10001),end=c(1000, 10000, 20000)),
            CNstatus=c("1+1", "2+1", "1+1"),
            nMajor=c(1, 2, 1),
            nMinor=c(1, 1, 1))
GR2=GRanges(seqnames=rep("1", 2),
            ranges=IRanges(start=c(500, 10001),end=c(10000, 25000)),
            CNstatus=c("2+1", "1+1"),
            nMajor=c(2, 1),
            nMinor=c(1, 1))
# in this example:
#    Region 500-1000 (size=501) is 1+1 for GR1 and 2+1 for GR2
#    Region 1001-20000 (size=19000) is identical between GR1 and GR2 (both 2+1 and 1+1)

computeISA(GR1, GR2, CNstatus="CNstatus") # Inter-sample agreement (%)
# 97.4309

computeMD(GR1, GR2, nMajor="nMajor", nMinor="nMinor") # Manhattan distance (bp)
# 501
```

### How do I go from a 1+1 copy-number state to a 2+0?
```R
myPaths=get_all_paths(c(1, 1), c(2, 0), WGD=TRUE)
print(myPaths)
# [1] "+1/+0;-0/-1"                 "+1/+0;-0/-1;WGD;-1/-0;-1/-0"
# [3] "-0/-1;+1/+0"                 "-0/-1;+1/+0;WGD;-1/-0;-1/-0"
# [5] "-0/-1;WGD"                   "-0/-1;WGD;-1/-0;WGD"        
# [7] "-0/-1;WGD;WGD;-1/-0;-1/-0"   "WGD;-0/-1;-1/-0;-0/-1;WGD"  
# [9] "WGD;-0/-1;-0/-1"             "WGD;-0/-1;-0/-1;-1/-0;WGD"

get_shortest_path(myPaths) # shortest path: 1 alteration (not including WGD)
# -0/-1;WGD
#         1

get_shortest_path(myPaths, wanted_WGD=0) # shortest path without any WGD: 2 alterations
# +1/+0;-0/-1
#           2
```

### How to split a DF?
```R
DF=data.frame(a=1:26, b=letters)
splitDF(DF, 3)
# $`0`
#   a b
# 1 1 a
# 2 2 b
# 3 3 c
# ...
# 
# $`1`
#     a b
# 9   9 i
# 10 10 j
# 11 11 k
# ...
# 
# $`2`
#     a b
# 18 18 r
# 19 19 s
# 20 20 t
# ...
```

### How to summrise segmented SNP information?
```R
DF=data.frame(chr=c(rep("chr1", 10),rep("chr2", 6)),
              pos=c(1:10*1e3, 1:6*1e3),
              logR=c(rep(0, 4), rep(0.54, 3), rep(0, 3), rep(-0.86, 3), rep(0, 3)),
              BAF=c(rep(0.5, 4), rep(0.34, 3), rep(0.5, 3), rep(0.09, 3), rep(0.5, 3)),
              row.names=paste0("SNP_", 1:16))
summarise_segmetation(DF, "chr", "pos", "pos", c("logR", "BAF"))
# $segments
#  chr start   end markers  logR  BAF segment
# chr1  1000  4000       4  0.00 0.50       1
# chr1  5000  7000       3  0.54 0.34       2
# chr1  8000 10000       3  0.00 0.50       3
# chr2  1000  3000       3 -0.86 0.09       4
# chr2  4000  6000       3  0.00 0.50       5
# 
# $IDs
# $IDs$`1`
# [1] "SNP_1" "SNP_2" "SNP_3" "SNP_4"
# 
# $IDs$`2`
# [1] "SNP_5" "SNP_6" "SNP_7"
# 
# $IDs$`3`
# [1] "SNP_8"  "SNP_9"  "SNP_10"
# 
# $IDs$`4`
# [1] "SNP_11" "SNP_12" "SNP_13"
# 
# $IDs$`5`
# [1] "SNP_14" "SNP_15" "SNP_16"
```

### Where do my packages come from?
```R
Rpackages()
#               Package   ...         Source
# abind           abind   ...           CRAN
# aCGH             aCGH   ...   Bioconductor
# affxparser affxparser   ...   Bioconductor
# affy             affy   ...   Bioconductor
# affyio         affyio   ...   Bioconductor
# anndata       anndata   ...           CRAN
# ...

# Show package dependencies
myPackages=RpackageDependencies()
head(myPackages$nodes)
print(myPackages$plot) # networkD3 plot
```