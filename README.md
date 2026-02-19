# myFun

myFun is a collection of my favourite R functions, packaged for simplicity.

[![R package version](https://img.shields.io/github/r-package/v/tlesluyes/myFun?color=blue)](manual/myFun.pdf) [![License](https://img.shields.io/github/license/tlesluyes/myFun?label=License&color=blue)](LICENSE) [![GitHub workflow status](https://img.shields.io/github/actions/workflow/status/tlesluyes/myFun/r.yml?logo=github&label=R%20CMD%20check)](https://github.com/tlesluyes/myFun/actions/)

## Installation

```R
# General dependencies
install.packages(c("devtools", "BiocManager"))
# Package dependencies
install.packages(c("doParallel", "foreach", "networkD3", "rvest"))
BiocManager::install(c("GenomeInfoDb", "GenomicRanges", "IRanges", "S4Vectors"))
# myFun
devtools::install_github("tlesluyes/myFun")
```

## Manual

All of the functions are described in [manual/myFun.pdf](manual/myFun.pdf).

## Some examples

### LogR, BAF, purity and ploidy values

What are the logR and BAF values of a 2+1 segment in a diploid tumour with 70% purity?

```R
computeLogR(2, 1, 0.70, 2)
# 0.433
computeBAF(2, 1, 0.70, 2)
# 0.37 0.63
```

A genomic segment has logR=0.5361 and BAF=0.3448/0.6552. Since we expect a 2+1 copy-number status, what would be the purity/ploidy fit for it?

```R
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

A pseudo-diploid sample (ploidy=2.4) has purity=74%. What is the re-estimated ploidy if I believe that the sample has purity=61% instead (e.g. clonal peak of CCF estimates from SNVs)?

```R
reestimate_ploidy(0.74, 2.4, 0.61, 0)
# 2.4852
```

A sample has purity=74% and ploidy=2.4 but the CNA profile needs to be doubled. What is the re-estimated purity?

```R
reestimate_purity(0.74, 2.4, "double")
# 0.5873
```

## Analysing CNA profiles

What are the chromosome sizes and bands in hg38?

```R
load_CHRsize("hg38")
head(CHRsize, n=3)
# chr      size     middle        sum        add
#   1 248956422  124478211  248956422          0
#   2 242193529  370053186  491149951  248956422
#   3 198295559  590297730  689445510  491149951

load_cytoband("hg38")
head(cytoband, n=3)
# chr    start      end   name gieStain  color start_adj  end_adj
#   1        1  2300000 p36.33     gneg grey90         1  2300000
#   1  2300001  5300000 p36.32   gpos25 grey80   2300001  5300000
#   1  5300001  7100000 p36.31     gneg grey90   5300001  7100000
```

Let's generate a playground with 5 random CNA profiles.

```R
CNAprofiles <- example_CNAs(n=5, assembly="hg38")
head(CNAprofiles[[1]], n=3)
# GRanges object with 3 ranges and 5 metadata columns:
#       seqnames      ranges strand |    nMajor    nMinor      gain      loss    CNstatus
#          <Rle>   <IRanges>  <Rle> | <numeric> <numeric> <logical> <logical> <character>
#   [1]        1 1-248956422      * |         1         1     FALSE     FALSE         1+1
#   [2]        2 1-242193529      * |         2         1      TRUE     FALSE         2+1
#   [3]        3 1-198295559      * |         1         1     FALSE     FALSE         1+1
#   -------
#   seqinfo: 24 sequences from hg38 genome
```

In the above example, `CNAprofiles` is a list of 5 `GRanges` objects. `seqinfo` for these objects has been generated already and can be obtained using `GRanges(..., seqinfo=get_Seqinfo("hg38"))`.

Also, genomic regions have been harmonized using the `harmonizeGRanges` function. This ensure that `granges` are the same across all samples so they are 1:1 matches. By default, it will remove any region which is not covered in all profiles but this behaviour can be changed using `keepHoles=TRUE` so metadata in missing regions will be set to `NA`.

From there, we may exclude genomic regions such as centromeres.

```R
centromeres <- get_centromeres("hg38")
head(centromeres, n=3)
# chr     start       end
#   1 121700001 125100000
#   2  91800001  96000000
#   3  87800001  94000000
centromeres <- makeGRangesFromDataFrame(centromeres)
CNAprofiles <- lapply(CNAprofiles, excludeGRanges, centromeres)
head(CNAprofiles[[1]], n=3)
# GRanges object with 3 ranges and 5 metadata columns:
# seqnames              ranges strand |    nMajor    nMinor      gain      loss    CNstatus
#    <Rle>           <IRanges>  <Rle> | <numeric> <numeric> <logical> <logical> <character>
#        1         1-121700000      * |         1         1     FALSE     FALSE         1+1
#        1 125100001-248956422      * |         1         1     FALSE     FALSE         1+1
#        2          1-91800000      * |         2         1      TRUE     FALSE         2+1
# ---
# seqinfo: 24 sequences from hg38 genome
```

This package offers two ways of measuring distances between CNA profiles.

1. Inter-sample agreement (ISA): the percentage of the genome with identical copy-number states. It's a similarity metrics.

2. Manhattan distance (MD): the allele-specific distance between two profiles. It's a distance metrics.

```R
# Compute metrics on the first two profiles
computeISA(CNAprofiles[[1]], CNAprofiles[[2]], CNstatus="CNstatus") # Inter-sample agreement (%)
# 51.65065
computeMD(CNAprofiles[[1]], CNAprofiles[[2]], nMajor="nMajor", nMinor="nMinor", convertMb=TRUE) # Manhattan distance (Mb)
# 1829.669

# Compute metrics on all profiles
computeISA_batch(CNAprofiles, CNstatus="CNstatus")
#           G1        G2        G3        G4        G5
# G1 100.00000  51.65065  40.86227  64.99530  48.34231
# G2  51.65065 100.00000  53.62919  57.94543  52.52448
# G3  40.86227  53.62919 100.00000  64.55561  69.43959
# G4  64.99530  57.94543  64.55561 100.00000  76.68070
# G5  48.34231  52.52448  69.43959  76.68070 100.00000
computeMD_batch(CNAprofiles, nMajor="nMajor", nMinor="nMinor", convertMb=TRUE)
#          G1       G2        G3        G4        G5
# G1    0.000 1829.669 1747.8111 1131.3126 1533.1636
# G2 1829.669    0.000 1612.6254 1214.7242 1323.6658
# G3 1747.811 1612.625    0.0000  988.2257  944.9499
# G4 1131.313 1214.724  988.2257    0.0000  773.5782
# G5 1533.164 1323.666  944.9499  773.5782    0.0000
```

**Warning**: because of the way the profiles are harmonized, computing these metrics on two profiles versus in batch mode will likely give slightly different results. This is because harmonizing all profiles altogether will remove any regions missing in at least one sample, although such regions may be considered when only comparing two profiles. It is therefore recommended to compute these metrics on two profiles (so only minimal regions are excluded) and not in batch mode (where more regions are excluded). However, when comparing a large number of profiles (e.g. >1,000), this process can be very time consuming and the batch mode can be used as an approximation in a much faster way. In the above example, the CNA profiles don't have any sample-specific gaps so the results are identical.

From there, we can calculate the fraction of gains and losses in the entire set of CNA profiles

```R
gains_losses <- occurrenceGRanges(CNprofiles, c("gain","loss"))
head(gains_losses, n=3)
# GRanges object with 3 ranges and 5 metadata columns:
# seqnames              ranges strand |  nSamples      gain    gain_p      loss    loss_p
#    <Rle>           <IRanges>  <Rle> | <integer> <integer> <numeric> <integer> <numeric>
#        1         1-121700000      * |         5         0         0         0         0
#        1 125100001-248956422      * |         5         0         0         0         0
#        2          1-91800000      * |         5         1        20         1        20
# ---
# seqinfo: 24 sequences from hg38 genome
```

We can now plot the landscape of gains and losses. Let's first get adjusted genomic positions in such a way that all chromosomes are adjacent to each other. Then, such adjusted position will be used in the plot.

```R
gains_losses <- data.frame(gains_losses)
gains_losses <- adjustPositions(gains_losses, CHRsize, chr_column="seqnames")
head(gains_losses[, c("seqnames", "start", "end", "start_adj", "end_adj")], n=3)
# seqnames     start       end start_adj   end_adj
#        1         1 121700000         1 121700000
#        1 125100001 248956422 125100001 248956422
#        2         1  91800000 248956423 340756422
par(mar=c(2.1, 4.2, 0.5, 0.5))
plot(NULL, ylim=c(-100, 100), xlim=c(1, CHRsize$sum[nrow(CHRsize)]), xaxt="n", yaxt="n", xlab="", ylab="", cex.lab=1.25, xaxs="i")
axis(1, at=CHRsize$middle, labels = CHRsize$chr, las=2)
rect(tmp$start_adj, 0, tmp$end_adj, tmp$gain_p, col="indianred1", border=NA)
rect(tmp$start_adj, 0, tmp$end_adj, -tmp$loss_p, col="cornflowerblue", border=NA)
abline(v=CHRsize$sum[1:(nrow(CHRsize)-1)], col="grey")
abline(h=0)
axis(2, at=seq(-100, 100, by=25), labels = paste0(abs(seq(-100, 100, by=25)), "%"), las=3)
mtext(c("Gain", "Loss"), side=2, at=c(50, -50), line=2.5, cex=2)
```

There are other functions of interest for checking `GRanges` objects such as `checkGR` for individual objects and `checkGRlist` for lists of GRanges objects. As for metadata, `cleanGRlistMetadata` will remove any region set to `NA` in a list of `GRanges` objects.

### Tumour evolution

How do I go from a 1+1 copy-number state to a 2+0?

```R
myPaths <- get_all_paths(c(1, 1), c(2, 0), WGD=TRUE)
print(myPaths)
# [1] "+1/+0;-0/-1"                 "+1/+0;-0/-1;WGD;-1/-0;-1/-0"
# [3] "-0/-1;+1/+0"                 "-0/-1;+1/+0;WGD;-1/-0;-1/-0"
# [5] "-0/-1;WGD"                   "-0/-1;WGD;-1/-0;WGD"        
# [7] "-0/-1;WGD;WGD;-1/-0;-1/-0"   "WGD;-0/-1;-1/-0;-0/-1;WGD"  
# [9] "WGD;-0/-1;-0/-1"             "WGD;-0/-1;-0/-1;-1/-0;WGD"
```

What is the shortest path?

```R
get_shortest_path(myPaths) # shortest path: 1 alteration (not including WGD)
# -0/-1;WGD
#         1
```

The above path includes a WGD event. What would be the shortest path without any WGD?

```R
get_shortest_path(myPaths, wanted_WGD=0) # shortest path without any WGD: 2 alterations
# +1/+0;-0/-1
#           2
```

### Other utilities

**How to split a DF?** (useful for parallelization)

```R
DF <- data.frame(a=1:26, b=letters)
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

**How to summrise segmented SNP information?**

```R
DF <- data.frame(chr=c(rep("chr1", 10), rep("chr2", 6)),
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

**Where do my packages come from?**

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
myPackages <- RpackageDependencies()
head(myPackages$nodes)
print(myPackages$plot) # networkD3 plot
```
