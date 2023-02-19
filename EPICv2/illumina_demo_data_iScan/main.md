ChAMP run on illumina EPICv2 iScan Demo Data
================
Yuan Tian (<champ450K@gmail.com>)
2023-02-19 00:12:51

This is a run markdown for ChAMP supporting on illumina HumanMethylation
EPICv2 array. To replicate this run, you need to install the ChAMPdata
(\>=2.31.1), and ChAMP version \>= 2.22.1. The latest ChAMPdata can be
obtained [here](https://github.com/YuanTian1991/ChAMPdata). And the
latest ChAMP is [here](https://github.com/YuanTian1991/ChAMP).

In this demo, I am using the defatul [demo data provided by
illumina](https://emea.support.illumina.com/array/array_kits/infinium-methylationepic-beadchip-kit/downloads.html).

``` r
library("ChAMP")
```

After the above library loading, you can check library version with
`sessionInfo()`.

## 1. Import IDAT Data

Like previous version, ChAMP requires only one CSV to be put in a
folder, along with IDAT files to be loaded.

``` r
library("illuminaio")

source("~/Personal/ChAMP-Dev/ChAMP/R/champ.import.R")

myImport <- champ.import("./DemoDataEPICv2/", arraytype = "EPICv2") # This can also be set just as `EPIC`
```

## 2. Filtering

After loading, filtering will be done, specifically, the SNP mask is
using the one provided by Zhou
[here](http://zwdzwd.github.io/InfiniumAnnotation).

``` r
library("impute")

source("~/Personal/ChAMP-Dev/ChAMP/R/champ.filter.R")
myFilter <- champ.filter(beta=myImport$beta, 
                         pd=myImport$pd, 
                         detP=myImport$detP,
                         beadcount=myImport$beadcount, 
                         ProbeCutoff=0.1,
                         arraytype = "EPICv2") # This can also be set just as `EPIC`
```

Below example shows if you need to do population-specific filtering.
Currently, EPICv2 annotation only support AFR, AMR, EAS, EUR and SAS.

``` r
myFilter <- champ.filter(beta=myImport$beta, 
                         pd=myImport$pd, 
                         detP=myImport$detP,
                         beadcount=myImport$beadcount,
                         population="EUR",
                         arraytype = "EPICv2") # This can also be set just as `EPIC`
```

## 3. Loading

Above two step can be merged into just one `champ.load()`. The defaul
champ.load have two method: `ChAMP` or `minfi`, currently minfi package
is not working for EPICv2, so only `ChAMP` works here.

``` r
source("~/Personal/ChAMP-Dev/ChAMP/R/champ.load.R")

myLoad <- champ.load("./DemoDataEPICv2/", arraytype = "EPICv2")
```

## 4. CpG.GUI

Then, we can use CpG.GUI() to explore a beta matrix, along with a vector
of phenotype.

``` r
source("~/Personal/ChAMP-Dev/ChAMP/R/CpG.GUI.R")

CpG.GUI()
```
