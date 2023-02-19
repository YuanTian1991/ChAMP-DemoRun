illumina EPICv2 iScan Demo \[Debug\]
================
Yuan Tian (<champ450K@gmail.com>)
19 February, 2023

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
library("shiny")
library("shinythemes")
library("plotly")

source("~/Personal/ChAMP-Dev/ChAMP/R/CpG.GUI.R")

CpG.GUI(rownames(myLoad$beta), arraytype = "EPICv2")
```

## 5. QC

QC.GUI() can be used for check beta matrix

``` r
library("shiny")
library("shinythemes")
library("plotly")
library("RColorBrewer")
library("dendextend")
library(isva)

source("~/Personal/ChAMP-Dev/ChAMP/R/EstDimRMTv2.R")
source("~/Personal/ChAMP-Dev/ChAMP/R/QC.GUI.R")

QC.GUI(beta=myLoad$beta, pheno = myLoad$pd$Sample_Group, arraytype = "EPICv2")
```

Also, champ.QC() function works, but I mostly only use QC.GUI() above.

``` r
source("~/Personal/ChAMP-Dev/ChAMP/R/champ.QC.R")

champ.QC(beta = myLoad$beta, pheno = myLoad$pd$Sample_Group)
```

## 5. Normalisation

Since minfi is not working now, currently, champ.norm() only support
`BMIQ` and `PBC` method. `SWAN` and `FunctionalNormalisation` is not
working for now.

``` r
library("doParallel")

source("~/Personal/ChAMP-Dev/ChAMP/R/champ.norm.R")
source("~/Personal/ChAMP-Dev/ChAMP/R/champ.BMIQ.R")


# I recommand BMIQ, which ONLY works for beta matrix.
myNorm <- champ.norm(myLoad$beta, method = "BMIQ", arraytype = "EPICv2")

# or PBC
source("~/Personal/ChAMP-Dev/ChAMP/R/DoPBC.R")
myNorm <- champ.norm(myLoad$beta, method = "PBC", arraytype = "EPICv2")
```

Below is the normalised data from QC.GUI()

``` r
QC.GUI(beta=myNorm, pheno = myLoad$pd$Sample_Group, arraytype = "EPICv2")
```

## 6. SVD Check

Than champ.SVD() should be used to check confounding effect.

``` r
source("~/Personal/ChAMP-Dev/ChAMP/R/champ.SVD.R")

champ.SVD(myNorm, pd = myLoad$pd)
```

The result shows there are confounding effect for Array, so I want to
adjust it with Combat.

## 7. Combat Batch Effect Adjustment

Below code works for batch correction, but I do admit that current
champ.runCombat() need to improve and failed constantly. A better method
(maybe scBatch?) need to be incoprated.

``` r
library("sva")

source("~/Personal/ChAMP-Dev/ChAMP/R/champ.runCombat.R")

myCombat <- champ.runCombat(beta=myNorm,pd=myLoad$pd,batchname=c("Array"))
```

## 7. DMP and GUI

Finally, below is the EWAS analysis, firstly do DMP analysis.

``` r
library("limma")
source("~/Personal/ChAMP-Dev/ChAMP/R/champ.DMP.R")

# Below code exist because the origin sample group contains space, replace it with `_` works.
tmp_pheno <- myLoad$pd$Sample_Group
tmp_pheno <- gsub(" ", "_", tmp_pheno)

myDMP <- champ.DMP(beta = myNorm,pheno=tmp_pheno, arraytype = "EPICv2")
```

Then, CpG.GUI() can be used to check the DMP result:

``` r
library("shiny")
library("shinythemes")
library("plotly")

source("~/Personal/ChAMP-Dev/ChAMP/R/DMP.GUI.R")

DMP.GUI(DMP=myDMP[[1]], beta=myNorm, pheno=myLoad$pd$Sample_Group)
```

It works nicely.

## 7. DMR and GUI

This is a challenge part. Firstly I believe ProbeLasso and DMRcate is
not working, the formal requires a special king of old annotation, which
I have not generated (but I don’t want to use the old format anymore, so
I will change the code to fit new annotation). DMRcate package can’t be
install on my computer…Like some user complained, I don’t know why yet.

So, the Bumphunter is the only choice here. Actually this is the one I
recommand.

``` r
library("bumphunter")

source("~/Personal/ChAMP-Dev/ChAMP/R/champ.DMR.R")
tmp_pheno <- myLoad$pd$Sample_Group
tmp_pheno <- gsub(" ", "_", tmp_pheno)

myDMR <- champ.DMR(beta=myNorm, pheno=tmp_pheno, method="Bumphunter", arraytype = "EPICv2")
```

Below is the DMP.GUI

``` r
library("shiny")
library("shinythemes")
library("plotly")

source("~/Personal/ChAMP-Dev/ChAMP/R/DMR.GUI.R")

DMR.GUI(DMR=myDMR, beta=myNorm, pheno=tmp_pheno, arraytype = "EPICv2")
```

It works well, just a bit slow, I should accelerate it a bit with
data.table….

## 7. Block and GUI

Similar as DMR, Block should works:

``` r
source("~/Personal/ChAMP-Dev/ChAMP/R/champ.Block.R")
tmp_pheno <- myLoad$pd$Sample_Group
tmp_pheno <- gsub(" ", "_", tmp_pheno)

myBlock <- champ.Block(beta=myNorm, pheno = tmp_pheno, arraytype = "EPICv2")
```

Block.GUI is not working for now.
