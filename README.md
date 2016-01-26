## compcodeR

**compcodeR** is an R package that provides extensive functionality for comparing results obtained by different methods for differential expression analysis of (mainly) RNAseq data. It also contains functions for simulating count data and interfaces to several packages for performing the differential expression analysis.

The release version of **compcodeR** can be installed via [Bioconductor](http://www.bioconductor.org/packages/release/bioc/html/compcodeR.html):

```
source("http://bioconductor.org/biocLite.R")
biocLite("compcodeR")
```

The most recent (development) version can also be installed from GitHub directly, using the devtools package:

```
library(devtools)
devtools::install_github("csoneson/compcodeR")
```