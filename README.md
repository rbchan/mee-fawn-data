# mee-fawn-data
Data and code for [Methods in Ecology and Evolution paper](https://doi.org/10.1111/2041-210X.13068) (*In press*) about modeling fawn survival and recruitment.


## Data
The data are in the `data/fawnData.gzip` file, which can be loaded into `R` using the commands:

```
setwd("CHANGE_ME/mee-fawn-data/data") ## Change the "CHANGE_ME" part of the path
load("fawnData.gzip")
```


## Code
To run the `R` code, you can issue these commands:

```
setwd("CHANGE_ME/mee-fawn-data/R") ## Change the "CHANGE_ME" part of the path
install.packages("coda")
source("fawn.R")                   ## Could take a few minutes
```

## Documentation
For an overview of the data and code, you can create a PDF using these `R` commands:

```
install.packages("knitr")
library(knitr)
setwd("CHANGE_ME/mee-fawn-data/R") ## Change the "CHANGE_ME" part of the path
knit("fawn.Rnw")                   ## Could take a few minutes
tools::texi2pdf("fawn.tex", clean=TRUE)
system("open fawn.pdf") ## If this doesn't open the PDF, open it manually
```


## DOI
[![DOI](https://zenodo.org/badge/141598623.svg)](https://zenodo.org/badge/latestdoi/141598623)


