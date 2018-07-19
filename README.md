# mee-fawn-data
Data and code for Methods in Ecology and Evolution paper (to appear) about modeling fawn survival and recruitment.


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
source("fawn.R") 
```

## Documentation
For an overview of the data and code, you can create a PDF using these `R` commands:

```
install.packages("knitr")
library(knitr) 
knit("fawn.Rnw")
tools::texi2pdf("fawn.tex", clean=TRUE)
system("open fawn.pdf") # If this doesn't open the PDF, you can do it manually
```


