## ----data----------------------------------------------------------------
load("../data/fawnData.gzip")
source("fawn-mcmc-dlife.R")

## ----fm1,eval=FALSE------------------------------------------------------
## fm1 <- scrOpenD(y.data=yDec,
##                 bday.range=bdRangeDec,
##                 x=jd9$x,
##                 M=150,
##                 niter=100,
##                 buffer=700,
##                 lifemod="exponential",
##                 Rage=180,
##                 birthWindow=c(1, 150),
##                 upper.life=365,
##                 oper=operDec,
##                 report=10,
##                 ## tune order: sig0, sig1, lam0, bday.bar, bay.var,
##                 ##             bday, omega0, omega1, life, s
##                 tune=c(25, 3.6, 0.023, 5, 4.4, 30,
##                        0.18, 0.25, 55, 300))

## ----str1----------------------------------------------------------------
str(fm1)

## ----coda----------------------------------------------------------------
library(coda)
mc1 <- mcmc(fm1$samples)
#plot(mc1[,1:4]) # First 4 parameters only
plot(mc1) # First 4 parameters only

## ----fm1p,eval=FALSE-----------------------------------------------------
## library(parallel)
## (nCores <- min(detectCores()-1,8))
## cl1 <- makeCluster(nCores)
## clusterExport(cl1, c("scrOpenD", "jd9", "y9",
##                      "yDec", "operDec", "bdRangeDec"))
## clusterSetRNGStream(cl1, iseed=98210)
## 
## out1p1e <- clusterEvalQ(cl1, {
##     fm1p <- scrOpenD(y.data=yDec,
##                       bday.range=bdRangeDec,
##                       x=jd9$x,
##                       M=150,
##                       niter=72000,
##                       buffer=700,
##                       lifemod="exponential",
##                       Rage=180,
##                       birthWindow=c(1, 150),
##                       upper.life=365,
##                       oper=operDec,
##                       report=0,
##                       tune=c(25, 3.6, 0.023, 5, 4.4, 30,
##                              0.18, 0.25, 55, 300))
##     return(fm1p)
## })

## ----mc1p,eval=FALSE-----------------------------------------------------
## mc1p <- as.mcmc.list(lapply(out1p1e, function(x) mcmc(x$samples)))

