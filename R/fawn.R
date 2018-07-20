## ----data----------------------------------------------------------------
load("../data/fawnData.gzip")
ls()

## ----gibbs-fn------------------------------------------------------------
source("fawn-mcmc-dlife.R")
args(scrOpenD)

## ----fm1,eval=TRUE,cache=TRUE,results='hide',size='small'----------------
fm1 <- scrOpenD(y.data=yDec,             ## fawn x trap x occasion encounter histories
                bday.range=bdRangeDec,   ## birth date ranges
                x=x,                     ## trap locations
                M=150,                   ## data augmentation size
                niter=100,               ## MCMC iterations. Should be >5000 in practice
                buffer=700,              ## buffer defining the state-space
                lifemod="exponential",   ## either 'exponential' or 'weibull'
                Rage=180,                ## age at recruitment)
                birthWindow=c(1, 150),   ## timeframe when births occur. 1=December first.
                upper.life=365,          ## lifetimes are right-censored after 365 days
                oper=operDec,            ## binary operational status matrix
                report=10,               ## display information every 10 iterations
                ## tune order: sig0, sig1, lam0, bday.bar, bay.var,
                ##             bday, omega0, omega1, life, s
                tune=c(25, 3.6, 0.023, 5, 4.4, 30,
                       0.18, 0.25, 55, 300))

## ----str1----------------------------------------------------------------
str(fm1)

## ----coda----------------------------------------------------------------
library(coda) ## Install the 'coda' package if you don't already have it
mc1 <- mcmc(fm1$samples)
#plot(mc1[,1:4]) # First 4 parameters only
plot(mc1)

## ----fm1p,eval=FALSE-----------------------------------------------------
## library(parallel)
## (nCores <- min(detectCores()-1,8))
## cl1 <- makeCluster(nCores)
## clusterExport(cl1, c("scrOpenD", "x",
##                      "yDec", "operDec", "bdRangeDec"))
## clusterSetRNGStream(cl1, iseed=98210)
## 
## out1p1e <- clusterEvalQ(cl1, {
##     fm1p <- scrOpenD(y.data=yDec,
##                       bday.range=bdRangeDec,
##                       x=x,
##                       M=150,
##                       niter=7000, ## 72000 iterations were used for the paper
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

