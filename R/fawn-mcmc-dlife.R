scrOpenD <- function(
    y.data,     ## nDetected x nTraps x nOccasions encounter history array
    bday.range, ## nDetected x 2 matrix of birth date ranges
    x,          ## nTraps x 2 matrix of trap coords
    M,          ## data augmentation size parameter
    niters,     ## MCMC iterations
    buffer,     ## buffer used to define the state-space
    lifemod=c("weibull", "exponential"),
    Rage=180,    ## Recruitment age
    birthWindow, ## Range defining possible bdays
    upper.life,  ## upper limit of lifetimes to be estimated.
                 ##   lifetimes>upper.life are right-censored.
    oper=NULL,   ## JxK matrix indicating if trap was operational
    report=10,   ## Report status every x iterations
    ## tune order: sigma, lam0,
    tune)
{

    lifemod <- match.arg(lifemod)
    n0 <- nrow(y.data)  # nDetected
    J <- ncol(y.data)   # nTraps
    K <- dim(y.data)[3] # nOccasions

    ## Initial values
    s <- cbind(runif(M, min(x[,1]), max(x[,1])),
               runif(M, min(x[,2]), max(x[,2])))
    sigma0 <- runif(1, 300, 400)
    sigma1 <- runif(1, 4, 10)
    lam0 <- runif(1, 0.4, 0.5)
    psi <- runif(1, 0.9, 0.99)
    bday.bar <- rnorm(1, 70, 5)
    bday.var <- runif(1, 8, 12) ## SD, not variance
    omega0 <- runif(1, 0.9, 1.1)
    omega1 <- runif(1, 100, 200)
    bday <- rnorm(M, bday.bar, bday.var)
    bdayF <- floor(bday+0.5)
    life <- ceiling(rweibull(M, omega0, omega1))
    dday <- bdayF+life-1

    first <- rep(Inf, M)
    last <- rep(-Inf, M)
    lowerBday <- birthWindow[1]
    upperBday <- birthWindow[2]
    bday.range.in <- bday.range
    if(any(bday.range.in[,1]<lowerBday)) {
        stop("bday.range.in[,1] values must be greater than or equal to birthWindow[1]")
    }
    if(any(bday.range.in[,2]>upperBday)) {
        stop("bday.range.in[,2] values must be less than or equal to birthWindow[2]")
    }
    bday.range <- cbind(rep(lowerBday, M), rep(upperBday, M))
    bday.range[1:n0,] <- bday.range.in
    for(i in 1:n0) {
        y.i <- y.data[i,,]
        traps.i <- which(rowSums(y.i)>0)
        s[i,] <- colMeans(x[traps.i,,drop=FALSE])
        dets.i <- which(colMeans(y.i)>0)
        first[i] <- min(dets.i)
        last[i] <- max(dets.i)
        bday[i] <- bday.range[i,1]
        if(bday[i]>first[i])
            stop("Bad bday")
        bdayF[i] <- floor(bday[i]+0.5)
        life[i] <- (last[i]+10)-bdayF[i]+1
    }
    dist2mat <- matrix(NA, M, J)
    for(j in 1:J) {
        dist2mat[,j] <- (s[,1]-x[j,1])^2 + (s[,2]-x[j,2])^2
    }
    dist2 <- array(dist2mat, c(M,J,K))

    ll.bday <- dnorm(bday, bday.bar, bday.var, log=TRUE)

    ## Lifetime modeled as outcome of categorical distribution
    ## with probabilities defined by parametric model
    if(missing(upper.life))
        upper.life <- max(365,K+abs(bday.range.in[,1]))+1
    life[life>upper.life] <- upper.life+1
    dday <- bdayF+life-1
    life.poss <- 1:upper.life
    nlife.poss <- length(life.poss)+1
    life.breaks <- seq(from=0, to=upper.life, by=1)
    up1 <- upper.life+1L
    if(lifemod=="weibull") {
        cdf.life <- pweibull(life.breaks, omega0, omega1)
    } else {
        if(lifemod=="exponential") {
            omega0 <- 1
            cdf.life <- pexp(life.breaks, 1/omega1)

        }
    }
    pi.life <- cdf.life[-1] - cdf.life[-up1]
    pi.life[upper.life+1] <- 1-sum(pi.life)
    ll.life <- log(pi.life[life])
    ll.life.sum <- sum(ll.life)

    ## Age on bdayF should be 1, not 0
    agemat <- matrix(1, M, K)
    for(i in 1:M) {
        if(bdayF[i]>K)
            next
        if(bdayF[i]>0) {
            age.index <- bdayF[i]:K
            agemat[i,age.index] <- 1:length(age.index)
        } else {
            age.at.k1 <- 2-bdayF[i]   ## Age on bdayF should be 1, not 0
            agemat[i,] <- seq(from=age.at.k1, by=1, length.out=K)
        }
    }
    age <- array(agemat[,rep(1:K, each=J)], c(M,J,K))

    sigma <- sigma0*exp(-sigma1/age)


    ## Make lam an array to be conformable with y, and to avoid loops
    if(is.null(oper)) {
        oper <- array(1, c(M, J, K))
    } else {
        oper <- array(oper, c(J, K, M)) ## oper should be JxK matix
        oper <- aperm(oper, c(3,1,2))   ## now it's a MxJxK array
    }
    lam <- lam0*exp(-dist2/(2*sigma^2))*oper
    y <- array(0L, c(M, J, K))
    y[1:nrow(y.data),,] <- y.data
    seen <- rowSums(y) > 0 #(y>0, 1, any)

    b <- ifelse(rowSums(y)>0, 1L, 0L)
    B <- sum(b)
    ## Make z an array to conform with y and avoid loops
    z <- array(0L, c(M, J, K))
    for(i in 1:M) {
        if(b[i]<1)
            next
        zitmp <- z[i,,]
        zitmp[,max(bdayF[i],1):min(dday[i],K)] <- 1L
        z[i,,] <- zitmp
    }

    R <- b %*% life>Rage ## nRecruits

    ll.b <- dbinom(b, 1, psi, log=TRUE)
    ll.y <- dpois(y, lam*z, log=TRUE)
    ll.y.cand <- ll.y
    ll.y.sum <- ll.y.cand.sum <- sum(ll.y)

    ## Proposal distribution for lifetime random variable
    ## Results in much better mixing than a normal or
    ##    log-normal proposal distribution
    ## poss should be a sequence of consecutive integers beginning at 1
    dprop.life <- function(poss, curr, tune, redist=0.1) {
        nposs <- length(poss)
        prop.life <- c(pnorm(poss[-1], curr, tune) -
                       pnorm(poss[-nposs], curr, tune), 0)
        prop.life[1] <- pnorm(poss[2], curr, tune)
        prop.life[nposs] <- 1-sum(prop.life[-nposs])
        ## Redistribute redist/2 proportion of the mass to all cells
        prop.life <- prop.life*(1-redist)
        prop.life <- prop.life+(redist/2)/nposs
        ## The other redist/2 propotion goes to the censor cell
        prop.life[nposs+1] <- 1-sum(prop.life)
        prop.life[curr] <- 0 # No point in proposing current value
        prop.life <- prop.life / sum(prop.life)
        return(prop.life)
    }

    out <- matrix(NA, nrow=niters, ncol=11)
    colnames(out) <- c("sigma0", "sigma1", "lam0", "bday.bar", "bday.var",
                       "omega0", "omega1", "psi", "B", "R", "deviance")

    bs <- bdays <- lifes <- matrix(NA, M, niters)
    ss <- array(NA, c(M, 2, niters))

    cat("\nstarting values =",
        round(c(sigma0, sigma1, lam0, bday.bar, bday.var,
                omega0, omega1, psi, B, R, -2*ll.y.sum), 3), "\n\n")

    reportit <- report>0

    for(iter in 1:niters) {

        if(reportit) {
            if(iter %% report == 0) {
                cat("\niter", iter, format(Sys.time(), "%H:%M:%S"), "\n")
                cat("   current =", round(out[iter-1,], 3), "\n")
##                plot(s)
            }
        }


        ## sample sigma0
        sigma0.cand <- rnorm(1, sigma0, tune[1])
        if(sigma0.cand > 0 & sigma0.cand < 1500) {
            prior <- prior.cand <- 0
            sigma.cand <- sigma0.cand*exp(-sigma1/age)
            lam.cand <- lam0*exp(-dist2/(2*sigma.cand^2))*oper
            ll.y.cand <- dpois(y, lam.cand*z, log=TRUE)
            ll.y.cand.sum <- sum(ll.y.cand)
            if(runif(1) < exp((ll.y.cand.sum + prior.cand) -
                              (ll.y.sum + prior))) {
                ll.y <- ll.y.cand
                ll.y.sum <- ll.y.cand.sum
                lam <- lam.cand
                sigma0 <- sigma0.cand
                sigma <- sigma.cand
            }
        }

        ## sample sigma1
        sigma1.cand <- rnorm(1, sigma1, tune[2])
        if(sigma1.cand > 0) {
            prior <- prior.cand <- 0
            sigma.cand <- sigma0*exp(-sigma1.cand/age)
            lam.cand <- lam0*exp(-dist2/(2*sigma.cand^2))*oper
            ll.y.cand <- dpois(y, lam.cand*z, log=TRUE)
            ll.y.cand.sum <- sum(ll.y.cand)
            if(runif(1) < exp((ll.y.cand.sum + prior.cand) -
                              (ll.y.sum + prior))) {
                ll.y <- ll.y.cand
                ll.y.sum <- ll.y.cand.sum
                lam <- lam.cand
                sigma1 <- sigma1.cand
                sigma <- sigma.cand
            }
        }

        ## sample lam0
        lam0.cand <- rnorm(1, lam0, tune[3])
        lok <- TRUE
        if((lam0.cand >= 0) & lok) {
            lam.cand <- lam0.cand*exp(-dist2/(2*sigma^2))*oper
            ll.y.cand <- dpois(y, lam.cand*z, log=TRUE)
            ll.y.cand.sum <- sum(ll.y.cand)
            prior <- prior.cand <- 0
            if(runif(1) < exp((ll.y.cand.sum + prior.cand) -
                              (ll.y.sum + prior))) {
                ll.y <- ll.y.cand
                ll.y.sum <- ll.y.cand.sum
                lam0 <- lam0.cand
                lam <- lam.cand
            }
        }

        ## sample b (and z)
        bUps <- 0
        for(i in 1:M) {
            if(seen[i])
                next
            bcand <- 1-b[i]
            zi.cand <- matrix(0L, J, K)
            if(dday[i]>0 & bdayF[i]<=K & bcand>0)
                zi.cand[,max(bdayF[i],1):min(dday[i],K)] <- 1L
            if(bcand>0) {
                ll.y.cand[i,,] <- dpois(y[i,,], lam[i,,]*zi.cand, log=TRUE)
            } else {
                ll.y.cand[i,,] <- 0
            }
            ll.b[i] <- dbinom(b[i], 1, psi, log=TRUE)
            ll.b.cand <- dbinom(bcand, 1, psi, log=TRUE)
            if(runif(1) < exp((sum(ll.y.cand[i,,])+ll.b.cand) -
                              (sum(ll.y[i,,])+ll.b[i]) )) {
                b[i] <- bcand
                z[i,,] <- zi.cand
                ll.y[i,,] <- ll.y.cand[i,,]
                ll.b[i] <- ll.b.cand
                bUps <- bUps+1
            }
        }

        B <- sum(b) # nBorn

        ## sample psi
        psi <- rbeta(1, 1+B, 1+M-B)

        ## Sample bday.bar with prior = N(0,100)
        bday.bar.cand <- rnorm(1, bday.bar, tune[4])
        lprior.bday.bar <- lprior.bday.bar.cand <- 0
        if(bday.bar.cand > 0 & bday.bar.cand < 150) { ## Unif(Dec 1, May 1) prior
            ll.bday <- dnorm(bday, bday.bar, bday.var, log=TRUE)
            ll.bday.cand <- dnorm(bday, bday.bar.cand, bday.var, log=TRUE)
            if(runif(1) < exp((sum(ll.bday.cand)+lprior.bday.bar.cand) -
                               (sum(ll.bday)+lprior.bday.bar))) {
                bday.bar <- bday.bar.cand
                ll.bday <- ll.bday.cand
            }
        }

        ## Sample bday.var with prior = G(0.001, 0.001)
        bday.var.cand <- rnorm(1, bday.var, tune[5])
        if(bday.var.cand>0 & bday.var.cand<100) { ## Unif(0,100) prior
            lprior.bday.var <- 0 #dgamma(bday.var, 0.001, 0.001, log=TRUE)
            lprior.bday.var.cand <- 0 #dgamma(bday.var.cand, 0.001, 0.001, log=TRUE)
            ll.bday <- dnorm(bday, bday.bar, bday.var, log=TRUE)
            ll.bday.cand <- dnorm(bday, bday.bar, bday.var.cand, log=TRUE)
            if(runif(1) < exp((sum(ll.bday.cand)+lprior.bday.var.cand) -
                               (sum(ll.bday)+lprior.bday.var))) {
                bday.var <- bday.var.cand
                ll.bday <- ll.bday.cand
            }
        }

        ## Sample birth date
        ## Must jointly update bday and dday
        bdayUps <- 0
        for(i in 1:M) {
            bday.cand <- rnorm(1, bday[i], tune[6])
            bdayF.cand <- floor(bday.cand)
            dday.cand <- floor(bdayF.cand+life[i]-1)
            if((bdayF.cand>first[i]))
                next
            if((bdayF.cand < bday.range[i,1]) || (bdayF.cand > bday.range[i,2]))
               next
            if(bdayF.cand <= K) {
                if(bdayF.cand>0) {
                    age.i <- rep(1, K)      ## Age at bdayF should be 1, not 0
                    age.ind <- bdayF.cand:K
                    age.i[age.ind] <- 1:length(age.ind)
                } else {
                    age.at.k1 <- 2-bdayF.cand   ## Age at bdayF should be 1, not 0
                    age.i <- seq(from=age.at.k1, by=1, length.out=K)
                }
            }
            agemat.i <- matrix(age.i, J, K, byrow=TRUE)
            zi.cand <- matrix(0L, J, K)
            if(dday.cand>0 & bdayF.cand<=K & b[i]>0)
                zi.cand[,max(bdayF.cand,1):min(dday.cand,K)] <- 1L
            sigma.i.cand <- sigma0*exp(-sigma1/agemat.i)
            lam.i.cand <- lam0*exp(-dist2[i,,]/(2*sigma.i.cand^2))*oper[i,,]
            ll.bday[i] <- dnorm(bday[i], bday.bar, bday.var, log=TRUE)
            ll.bday.cand <- dnorm(bday.cand, bday.bar, bday.var, log=TRUE)
            if(b[i]>0) {
                if((bdayF[i]<0) & (bdayF.cand<0)) {
                    ll.y.cand[i,,] <- ll.y[i,,]
                } else {
                    ll.y.cand[i,,] <- dpois(y[i,,], lam.i.cand*zi.cand, log=TRUE)
                }
            } else {
                ll.y.cand[i,,] <- ll.y[i,,] <- 0
            }
            if(runif(1) < exp((sum(ll.y.cand[i,,]) + ll.bday.cand) -
                               (sum(ll.y[i,,]) + ll.bday[i]))) {
                bday[i] <- bday.cand
                bdayF[i] <- bdayF.cand
                dday[i] <- dday.cand
                ll.bday[i] <- ll.bday.cand
                ll.y[i,,] <- ll.y.cand[i,,]
                z[i,,] <- zi.cand
                age[i,,] <- agemat.i
                lam[i,,] <- lam.i.cand
                sigma[i,,] <- sigma.i.cand
                bdayUps <- bdayUps+b[i]
            }
        }

        ## Sample omega0
        if(lifemod=="weibull") {
            omega0.cand <- rlnorm(1, log(omega0), tune[7])
            if(omega0.cand < 5) { # Unif(0,5) prior
                cdf.life.cand <- pweibull(life.breaks, omega0.cand, omega1)
                pi.life.cand <- cdf.life.cand[-1] - cdf.life.cand[-up1]
                pi.life.cand[up1] <- 1-sum(pi.life.cand)
                ll.life.cand <- log(pi.life.cand[life])
                ll.life.sum.cand <- sum(ll.life.cand)
                ll.omega0 <- dlnorm(omega0, log(omega0.cand), tune[7], log=TRUE)
                ll.omega0.cand <- dlnorm(omega0.cand, log(omega0), tune[7], log=TRUE)
                if(runif(1) < exp((ll.life.sum.cand+ll.omega0) -
                                  (ll.life.sum+ll.omega0.cand))) {
                    omega0 <- omega0.cand
                    pi.life <- pi.life.cand
                    ll.life <- ll.life.cand
                    ll.life.sum <- ll.life.sum.cand
                }
            }
        }

        ## Sample omega1
        omega1.cand <- rlnorm(1, log(omega1), tune[8])
        if(omega1.cand < 5000) { # Unif(0,5000) prior
            if(lifemod=="weibull") {
                cdf.life.cand <- pweibull(life.breaks, omega0, omega1.cand)
            } else {
                if(lifemod=="exponential") {
                    cdf.life.cand <- pexp(life.breaks, 1/omega1.cand)
                }
            }
            pi.life.cand <- cdf.life.cand[-1] - cdf.life.cand[-up1]
            pi.life.cand[up1] <- 1-sum(pi.life.cand)
            ll.life.cand <- log(pi.life.cand[life])
            omega1.prior <- 0 #dgamma(omega1, 0.001, 0.001, log=TRUE)
            omega1.cand.prior <- 0#dgamma(omega1.cand, 0.001, 0.001, log=TRUE)
            ll.life.sum.cand <- sum(ll.life.cand)
            ll.omega1 <- dlnorm(omega1, log(omega1.cand), tune[8], log=TRUE)
            ll.omega1.cand <- dlnorm(omega1.cand, log(omega1), tune[8], log=TRUE)
            if(runif(1) < exp((ll.life.sum.cand+omega1.cand.prior+ll.omega1) -
                              (ll.life.sum+omega1.prior+ll.omega1.cand))) {
                omega1 <- omega1.cand
                pi.life <- pi.life.cand
                ll.life <- ll.life.cand
                ll.life.sum <- ll.life.sum.cand
            }
        }

        ## Sample lifetime (and hence death date)
        lifeUps <- 0
        for(i in 1:M) {
            ## Asymmetric proposal distribution
            prop.lcand <- dprop.life(life.poss, life[i], tune[9])
            life.cand <- sample.int(nlife.poss, 1, prob=prop.lcand)
            lprop.life.cand <- log(prop.lcand[life.cand])
            prop.lcurr <- dprop.life(life.poss, life.cand, tune[9])
            lprop.life <- log(prop.lcurr[life[i]])
            if((life.cand < 1) || (life.cand > (up1)))
                next
            dday.cand <- bdayF[i] + life.cand - 1
            if((dday.cand<bdayF[i]) || (dday.cand<last[i]))
                next
            zi.cand <- matrix(0L, J, K)
            if((dday.cand>0) & (bdayF[i]<=K) & (b[i]>0)) {
                zi.cand[,max(bdayF[i],1):min(dday.cand,K)] <- 1L
            }
            ll.life.cand <- log(pi.life[life.cand])
            if(b[i]>0) {
                if((dday[i]>K) & (dday.cand>K)) {
                    ll.y.cand[i,,] <- ll.y[i,,]
                } else {
                    ll.y.cand[i,,] <- dpois(y[i,,], lam[i,,]*zi.cand, log=TRUE)
                }
            } else {
                ll.y.cand[i,,] <- ll.y[i,,] <- 0
            }
            if(runif(1) < exp((sum(ll.y.cand[i,,]) + ll.life.cand + lprop.life) -
                               (sum(ll.y[i,,]) + ll.life[i] + lprop.life.cand))) {
                life[i] <- life.cand
                dday[i] <- dday.cand
                ll.life[i] <- ll.life.cand
                ll.y[i,,] <- ll.y.cand[i,,]
                z[i,,] <- zi.cand
                lifeUps <- lifeUps+b[i]
            }
        }

        R <- b%*%(life>Rage) ## nRecruits
        ll.life.sum <- sum(ll.life)

        ## sample s
        sups <- 0
        for(i in 1:M) {   # note this is "M" in general
            scand1 <- rnorm(1, s[i,1], tune[10])
            scand2 <- rnorm(1, s[i,2], tune[10])
            dtmp2 <- (scand1 - x[,1])^2 +
                     (scand2 - x[,2])^2
            inbuffer <- any(sqrt(dtmp2) < buffer)
            if(!inbuffer)
                next
            lam.cand <- lam0*exp(-dtmp2/(2*sigma[i,,]^2) )*oper[i,,]
            if(b[i]>0) {
                ll.y.cand[i,,] <- dpois(y[i,,], lam.cand*z[i,,], log=TRUE)
            } else {
                ll.y.cand[i,,] <- ll.y[i,,] <- 0
            }
            if(runif(1) < exp(sum(ll.y.cand[i,,]) - sum(ll.y[i,,]))) {
                ll.y[i,,] <- ll.y.cand[i,,]
                s[i,] <- c(scand1, scand2)
                lam[i,,] <- lam.cand
                dist2[i,,] <- dtmp2
                sups <- sups+b[i] # Only update for real guys
            }
        }

        ll.y.sum <- sum(ll.y)
        deviance <- -2*(ll.y.sum)

        if(reportit) {
            if(iter %% report == 0) {
                cat("   Acceptance rates\n")
                cat("     b =", bUps/sum(!seen), "\n")
                cat("     bday =", bdayUps/B, "\n")
                cat("     life =", lifeUps/B, "\n")
                cat("     s =", sups/B, "\n")
                cat("  bday dist = \n")
                print(summary(bday[b==1]))
                cat("  lifetime dist = \n")
                print(summary(life[b==1]))
            }
        }

        out[iter,] <- c(sigma0, sigma1, lam0, bday.bar, bday.var,
                        omega0, omega1, psi, B, R, deviance )
        bs[,iter] <- b
        bdays[,iter] <- bday
        lifes[,iter] <- life
        ss[,,iter] <- s
    }

    last <- list(sigma=sigma, sigma0=sigma0, sigma1=sigma1,
                 lam0=lam0, psi=psi, b=b, z=z, s=s,
                 bday=bday, life=life, dday=dday)
    rng <- .Random.seed
    ret <- list(samples=out, b.post=bs, bday.post=bdays, life.post=lifes, s.post=ss,
                last=last, RNG=rng)
    return(ret)
}
