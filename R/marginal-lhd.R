# Functions for preparing and parsing MCMCTree files for
# marginal likelihood calculation

# TODO: Error cheking, check for ctlf and betaf file existence
#' Prepare mcmctree or bpp control files for marginal likelihood calculation
#'
#' @param beta numeric vector of beta values
#' @param ctlf character, mcmctree or bpp control file template
#' @param betaf character, file onto which to write selected beta values
#'
#' @details
#' This function generates a set of \code{n} directories each
#' containing a modified \code{ctlf} control file with the appropriate beta
#' value to run mcmctree (or bpp) to obtain MCMC samples under the required
#' power-posterior distribution. For the general theory of marginal likelihood
#' calculation with power posteriors see Yang (2014).
#'
#' The beta values are printed to \code{betaf}.
#'
#' @references
#' Yang Z (2014) \emph{Molecular Evolution: A Statistical Approach}. Oxford
#' University Press. Pages 256--260.
#'
#' @seealso
#' \code{\link{make.beta}}, \code{\link{stepping.stones}}
#'
#' @author Mario dos Reis
#'
#' @export
make.bfctlf <- function(beta, ctlf="mcmctree.ctl", betaf="beta.txt") {
  b <- beta
  b[b==0] <- 1e-300
  n <- length(b)

  bf <- "beta"
  for (i in 1:n) {
    dir.create(as.character(i))
    cat(paste("BayesFactorBeta = ", b[i], "\n", sep=""), file=bf)
    newf <- paste(i, "/", ctlf, sep="")
    file.append(file1=newf, file2=c(ctlf, bf))
  }
  unlink(bf)
  write(b, betaf, ncol=1)
}

#' Make beta values for marginal likelihood calculation
#'
#' @description
#' Make appropriate beta values
#'
#' @param n numeric, number of beta points
#' @param method character, the method to choose the beta points, see details
#' @param a numeric, exponent for stepping stones beta generation, see details
#'
#' @details
#' If \code{method = "step-stones"}, the beta values are given by the
#' formula
#'
#' \deqn{\beta_{i}=\left(\frac{i-1}{n}\right)^{a}.}
#'
#' Values of \code{a} between 5 to 8 appear appropriate. Large \code{a} values
#' produce beta values close to zero.
#'
#' If \code{method = "gauss-quad"}, the beta values are calculated according to
#' the \code{n} Gauss-Legendre quadrature rule (see Rannala and Yang, 2017).
#'
#' @return
#' Numeric vector with \code{n} beta values
#'
#' @seealso
#' The generated beta values are suitable input for \code{\link{make.bfctlf}}.
#'
#' @author Mario dos Reis
#'
#' @references
#' Rannala B and Yang Z (2017) Efficient Bayesian species tree inference under
#' the multispecies coalescent. \emph{Systematic Biology}, 66: 823--842.
#'
#' @export
make.beta <- function(n, method=c("step-stones", "gauss-quad"),  a=5) {
  method <- match.arg(method)
  if (method == "step-stones")
    b <- .stepping.stones.beta(n, a)
  if (method == "gauss-quad")
    b <- .gauss.quad.beta(n)
  return (b)
}

#' Estimate marginal likelihood by stepping stones
#'
#' @description
#' Estimate the marginal likelihood using the stepping stones
#' method from a sample of \code{n} power posterior MCMC chains sampled
#' with mcmctree (or bpp).
#'
#' @param mcmcf character, mcmc output file name
#' @param betaf character, file with beta values
#' @param se logical, whether to calculate the standard error
#'
#' @details
#' The MCMC samples should be stored in a directory structure created
#' by \code{make.bfctlf} with \code{method = "step-stones"}. The function will
#' read the stored log-likelihood values and calculate the log-marginal
#' likelihood.
#'
#' If \code{se = TRUE}, an approximation based on the Delta method is used to
#' calculate the standard error (see Xie et al. 2011). Warnings are given if the
#' approximation appears unreliable.
#'
#' @return
#' A list with components \code{logml}, the log-marginal likelihood estimate;
#' \code{se}, the standard error of the estimate; \code{mean.logl}, the mean of
#' log-likelihood values sampled for each beta; and \code{b}, the beta values
#' used.
#'
#' @references
#' Xie et al. (2011) Improving marginal likelihood estimation for Bayesian
#' phylogenetic model selection. \emph{Systematic Biology}, 60: 150--160.
#'
#' @seealso
#' \code{\link{make.bfctlf}} to prepare directories and mcmctree or bpp control
#' files to calculate the power posterior.
#'
#' @author Mario dos Reis
#'
#' @export
stepping.stones <- function(mcmcf="mcmc.txt", betaf="beta.txt", se=TRUE) {
  b <- c(scan(betaf), 1)
  n <- length(b) - 1
  lnLs <- list()
  Ls <- list()
  C <- zr <- ess <- vzr <- mlnl <- numeric(n)
  bdiff <- diff(b)

  for (i in 1:n) {
    lnLs[[i]] <- na.omit(read.table(paste(i, "/", mcmcf, sep=""), header=TRUE, fill=TRUE)$lnL)
    mlnl[i] <- mean(lnLs[[i]])
    lnLs[[i]] <- bdiff[i] * lnLs[[i]]
    C[i] <- max(lnLs[[i]])
    Ls[[i]] <- exp(lnLs[[i]] - C[i])
    zr[i] <- mean(Ls[[i]])
  }

  if (se) {
    for (i in 1:n) {
      ess[i] <- coda::effectiveSize(Ls[[i]])
      vzr[i] <- var(Ls[[i]]) / ess[i]
      # the delta approximation does not work well if vzr/zr^2 > 0.1
      if (vzr[i] / zr[i]^2 > 0.1)
        warning ("unreliable se: var(r_k)/r_k^2 = ", vzr[i] / zr[i]^2, " > 0.1 for b = ", b[i])
    }
    vmlnl <- sum(vzr / zr^2)
  } else {
    vmlnl <- NA
  }

  lnml <- sum(log(zr) + C)
  return ( list(logml=lnml, se=sqrt(vmlnl), mean.logl=mlnl, b=b[1:n]) )
}

#' Estimate marginal likelihood by thermodynamic integration
#'
#' @description
#' Estimate marginal likelihood by thermodynamic integration and Gauss-Legendre
#' quadrature from a sample of \code{n} power posterior MCMC chains sampled
#' with mcmctree (or bpp).
#'
#' @param mcmcf character, mcmc output file name
#' @param betaf character, file with beta values
#' @param se logical, whether to calculate the standard error
#'
#' @details
#' The MCMC samples should be stored in a directory structure created
#' by \code{make.bfctlf} with \code{method = "gauss-quad"}. The function
#' will read the stored log-likelihood values and calculate the log-marginal
#' likelihood.
#'
#' Numerical integration is done using Gauss-Legendre quadrature. See Rannala
#' and Yang (2017) for details (also dos Reis et al. 2017, Appendix 2).
#'
#' @return
#' A list with components \code{logml}, the log-marginal likelihood estimate;
#' \code{se}, the standard error of the estimate; \code{mean.logl}, the mean of
#' log-likelihood values sampled for each beta; and \code{b}, the beta values
#' used.
#'
#' @references
#' Rannala B and Yang Z. (2017) Efficient Bayesian species tree inference under
#' the multispecies coalescent. \emph{Systematic Biology} 66: 823-842.
#'
#' dos Reis et al. (2017) Using phylogenomic data to explore the effects of
#' relaxed clocks and calibration strategies on divergence time estimation:
#' Primates as a test case. \emph{bioRxiv}
#'
#' @seealso
#' \code{\link{make.bfctlf}} to prepare directories and mcmctree or bpp control
#' files to calculate the power posterior.
#'
#' @author Mario dos Reis
#'
#' @export
gauss.quad <- function(mcmcf="mcmc.txt", betaf="beta.txt", se=TRUE) {
  b <- scan(betaf)
  n <- length(b)
  lnLs <- list()
  w <- glqrules[[n]]$w  # the w's are symmetrical
  ess <- vv <- numeric(n)

  for (i in 1:n) {
    lnLs[[i]] <- na.omit(read.table(paste(i, "/", mcmcf, sep=""), header=TRUE, fill=TRUE)$lnL)
  }

  mlnl <- sapply(lnLs, mean)
  lnml <- sum( mlnl * w / 2 )

  if (se) {
    for (i in 1:n) {
      ess[i] <- coda::effectiveSize(lnLs[[i]])
      vv[i] <- var(lnLs[[i]]) / ess[i]
    }
    vmlnl <- sum(w^2 * vv) / 4
  } else {
    vmlnl <- NA
  }

  return ( list(logml=lnml, se=sqrt(vmlnl), mean.logl=mlnl, b=b) )
}

#' Calculate Bayes factors and posterior model probabilities
#'
#' @param ... list of marginal likelihood objects, see details
#' @param prior numeric, the prior model probabilities
#' @param boot logical, whether to perform parametric boostrap of probabilities
#' @param n numeric, number of bootstrap samples
#' @param prob numeric, the probability used to calculate the boostrap CI
#'
#' @details
#' Input is a list of marginal likelihood objects, with each object generated by
#' either \code{stepping.stones()} or \code{gauss.quad()}. If \code{boot = 
#' TRUE}, parametric bootstrap is performed by assuming the log-marginal 
#' likelihood estimates are normally distributed with standard deviation equal 
#' to the standard error. The re-sampled \code{n} marginal log-likelihoods are
#' used to estimate re-sampled posterior probabilities and to calculate an 
#' equal-tail bootstrap confidence interval for these.
#' 
#' Note that the length of \code{prior} should be the same as the number of
#' models being compared. The \code{prior} is rescaled so that 
#' \code{sum(prior) == 1}.
#' 
#' @examples 
#' # See Table 5 in dos Reis et al. (2018, Syst. Biol., 67: 594-615)
#' # Bayesian selection of relaxed clock models for the 1st and 2nd sites 
#' # of mitochondrial protein-coding genes of primates
#' # Models: strick clock, independent-rates, and autocorrelated-rates
#' sc <- list(); sc$logml <- -16519.03; sc$se <- .01
#' ir <- list(); ir$logml <- -16480.58; ir$se <- .063
#' ar <- list(); ar$logml <- -16477.82; ar$se <- .035
#' bayes.factors(sc, ir, ar)
#' bayes.factors(sc, ir, ar, prior=c(.25,.5,.25))
#' bayes.factors(sc, ir, ar, prior=c(0,1,0))
#'
#' @return
#' A list with elements \code{bf} and \code{logbf}, the Bayes factors and
#' log-Bayes factors; \code{pr}, the posterior model probabilities; \code{prior}
#' the prior model probabilities and, if \code{boot = TRUE}, \code{pr.ci} the
#' equal-tail bootstrap confidence interval.
#'
#' @author Mario dos Reis
#'
#' @export
bayes.factors <- function(..., prior=NULL, boot=TRUE, n=1e4, prob=0.95) {
  model <- list(...)
  N <- length(model)
  logml <- numeric(N)
  
  if (is.null(prior)) prior <- rep(1, N)

  for (i in 1:N) logml[i] <- model[[i]]$logml

  logbf <- logml - max(logml)
  bf <- exp(logbf)
  pr <- bf * prior / sum(bf * prior)
  
  rtn <- list(bf=bf, logbf=logbf, pr=pr, prior=prior / sum(prior))
  
  # calculates CI analytically when N = 2 models
  # does not appear necessary as this converges to the bootstrap
  # estimate when n is very large
  # (needs to be tested for various priors)
  # if (!boot & N == 2) { 
  #   logbfsd <- sqrt(model[[1]]$se^2 + model[[2]]$se^2)
  #   mm <- matrix(0, ncol=2, nrow=2)
  #   prob = 1 - prob
  #   qn <- qnorm(prob/2, 0, 1)
  #   mm[1,] <- c(logbf[1], logbf[2] + qn * logbfsd) + log(prior)
  #   mm[2,] <- c(logbf[1], logbf[2] - qn * logbfsd) + log(prior)
  #   bfm <- exp(mm - apply(mm, 1, max))
  #   prm <- t(bfm / apply(bfm, 1, sum))
  #   prm[1,] <- rev(prm[1,])
  #   colnames(prm) <- c(prob/2, 1 - prob/2)
  #   
  #   rtn$pr.ci <- prm
  # }
  
  # calculate parametric bootstrap CIs for posterior probabilities
  if (boot) {
    mm <- matrix(0, ncol=N, nrow=n)
    for (i in 1:N) {
      mm[,i] <- rnorm(n, mean = logml[i] + log(prior[i]), sd = model[[i]]$se)
    }
    bfm <- exp(mm - apply(mm, 1, max))
    prm <- bfm / apply(bfm, 1, sum)
    prob = 1 - prob
    prci <- apply(prm, 2, quantile, probs=c(prob/2, 1 - prob/2), na.rm=TRUE)
    #colnames(prm) <- c(prob/2, 1 - prob/2)
    
    rtn$pr.ci <- t(prci)
  }

  return ( rtn )
}

# gauss-laguerre beta generator function
.gauss.quad.beta <- function(n) {
  x <- -glqrules[[n]]$x # reverse the x's
  b <- (1 + x) / 2
  return (b)
}

# stepping stones beta generator function
.stepping.stones.beta <- function(n, a) {
  return ( ((1:n-1)/n)^a )
}

#' Generate block bootstrap replicates of sampled power likelihoods
#' 
#' @param R numeric, number of bootstrap replicates
#' @param p numeric, block length, giving as a proportion of the MCMC sample size
#' @param mcmcf character, mcmc output file name
#' @param betaf character, file with beta values
#' @param preff character, prefix for files storing boot replicates
#' 
#' @details 
#' Block bootstrap replicates are generated using the stationary boostrap method
#' of Politis and Romano (1994). The replicates are stored in files named using
#' \code{preff} and the replicate number. For example, if \code{preff = "lnL"}
#' (the default) then the files are lnL0.txt, lnL1.txt, lnL2.txt, ..., etc, with
#' lnL0.txt corresponding to the original log-likelihood sample. Replicates are
#' stored within the directories corresponding to the appropriate beta values.
#' The collection of files can grow large quickly so you may want to use a small
#' number of replicates (say R = 10 to R = 100).
#' 
#' This function uses code extracted from the boot package by Canty and Ripley.
#' 
#' @seealso 
#' \link{stepping.stones.boot} and \link{tsboot} (from the boot package).
#' 
#' @references 
#' Politis and Romano (1994) The stationary boostrap. \emph{J. Am. Stat. Assoc.}, 
#' 89: 1303--1313.
#' 
#' @export
block.boot <- function(R, p=0.1, mcmcf="mcmc.txt", betaf="beta.txt", preff="lnL") {
  nb <- length(scan(betaf))
  lnLs <- list()
  
  # iterate over directories corresponding to the beta points
  for (i in 1:nb) {
    lnLs[[i]] <- na.omit(read.table(paste(i, "/", mcmcf, sep=""), header=TRUE, fill=TRUE)$lnL)
    write.table(data.frame(lnL=lnLs[[i]]), file = paste(i, "/", preff, "0.txt", sep=""))
    n <- n.sim <- length(lnLs[[i]])
    
    # iterate over bootstrap replicates
    for (j in 1:R) {
      ta.out <- ts.array(n=n, n.sim=n.sim, R=1, l=n*p, sim="geom", endcorr=TRUE)
      ends <- cbind(ta.out$starts[1,], ta.out$lengths[1,])
      inds <- apply(ends, 1L, make.ends, n)
      inds <- if (is.list(inds)) matrix(unlist(inds)[1L:n.sim], n.sim, 1L)
      else matrix(inds, n.sim, 1L)
      write.table(data.frame(lnL=(lnLs[[i]])[inds]), file = paste(i, "/", preff, j, ".txt", sep=""))
    }
  }
}

#' Estimate marginal likelihood from bootstrap replicates
#' 
#' @param R numeric, number of bootstrap replicates used
#' @param betaf character, file with beta values
#' @param preff character, prefix for files storing boot replicates
#' 
#' @details 
#' \code{stepping.stones.boot} and \code{gauss.quad.boot} are used to calculate
#' the marginal likelihoods on bootstrap replicates using the stepping stones 
#' and gaussian quadrature methods respectively. The replicates must have
#' been generated using \link{block.boot}.  
#' 
#' @return 
#' A list with components \code{logml}, the original log-marginal likelihood
#' estimate, \code{logmlR}, the vector of log-marginal likelihood estimates on
#' the boostrap replicates, \code{se} and \code{ci}, the standard error and 95\%
#' credibility interval of \code{logml} calculated on the bootstrap replicates,
#' and \code{b}, the beta values used.
#' 
#' @seealso 
#' \link{block.boot}, \link{stepping.stones} and \link{gauss.quad}.
#' 
#' @name logL.boot
NULL

#' @rdname logL.boot
#' @export
stepping.stones.boot <- function(R, betaf="beta.txt", preff="lnL") {
  return( .logL.boot(R, betaf, preff, method="step-stones") )
}

#' @rdname logL.boot
#' @export
gauss.quad.boot <- function(R, betaf="beta.txt", preff="lnL") {
  return( .logL.boot(R, betaf, preff, method="gauss-quad") )
}

.logL.boot <- function(R, betaf, preff, method) {
  
  if (method == "step-stones") lnL.fun <- stepping.stones
  if (method == "gauss-quad") lnL.fun <- gauss.quad
  
  b <- scan(betaf)
  nb <- length(b)
  lnLR <- numeric(R)
  
  lnLf <- paste(paste(preff, 0, ".txt", sep=""))
  lnL0 <- stepping.stones(lnLf, betaf, FALSE)$logml
  
  for (j in 1:R) {
    lnLf <- paste(paste(preff, j, ".txt", sep=""))
    lnLR[j] <- lnL.fun(lnLf, betaf, FALSE)$logml
  }
  return(list(logml=lnL0, logmlR=lnLR, se=sd(lnLR), ci=quantile(lnLR, c(.025, .975)), b=b))
}

# Block bootstrap functions: code lifted from boot package,
# file bootfuns.q by Canti and Ripley in
# https://github.com/cran/boot
# copyright (C) 1997-2001 Angelo J. Canty
# corrections (C) 1997-2011 B. D. Ripley
# Note, endcorr should be set to TRUE when using sim = "geom", and n.sim = n
ts.array <- function(n, n.sim, R, l, sim, endcorr)
{
  #
  #  This function finds the starting positions and lengths for the
  #  block bootstrap.
  #
  #  n is the number of data points in the original time series
  #  n.sim is the number require in the simulated time series
  #  R is the number of simulated series required
  #  l is the block length
  #  sim is the simulation type "fixed" or "geom".  For "fixed" l is taken
  #	to be the fixed block length, for "geom" l is the average block
  #	length, the actual lengths having a geometric distribution.
  #  endcorr is a logical specifying whether end-correction is required.
  #
  #  It returns a list of two components
  #  starts is a matrix of starts, it has R rows
  #  lens is a vector of lengths if sim="fixed" or a matrix of lengths
  #	corresponding to the starting points in starts if sim="geom"
  endpt <- if (endcorr) n else n-l+1
  cont <- TRUE
  if (sim == "geom") {
    len.tot <- rep(0,R)
    lens <- NULL
    while (cont) {
      #            inds <- (1L:R)[len.tot < n.sim]
      temp <- 1+rgeom(R, 1/l)
      temp <- pmin(temp, n.sim - len.tot)
      lens <- cbind(lens, temp)
      len.tot <- len.tot + temp
      cont <- any(len.tot < n.sim)
    }
    dimnames(lens) <- NULL
    nn <- ncol(lens)
    st <- matrix(sample.int(endpt, nn*R, replace = TRUE), R)
  } else {
    nn <- ceiling(n.sim/l)
    lens <- c(rep(l,nn-1), 1+(n.sim-1)%%l)
    st <- matrix(sample.int(endpt, nn*R, replace = TRUE), R)
  }
  list(starts = st, lengths = lens)
}

make.ends <- function(a, n)
{
  #  Function which takes a matrix of starts and lengths and returns the
  #  indices for a time series simulation. (Viewing the series as circular.)
  mod <- function(i, n) 1 + (i - 1) %% n
  if (a[2L] == 0) numeric()
  else  mod(seq.int(a[1L], a[1L] + a[2L] - 1, length.out = a[2L]), n)
}