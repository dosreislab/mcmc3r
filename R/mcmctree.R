# Functions for preparing and parsing MCMCTree files for
# marginal likelihood calculation

# TODO: Error cheking, check for ctlf and betaf file existence
#' Prepare mcmctree or bpp control files for marginal likelihood calculation
#' 
#' @param beta numeric vector of beta values 
#' @param ctlf character, mcmctree control file template
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
#' University Press. Page xx.
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
#' 
#' @details 
#' The MCMC samples should be stored in a directory structure created
#' by \code{make.bfctlf} with \code{method = "step-stones"}. The function will
#' read the stored log-likelihood values and calculate the log-marginal 
#' likelihood. If \code{se = TRUE}, then the standard error of the estimated is
#' also given (this requires the \code{coda} package).
#' 
#' @return 
#' A list with components \code{logml}, the log-marginal likelihood estimate; 
#' \code{se}, the standard error of the estimate; \code{mean.logl}, the mean of
#' log-likelihood values sampled for each beta; and \code{b}, the beta values
#' used.
#' 
#' @seealso 
#' \code{\link{make.bfctlf}} to prepare directories and mcmctree control files
#' to calculate the power posterior.
#' 
#' @author Mario dos Reis
#' 
#' @export
stepping.stones <- function(mcmcf="mcmc.txt", betaf="beta.txt") {
  b <- c(scan(betaf), 1)
  n <- length(b) - 1
  lnLs <- list()
  Ls <- list()
  C <- zr <- ess <- vzr <- mlnl <- numeric(n)
  bdiff <- diff(b)
  
  for (i in 1:n) {
    lnLs[[i]] <- read.table(paste(i, "/", mcmcf, sep=""), header=TRUE)$lnL
    mlnl[i] <- mean(lnLs[[i]])
    lnLs[[i]] <- bdiff[i] * lnLs[[i]]
    C[i] <- max(lnLs[[i]])
    Ls[[i]] <- exp(lnLs[[i]] - C[i])
    zr[i] <- mean(Ls[[i]])
  }
    

  for (i in 1:n) {
    ess[i] <- coda::effectiveSize(Ls[[i]])
    vzr[i] <- var(Ls[[i]]) / ess[i]
  }
  vmlnl <- sum(vzr / zr^2)
  
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
#' 
#' @details
#' The MCMC samples should be stored in a directory structure created
#' by \code{make.bfctlf} with \code{method = "gauss-quad"}. The function 
#' will read the stored log-likelihood values and calculate the log-marginal 
#' likelihood. If \code{se = TRUE}, then the standard error of the estimated 
#' is also given (this requires the \code{coda} package). 
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
#' \code{\link{make.bfctlf}} to prepare directories and mcmctree control files
#' to calculate the power posterior.
#' 
#' @author Mario dos Reis
#' 
#' @export
gauss.quad <- function(mcmcf="mcmc.txt", betaf="beta.txt") {
  b <- scan(betaf)
  n <- length(b)
  lnLs <- list()
  w <- glqrules[[n]]$w  # the w's are symmetrical
  ess <- vv <- numeric(n)
  
  for (i in 1:n) {
    lnLs[[i]] <- read.table(paste(i, "/", mcmcf, sep=""), header=TRUE)$lnL
  }
  
  mlnl <- sapply(lnLs, mean)
  lnml <- sum( mlnl * w / 2 )
  
  
  for (i in 1:n) {
    ess[i] <- coda::effectiveSize(lnLs[[i]])
    vv[i] <- var(lnLs[[i]]) / ess[i]
  }
  vmlnl <- sum(w^2 * vv) / 4
  
  return ( list(logml=lnml, se=sqrt(vmlnl), mean.logl=mlnl, b=b) )
}

#' Calculate Bayes factors and posterior model probabilities
#' 
#' @param ... list of marginal likelihood objects, see details
#' 
#' @details 
#' Input is a list of marginal likelihood objects, with each object generated
#'  by either \code{stepping.stones()} or \code{gauss.quad()}.
#' 
#' @return 
#' A list with elements \code{bf}, the Bayes factors; and \code{pr}, the 
#' posterior model probabilities.
#' 
#' @author Mario dos Reis
#' 
#' @export
bayes.factors <- function(...) {
  model <- list(...)
  n <- length(model)
  logml <- numeric(n)
  
  for (i in 1:n) logml[i] <- model[[i]]$logml
  
  bf <- exp(logml - max(logml))
  pr <- bf / sum(bf)
  
  return ( list(bf=bf, pr=pr) )
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

#' # Calculate the beta values for Bayes Factor calculation
#' # Points are chosen according to Gauss-Legendre quadrature rules
#' # n: number of beta points
#' #' @export
#' BFbeta <- function(n) {
#'   qr <- glqrules[[n]]
#'   b <- (1 + qr$x) / 2
#'   return(list(x=qr$x, b=b, w=qr$w))
#' }

# make.bfctlf <- function(file="mcmctree.ctl", n) {
#   df <- BFbeta(n)
#   betaf <- "beta"
#   for (i in 1:n) {
#     dir.create(as.character(i))
#     cat(paste("BayesFactorBeta = ", df$b[i], "\n", sep=""), file=betaf)
#     newf <- paste(i, "/", file, sep="")
#     file.append(file1=newf, file2=c(file, betaf))
#   }
#   unlink(betaf)
#   write.table(df, row.names=FALSE, file="beta.txt")
#   # TODO: print w, x, b to a file
#   # TODO: add option for replicate runs
# }

#' # Gauss-Legendre quadrature
#' #' @export
#' glq <- function(logml, b.df) sum(logml * b.df$w) / 2
#' 
#' # Bayes factors
#' # ml: marginal likelihood of models. length(ml) > 1
#' # returns BF and posterior Pr (under uniform prior)
#' #' @export
#' BF <- function(logml) {
#'   if (length(logml) < 2) stop("Provide at least two log-marginal likelihoods.")
#'   logml0 <- max(logml)
#'   lbf <- logml - logml0
#'   bf <- exp(lbf)
#'   sbf <- sum(bf)
#'   return(list(BF=bf, Pr=bf/sbf))
#' }