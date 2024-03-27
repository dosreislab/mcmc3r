#' Birth-death process with species sampling
#' 
#' @param x numeric, x values for which the density will be plotted. This 
#' `dbd` function is passed to R function `curve` via variable `expr`. 
#' The values you type as arguments for variables `from` and `to` in
#' function `curve` will be used by variable `x` in this function 
#' (see examples for details).
#' @param lambda numeric, birth rate.
#' @param mu numeric, death rate.
#' @param rho numeric, proportion of species sampled.
#' @param t1 numeric, age of the phylogeny's root.
#' 
#' @description
#' Kernel density function for the birth-death process with species sampling.
#' 
#' @details 
#' MCMCtree uses the BD kernel to generate the prior on node ages for those
#' nodes without fossil calibrations. You can look at the examples below for
#' some suggestions. Note rho must be between 0 and 1. The special case mu =
#' lambda, rho=0 gives a uniform density. See Yang and Rannala (2006) for full
#' details.
#' 
#' @return 
#' A numeric vector of probability densities.
#' 
#' @references 
#' Yang and Rannala. (2006) Bayesian Estimation of Species Divergence Times
#' Under a Molecular Clock Using Multiple Fossil Calibrations with Soft Bounds.
#' \emph{Mol. Biol. Evol.}, 23: 212–226.
#' 
#' Yang (2014) Molecular Evolution: A Statistical Approach. \emph{Oxford 
#' University Press} 
#' 
#' @examples
#' # Reproduce Fig. 10.10 from Yang (2014)
#' # (a) lambda = mu = 1, rho = 0 (uniform density):
#' curve(dBD(x, 1, 1, 0), xlim=c(0, 1), ylim=c(0, 4), xaxs="i", yaxs="i")
#' 
#' # (b) lambda = 10, mu = 5, rho = 0.01 (old node ages, useful for diversified 
#' # sampling):
#' curve(dBD(x, 10, 5, .01), from=0, to=1, lty=2, add=TRUE)
#' 
#' # (c) lambda = 10, mu = 5, rho = 0.001 (old node ages, useful for diversified 
#' # sampling):
#' curve(dBD(x, 10, 5, .001), from=0, to=1, lty=3, add=TRUE)
#' 
#' # (d) lambda = 10, mu = 5, rho = 0.99 (young node ages, useful for dense
#' # sampling of diverse phylogenies):
#' curve(dBD(x, 10, 5, .99), from=0, to=1, lty=4, add=TRUE)
#' 
#' @author
#' Mario dos Reis
#' 
#' @export
dBD <- function(x, lambda, mu, rho, t1 = 1) {
  d <- numeric(length(x))
  d[x < 0] <- 0
  d[x > t1] <- 0
  ii <- which(x >= 0 & x <= t1)
  
  
  if (mu == lambda) {
    d[ii] <-  (1 + rho * lambda) / (t1 * (1 + rho * lambda * x[ii])^2)
  }
  else {
    a <- mu - lambda
    p.0t <- function(x) {
      -rho * a / (rho * lambda + (lambda * (1 - rho) - mu) * exp (a * x))
    }
    p1.t <- p.0t(x[ii])^2 * exp(a * x[ii]) / rho
    v.t1 <- 1 - p.0t(t1) * exp(a * t1) / rho
    d[ii] <- lambda * p1.t / v.t1
  }
  
  return ( d )
}