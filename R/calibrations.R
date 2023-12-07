#' Calibration densities
#' 
#' @description Density, distribution, and quantile functions for calibrations 
#' used in MCMCtree.
#' 
#' @param x  numeric, vector of quantiles
#' @param q, numeric, quantile
#' @param prob, numeric probability
#' @param tL numeric, minimum age
#' @param tU numeric, maximum age
#' @param p  numeric, mode parameter for truncated Cauchy
#' @param c  numeric, tail decay parameter for truncated Cauchy
#' @param pL numeric, minimum probability bound
#' @param pU numeric, maximum probability bound
#' 
#' @details Calculates the density, distribution and quantile functions for the
#'   minimum (dL) calibration, and the density function for the joint (dB) and
#'   maximum (dU) calibration bounds as implemented in MCMCtree. See Yang and
#'   Rannala (2007) and Inoue et al. (2010) for details.
#' 
#' @return A vector of density, probability, or quantile values as appropriate.
#' 
#' @references 
#' Yang and Rannala. (2006) Bayesian Estimation of Species Divergence Times
#' Under a Molecular Clock Using Multiple Fossil Calibrations with Soft Bounds.
#' \emph{Mol. Biol. Evol.}, 23: 212–226.
#'
#' Inoue, Donoghue and Yang (2010) The Impact of the Representation of Fossil
#' Calibrations on Bayesian Estimation of Species Divergence Times. \emph{Syst.
#' Biol.}, 59: 74–89.
#' 
#' @examples 
#' # Plot a minimum bound calibration density:
#' curve(dL(x, 1), from=0, to=10, n=5e2)
#' 
#' # Plot a joint bounds calibration density:
#' curve(dB(x, 1, 6), from=0, to=10, n=5e2)
#' 
#' # Plot a maximum bound calibration density:
#' curve(dU(x, 6), from=0, to=10, n=5e2)
#' 
#' # Probability and quantile function for minimum bound (or truncated-Cauchy):
#' qv <- pL(0:20, tL=1)
#' qL(qv, tL=1)
#' 
#' @author Mario dos Reis
#'
#' @name calibrations
NULL

# L bound calibration:
# density function:
#' @rdname calibrations
#' @export
dL <- function(x, tL, p=0.1, c=1, pL=0.025) {
  a <- 1 - pL
  A <- .5 + atan(p/c) / pi
  theta <- a/pL * 1/(pi*A*c*(1 + (p/c)^2))
  
  dx <- numeric(length(x))
  
  i <- which(x < 0)
  j <- which(0 <= x & x < tL)
  k <- which(x >= tL)
  
  dx[i] <- 0
  dx[j] <- pL * theta/tL * (x[j]/tL)^(theta-1)
  dx[k] <- a * dcauchy(x[k], location=tL*(1+p), scale=c*tL) / A
  
  return(dx)
}

# probability function
#' @rdname calibrations
#' @export
pL <- Vectorize(
  function(q, tL, p=0.1, c=1, pL=0.025) {
    integrate(dL, lower=0, upper=q, tL=tL, p=p, c=c, pL=pL)$value
  }
)

# quantile function
#' @rdname calibrations
#' @export
qL <- Vectorize(
  function(prob, tL, p=0.1, c=1, pL=0.025) {
    uniroot(function(x) (pL(x, tL=tL, p=p, c=c, pL=pL) - prob), 
            int = c(0, tL * 10), extendInt="upX")$root
  }
)

# TODO: Add functions pB, qB, pU, and qU

# Joint bounds: (see Yang and Rannala 2006)
# We use exponentials here for easier plotting
#' @rdname calibrations
#' @export
dB <- function(x, tL, tU, pL=.025, pU=.025) {
  h <- (1-pL-pU) / (tU - tL)
  l2 <- h/pU
  il <- (x < tL)
  iu <- (x > tU)
  y <- rep(h, length(x))
  theta <- h * tL / pL
  y[which(il)] <- pL * theta/tL * (x[il]/tL)^(theta-1)
  y[which(iu)] <- pU * dexp(x[iu] - tU, l2)
  return (y)
}

# Maximum bound: (see Yang and Rannala 2006)
#' @rdname calibrations
#' @export
dU <- function(x, tU, pU=0.025) {
  dB(x=x, tL=0, tU=tU, pL=0, pU=pU)
}

# # This file originates from the phylogenomic mammal paper (dos Reis et al. 2012, PRSB)
# # TODO: Export functions and add documentation
# 
# u975 <- function(tL, p=0.1, c=1) {
#   A <- .5 + atan(p/c) / pi
#   return (tL*(1 + p + c*tan(pi*(.5 - .025*A/0.975))))
# }
# 
# # rename function
# c95 <- u975
# 
# ## c95(1, .5, c(.2, .5, 1, 2))
# 
# ## # Fossil calibrations:
# ## # Basal marsupial:
# ## c95(.615, .1, 1) # 15.03
# 
# ## # Paenungulate (elephant, manati):
# ## c95(.556, .1, 1) # 13.58
# 
# ## c95(.615, .1, 1) # 15.03
# ## c95(.524, .1, 1) # 12.8
# ## c95(.625, .1, 1) # 15.3
# ## c95(.486, .1, 1) # 11.9
# ## c95(.556, .1, 1) # 13.6
# ## c95(.337, .1, 1) # 8.2
# ## c95(.0725, .1, 1) # 1.77
# ## c95(.484, .1, 1) # 11.82

# # cumulative probability function:
# pL <- function(q, tL, p=0.1, c=1, b=.025) {
#   # TODO: redo using cauchy cumulative function F(x)
#   return (integrate(dL, lower=0, upper=q, tL=tL, p=p, c=c, b=b))
# }

