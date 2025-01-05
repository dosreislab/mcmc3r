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
#' @details 
#' Calculates the density, distribution and quantile functions for the minimum,
#' \code{dL}, joint, \code{dB}, and maximum, \code{dU}, calibration bounds as
#' implemented in MCMCtree (Yang and Rannala, 2006; Inoue et al. 2010). The
#' minimum bound is implemented using a truncated Cauchy distribution (Inoue et
#' al. 2010).
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
#' # Cumulative distribution:
#' curve(pL(x, 1), from=0, to=10, n=5e2)
#' 
#' # Plot a joint bounds calibration density:
#' curve(dB(x, 1, 6), from=0, to=10, n=5e2)
#' # Cummulative distribution:
#' curve(pB(x, 1, 6), from=0, to=10, n=5e2)
#' 
#' # Plot a maximum bound calibration density:
#' curve(dU(x, 6), from=0, to=10, n=5e2)
#' # Cummulative distribution:
#' curve(pU(x, 6), from=0, to=10, n=5e2)
#' 
#' # Check quantile function for minimum bound (or truncated-Cauchy):
#' qv <- 0:20; pvL <- pL(qv, tL=1)
#' # calculate quantiles back from probability vector:
#' # (note numerical error)
#' plot(qv, qL(pvL, tL=1)); abline(0, 1)
#' 
#' # Check quantile function for joint bounds:
#' pvB <- pB(qv, tL=2, tU=10, pL=.02, pU=.1)
#' # calculate quantiles back:
#' plot(qv, qB(pvB, tL=2, tU=10, pL=.02, pU=.1)); abline(0, 1)
#' 
#' # Check quantile function for upper bound:
#' pvU <- pU(qv, tU=15, pU=.15)
#' # calculate quantiles back:
#' plot(qv, qU(pvU, tU=15, pU=.15)); abline(0, 1)
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
  A <- .5 + atan(p/c) / pi # this is 1 - F(x)
  theta <- a/pL * 1/(pi*A*c*(1 + (p/c)^2))
  
  dx <- numeric(length(x))
  
  i <- (x < 0)
  j <- (0 <= x & x < tL)
  k <- (x >= tL)
  
  dx[i] <- 0
  dx[j] <- pL * theta/tL * (x[j]/tL)^(theta-1)
  dx[k] <- a * dcauchy(x[k], location=tL*(1+p), scale=c*tL) / A
  
  return(dx)
}

#' @rdname calibrations
#' @export
pL <- function(q, tL, p=0.1, c=1, pL=0.025) {
  a <- 1 - pL
  A <- .5 + atan(p/c) / pi
  theta <- a/pL * 1/(pi*A*c*(1 + (p/c)^2))
  
  px <- numeric(length(q))
  
  i <- (q < 0)
  j <- (0 <= q & q <= tL)
  k <- (q > tL)
  
  px[i] <- 0
  px[j] <- pL * (q[j] / tL)^theta
  px[k] <- pL + a * (pcauchy(q[k], location=tL*(1+p), scale=c*tL) - 1 + A) / A
  
  return(px)
}

# probability function
# @rdname calibrations
# @export
# old.pL <- Vectorize(
#   function(q, tL, p=0.1, c=1, pL=0.025) {
#     integrate(dL, lower=0, upper=q, tL=tL, p=p, c=c, pL=pL)$value
#   }
# )

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
  theta <- h * tL / pL
  l2 <- h/pU
  
  i0 <- (x <= 0)
  il <- (0 < x & x < tL)
  iu <- (x > tU)
  
  y <- numeric(length(x))
  y[!(i0 | il | iu)] <- h
  y[il] <- pL * theta/tL * (x[il]/tL)^(theta-1)
  y[iu] <- pU * dexp(x[iu] - tU, l2)
  
  return (y)
}

#' @rdname calibrations
#' @export
pB <- function(q, tL, tU, pL=.025, pU=.025) {
  h <- (1-pL-pU) / (tU - tL)
  theta <- h * tL / pL
  l2 <- h/pU
  
  i <- (q < 0)
  j <- (0 <= q & q < tL)
  k <- (tL <= q & q <= tU)
  l <- (q > tU)
  
  px <- numeric(length(q))
  px[i] <- 0
  px[j] <- pL * (q[j] / tL)^theta
  px[k] <- (q[k] - tL) * h + pL
  px[l] <- 1 - pU + pU * pexp(q[l] - tU, rate=l2)
  
  return(px)
}

#' @rdname calibrations
#' @export
qB <- Vectorize(
  function(prob, tL, tU, pL=0.025, pU=0.025) {
    uniroot(function(x) (pB(x, tL=tL, tU=tU, pL=pL, pU=pU) - prob), 
            int = c(0, tU * 2), extendInt="upX")$root
  }
)

# Maximum bound: (see Yang and Rannala 2006)
#' @rdname calibrations
#' @export
dU <- function(x, tU, pU=0.025) {
  dB(x=x, tL=0, tU=tU, pL=0, pU=pU)
}

#' @rdname calibrations
#' @export
pU <- function(q, tU, pU=0.025) {
  pB(q=q, tL=0, tU=tU, pL=0, pU=pU)
}

#' @rdname calibrations
#' @export
qU <- function(prob, tU, pU=0.025) {
  qB(prob=prob, tL=0, tU=tU, pL=0, pU=pU)
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

