# This file originates from the phylogenomic mammal paper (dos Reis et al. 2012, PRSB)

u975 <- function(tL, p=0.1, c=1) {
  A <- .5 + atan(p/c) / pi
  return (tL*(1 + p + c*tan(pi*(.5 - .025*A/0.975))))
}

# rename function
c95 <- u975

## c95(1, .5, c(.2, .5, 1, 2))

## # Fossil calibrations:
## # Basal marsupial:
## c95(.615, .1, 1) # 15.03

## # Paenungulate (elephant, manati):
## c95(.556, .1, 1) # 13.58

## c95(.615, .1, 1) # 15.03
## c95(.524, .1, 1) # 12.8
## c95(.625, .1, 1) # 15.3
## c95(.486, .1, 1) # 11.9
## c95(.556, .1, 1) # 13.6
## c95(.337, .1, 1) # 8.2
## c95(.0725, .1, 1) # 1.77
## c95(.484, .1, 1) # 11.82

# L bound calibration:
# density function:
dL <- function(x, tL, p=0.1, c=1, b=0.025) {
  a <- 1 - b
  A <- .5 + atan(p/c) / pi
  theta <- a/b * 1/(pi*A*c*(1 + (p/c)^2))
  
  dx <- numeric(length(x))
  
  i <- which(x < 0)
  j <- which(0 <= x & x < tL)
  k <- which(x >= tL)
  
  dx[i] <- 0
  dx[j] <- b * theta/tL * (x[j]/tL)^(theta-1)
  dx[k] <- a * dcauchy(x[k], location=tL*(1+p), scale=c*tL) / A
  
  return(dx)
}

# cumulative probability function:
pL <- function(q, tL, p=0.1, c=1, b=.025) {
  # redo using cauchy cumulative function F(x)
  return (integrate(dL, lower=0, upper=q, tL=tL, p=p, c=c, b=b))
}

# Joint bounds: (see Yang and Rannala 2006)
# We use exponentials here for easier plotting
dB <- function(x, tL, tU, pL=.025, pU=.025) {
  h <- (1-pL-pU) / (tU - tL)
  l1 <- h/pL
  l2 <- h/pU
  il <- (x < tL)
  iu <- (x > tU)
  y <- rep(h, length(x))
  y[which(il)] <- pL * dexp(-x[il] + tL, l1)
  y[which(iu)] <- pU * dexp(x[iu] - tU, l2)
  return (y)
}