# Functions for preparing and parsing MCMCTree files for
# Bayes factor calculations

# Calculate the beta values for Bayes Factor calculation
# Points are chosen according to Gauss-Legendre quadrature rules
# n: number of beta points
#' @export
BFbeta <- function(n) {
  qr <- glqrules[[n]]
  b <- (1 + qr$x) / 2
  return(list(x=qr$x, b=b, w=qr$w))
}

# Prepare mcmctree control files with appropriate beta values
# for Bayes Factors calculations
#' @export
make.bfctlf <- function(file="mcmctree.ctl", n) {
  df <- BFbeta(n)
  betaf <- "beta"
  for (i in 1:n) {
    dir.create(as.character(i))
    cat(paste("BayesFactorBeta = ", df$b[i], "\n", sep=""), file=betaf)
    newf <- paste(i, "/", file, sep="")
    file.append(file1=newf, file2=c(file, betaf))
  }
  unlink(betaf)
  write.table(df, row.names=FALSE, file="beta.txt")
  # TODO: print w, x, b to a file
  # TODO: add option for replicate runs
}

# Gauss-Legendre quadrature
#' @export
glq <- function(logml, b.df) sum(logml * b.df$w) / 2

# Bayes factors
# ml: marginal likelihood of models. length(ml) > 1
# returns BF and posterior Pr (under uniform prior)
#' @export
BF <- function(logml) {
  if (length(logml) < 2) stop("Provide at least two log-marginal likelihoods.")
  logml0 <- max(logml)
  lbf <- logml - logml0
  bf <- exp(lbf)
  sbf <- sum(bf)
  return(list(BF=bf, Pr=bf/sbf))
}
