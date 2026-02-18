#' Calculate the approximate log-likelihood
#' 
#' @param b numeric, the branch length
#' @param b.mle numeric, the branch length MLE
#' @param g.mle numeric, the gradient at the branch length MLE
#' @param H.mle numeric, the second derivative at the branch length MLE
#' @param transform character, the transform to be used (see details)
#'
#' @details Calculates the approximate log-likelihood for two species using
#'   Taylor expansion of the log-likelihood around the maximum likelihood
#'   estimate (MLE) of the branch length, using the branch length transforms
#'   described in dos Reis and Yang (2011).
#'  
#' @examples
#' # This reproduces Fig. 1 in dos Reis and Yang (2011):
#' 
#' # Exact likelihood for two-species function under the JC69 model:
#' JC69.lnL <- function(b, x, n) {
#'   p <- .75 - .75 * exp(-4 * b / 3)
#'   x * log(p) + (n - x) * log (1 - p)
#' }
#' 
#' # Exact gradient under the JC69 model:
#' JC69.g <- function(b, x, n) {
#'   p <- .75 - .75 * exp(-4 * b / 3)
#'   (x / p - (n - x) / (1 - p)) * exp(-4 * b / 3)
#' }
#' 
#' # Exact Hessian under the JC69 model:
#' JC69.H <- function(b, x, n) {
#'   p <- .75 - .75 * exp(-4 * b / 3)
#'   (-x/p^2 - (n-x)/(1-p)^2) * exp(-8 * b/3) - 4/3 * (x/p - (n-x)/(1-p)) * exp(-4 * b/3)
#' } 
#' 
#' # MLE of branch length under the JC69.modle:
#' JC69.b <- function(x, n) {
#'   -.75 * log(1 - 4/3 * x/n)
#' }
#' 
#' par(mfrow=c(1,3))
#' # Case 1: MLE of b is zero
#' x. <- 0; n <- 100; 
#' b <- JC69.b(x., n); lnL.max <- x. * log(1e-6/n) + (n - 1e-6) * log(1 - 1e-6/n)
#' curve(JC69.lnL(x, x., n) - lnL.max, from=0, to=.6, ylim=c(-50, 0), xlab="b", ylab="lnL", las=1)
#' g <- JC69.g(1e-6, x., n); H <- JC69.H(1e-6, x., n)
#' curve(lnL.app(x, 0, g, H, "NT"), lty=2, add=TRUE)
#' curve(lnL.app(x, 0, g, H, "SQRT"), lty=2, add=TRUE)
#' 
#' # Case 2: MLE of b is between zero and infinity
#' x. <- 37; n <- 100
#' b <- JC69.b(x., n); lnL.max <- x. * log(x./n) + (n - x.) * log(1 - x./n)
#' g <- JC69.g(b, x., n); H <- JC69.H(b, x., n)
#' curve(JC69.lnL(x, x., n) - lnL.max, from=0, to=2, ylim=c(-50, 0), xlab="b", ylab="lnL", las=1)
#' curve(lnL.app(x, b, g, H, "NT"), lty=2, add=TRUE)
#' curve(lnL.app(x, b, g, H, "SQRT"), lty=2, add=TRUE)
#' abline(v=b, lty=3)
#' 
#' # Case 3: MLE of b is very large
#' x. <- 74; n <- 100; 
#' b <- JC69.b(x., n); lnL.max <- x. * log(x./n) + (n - x.) * log(1 - x./n)
#' curve(JC69.lnL(x, x., n) - lnL.max, from=0, to=4, ylim=c(-50, 0), xlab="b", ylab="lnL", las=1)
#' g <- JC69.g(b, x., n); H <- JC69.H(b, x., n)
#' curve(lnL.app(x, b, g, H, "NT"), lty=2, add=TRUE)
#' curve(lnL.app(x, b, g, H, "SQRT"), lty=2, add=TRUE)
#' abline(v=b, lty=3)
#'
#' @references dos Reis and Yang (2011) Approximate likelihood calculation on a
#' phylogeny for Bayesian estimation of divergence times. \emph{Molecular
#' Biology and Evolution}, 28: 2161–2172.
#'
#' @export
lnL.app <- function(b, b.mle, g.mle, H.mle, transform=c("NT", "SQRT", "ARCSIN")) {
  transform <- match.arg(transform)
  if (transform == "NT") {
    u <- b; u.mle <- b.mle
    gu.mle <- g.mle; Hu.mle <- H.mle
  } else if (transform == "SQRT") {
    u <- sqrt(b); u.mle <- sqrt(b.mle)
    gu.mle <- 2 * g.mle * u.mle; Hu.mle <- 2 * g.mle + 4 * H.mle * b.mle
  } else if (transform == "ARCSIN") {
    u <- 2 * arcsin(sqrt(.pf(b))); u.mle <- 2 * arcsin(sqrt(.pf(b.mle)))
    db.du <- cos(u.mle/2) * sin(u.mle/2) / (1 - 4 * sin(u.mle/2)^2 / 3)
    d2b.du2 <- .5 * (cos(u.mle/2)^2 - sin(u.mle/2)^2) / (1 - 4 * sin(u.mle/2)^2 / 3) +
      4 * cos(u.mle/2)^2 * sin(u.mle/2)^2 / (3 * (1 - 4 * sin(u.mle/2)^2 / 3))
    gu.mle <- g.mle * db.du
    Hu.mle <- g.mle * d2b.du2 + H.mle * db.du^2
  }
  (u - u.mle) * gu.mle + .5 * (u - u.mle)^2 * Hu.mle
}

# Internal: calculate proportion of different sites given JC69 model
.pf <- function(b.) .75 - .75 * exp(-4 * b. / 3)