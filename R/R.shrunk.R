#' Shrunk correlation matrix estimate
#'
#' Shrunk estimate of the correlation matrix obtained after applying
#' the equation \eqn{{R}^{*} = \delta I + (1-\delta){R'}}{R.sh = delta I +
#' (1-delta)R.unb}
#' The shrinkage parameter used to estimate the matrix is
#' is delta = 0.01 so the resulting estimate is closer to the
#' unbiased correlation matrix.
#'
#' @format A matrix of size n x n, where n = 87
#' (morphological traits, the 87 coordinates):
#' \describe{
#'   \item{n}{Number of coordinates for which the correlation values
#'   have been calculated, 87}
#' }
"R.shrunk"
