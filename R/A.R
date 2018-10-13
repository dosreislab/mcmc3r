#' A matrix
#'
#' Matrix which, once multiplied by its transpose,
#' yields to the inverse of the estimate of the shrinkage correlation
#' matrix, R.sh. The latter is used to transform the data set as
#' it corrects the morphological data set for the corresponding
#' character correlation
#'
#' @format A matrix of size n x n, where n = 87
#' (morphological traits: 29 landmarks x 3D coordinates):
#' \describe{
#'   \item{n}{Number of traits for which the correlation values
#'   have been calculated, 87}
#' }
"A"
