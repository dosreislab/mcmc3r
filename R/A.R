#' A matrix
#'
#' A matrix which, once multiplied by its transpose,
#  yields to the inverse of the estimate of the shrunk correlation
#  matrix. This matrix is later used to transform the data set as
#  it corrects the morphological data set for the corresponding
#  correlation
#'
#' @format A matrix of size n x n, where n = 87
#' (morphological traits, the 87 coordinates):
#' \describe{
#'   \item{n}{Number of coordinates for which the correlation values
#'   have been calculated, 87}
#' }
"A"
