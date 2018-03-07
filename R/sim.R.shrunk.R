#' Estimate of the shrunk correlation matrix for a population
#'
#' This estimate has been calculated from the previously
#' generated object "sim.population" as the population
#' sample (20 specimens).
#' We used the function calc.pop.cor to generate it and selected
#' the item "R.shrunk" from the list returned.
#'
#' @format A matrix of size n x n, where n = 100
#' (morphological traits, the 100 simulated traits):
#' \describe{
#'   \item{n}{Number simulated traits for which the correlation values
#'   have been calculated, 100}
#' }
"sim.R.shrunk"
