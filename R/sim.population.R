#' Simulated population
#'
#' Simulated population with a sample of s = 20 specimens
#  under a normal distribution with mean = 0, variance c = 0.25,
#  i.e. x ~ N(0,0.25), and R = sim.R (constant correlation model with
#  rho = 0.50 and n = 100 characters, thus size is n x n).
#'
#' @format A matrix of size n x n, where n = 100
#' (morphological traits, the 100 simulated traits):
#' \describe{
#'   \item{n}{Number of simulated traits for which the correlation values
#'   have been calculated, 100}
#' }
"sim.population"
