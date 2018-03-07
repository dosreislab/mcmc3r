#' Estimate the population variance of a population sample.
#'
#' For this estimate, we used the previously generated object
#' "sim.population" as the population sample (20 specimens).
#' We used the function calc.pop.cor, which generates the estimate
#' of the shrunk correlation matrix but additionally returns the
#' variance of the population matrix. Note that these two variables
#' are to be used within the write.morpho function when generating
#' the alignment file to be input in MCMCTree.
#'
#' @format A vector of size n, where n = 100
#' (morphological traits, the 100 simulated traits), with the
#' estimated population variances:
#' \describe{
#'   \item{n}{Number of estimated population variances, 100}
#' }
"sim.var.pop"
