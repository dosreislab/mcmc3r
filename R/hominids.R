#' A BPP A00 MCMC sample for an hominid phylogeny
#'
#' @description This dataset contains the results from the BPP A00 analysis of
#'   hominid evolution from Angelis and dos Reis (2015).
#'
#' @format \code{hominids} is a list with elements \code{mcmc}, a dataframe with
#'   20,000 rows and 8 columns, and \code{tree}, an object of class
#'   \code{phylo} from the \code{ape} package.
#'
#'   \code{mcmc} is a posterior sample from a BPP A00 MCMC analysis containing
#'   the relative divergence times (tau's) and nucleotide diversities (theta's)
#'   for the four species ape (hominid) phylogeny.
#'
#'   \code{tree} contains the phylogeny with node ages given as the posterior
#'   means of the tau's in \code{mcmc}.
#'
#' @source K. Angelis and M. dos Reis (2015) \emph{The impact of ancestral
#'   population size and incomplete lineage sorting on Bayesian estimation of
#'   species divergence times.} Curr. Zool., 61: 874--885.
#'
#' @seealso \code{\link{microcebus}}
"hominids"
