#' Convert an MCMC sample from MCMCTree or BPP to a list of trees
#'
#' @param tree an object of class phylo
#' @param mcmc data frame with an MCMC sample from MCMCTree or a BPP A00
#'   analysis
#' @param time.name character vector of length one, this is usually
#'   \code{"tau_"} for BPP, or \code{"t_"} for MCMCTree.
#' @param thin numeric, the fraction of MCMC samples to keep
#'
#' @details \code{tree} must be rooted and strictly bifurcating, and it must
#'   match (i.e, it must use the same Newick representation) as the tree used by
#'   BPP or MCMCTree to obtain the MCMC sample. The function uses the node ages
#'   in \code{mcmc} to calculate branch lengths and generate a list of trees
#'   (with the same topology as \code{tree}), one tree per (thinned) MCMC
#'   sample. The tips of the phylogeny are assumed to have age zero.
#'
#' @return An object of class multiPhylo (i.e., a list of trees).
#'
#' @examples
#' # hominid MCMC sample from BPP:
#' data(hominids)
#' hominid.trees <- mcmc2multiphylo(hominids$tree, hominids$mcmc, "tau_", 0.001)
#'
#' \dontrun{
#' # If you have the ape package installed, you can output the trees in 
#' # Newick format
#' ape::write.tree(hominid.trees) 
#' }
#'
#' @author Mario dos Reis
#'
#' @export
# TODO: Implement tips of different ages
# TODO: Verify that ape always reads the trees in a way that is compatible with
# MCMCtree's
mcmc2multiphylo <- function(tree, mcmc, time.name, thin) {
  tts <- list()
  N <- dim(mcmc)[1]
  ti <- grep(time.name, names(mcmc))
  ii <- floor(seq(from=1, to=N, length.out = N * thin))
  n <- length(ii)

  for (i in 1:n) {
    tts[[i]] <- .mcmc2multiphylo.treef(tree, unlist(mcmc[ii[i],ti]))
  }
  class(tts) <- "multiPhylo"

  return (tts)
}

# This helper function does pre-order tree traversal to calculate branch lengths
# from node ages.
.mcmc2multiphylo.treef <- function(tree, ages) { # TODO: Tip ages
  # tree is assummed rooted and perfectly bifurcating
  ns <- tree$Nnode + 1
  root <- ns + 1
  tip.ages <- rep(0, ns) # assumes tips are extant
  all.ages <- c(tip.ages, ages)
  blens <- numeric(2*ns - 2)

  porder <- function(node) { # parent is NULL if root
    #print(paste("i'm", node))
    desc <- which(tree$edge[,1] == node) # my descendant branches

    if (node != root) {
      parent <- tree$edge[which(tree$edge[,2] == node)]
      mybranch <- which(tree$edge[,2] == node)
      blens[mybranch] <<- all.ages[parent] - all.ages[node]
      #print(blens)
    }

    # If I'm a tip:
    if (!length(desc)) { return(0) }
    # Otherwise I'm an internal node:
    else {
      porder(tree$edge[desc[1], 2]) # visit my left daughter
      porder(tree$edge[desc[2], 2]) # visit my right daughter
    }

    return(0)
  }

  porder(root)

  tt <- tree
  tt$edge.length <- blens

  return(tt)
}
