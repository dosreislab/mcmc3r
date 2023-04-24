#' Ancestral character reconstruction from an MCMC sample
#' 
#' @description  
#' Obtain the ancestral reconstruction of characters from an MCMC sample using
#' Felsenstein (1973) model of continuous character evolution.
#' 
#' @param tree a rooted, strictly bifurcating phylogeny
#' @param M s x k matrix of k continuous morphological measurements for s species
#' @param mcmc data frame with MCMC output from MCMCtree
#' @param time.name character vector of length one
#' @param rate.name character vector of length one
#' @param tip.ages numeric, the ages of the terminal taxa in the tree
#' 
#' @details
#' 
#' The function first calculates the mean of divergence times and morphological
#' rates in the MCMC sample. These are used to reconstruct the branch lengths in
#' units of morphological evolution, and then Eq. (7) in Felsenstein (1973) is
#' used to calculate the ancestral reconstruction on the tree. 
#'
#' Note \code{time.name} is the name format used for the node ages in the MCMC
#' dataframe, usually of the form \code{time.name = "t_"}. Similarly
#' \code{rate.name} is the name format used in the MCMC sample for the rates, of
#' the form \code{rate.name = "r_g1_"} where the subscript in "g" must be the
#' partition number containing the morphological rates. Note \code{tree} must be
#' the same used by MCMCtree when generating the MCMC sample, right down to its
#' Newick representation. Taxon names in \code{tree} and \code{M} must match.
#' 
#' @author Mario dos Reis
#' 
#' @seealso \link{write.morpho}
#' 
#' @references 
#' Felsenstein J (1973) Maximum-likelihood estimation of evolutionary trees from
#' continuous characters. \emph{Am J Hum Genet,} 25: 471--492.
#' 
#' @return 
#' A n x k matrix with the ancestral reconstruction for the k characters at the
#' n internal nodes of the phylogeny.
#' 
#' @examples 
#' data(carnivores) 
#' back.ages = as.numeric(unlist(strsplit(carnivores$tree$tip.label, "\\^"))[seq(from=2, to=2*19, by=2)])
#' back.ages = max(back.ages) - back.ages 
#' C <- carnivores$C.proc
#' rownames(C) <- rownames(carnivores$M)
#' recM = mcmc2anc(carnivores$tree, C, mcmc=carnivores$mcmc, time.name="t_", 
#'        rate.name="r_g1_", tip.ages=back.ages)
#'
#' x <- seq(1, 87, by=3); y <- seq(2, 87, by=3)
#' # Plot landmarks for the 19 carnivores:
#' plot(carnivores$C.proc[,x], carnivores$C.proc[,y], pch='+', cex=.5)
#' # plot ancestral reconstruction at the root (node 20):
#' points(recM["20",x], recM["20",y], pch=19, col="red")
#' # mean shape
#' mS <- apply(carnivores$C.proc, 2, mean)
#' points(mS[x], mS[y], pch=19, col="blue")
#' 
#' @export
mcmc2anc <- function(tree, M, mcmc, time.name, rate.name, tip.ages=NULL) {
  # tt must be rooted and strictly bifurcating
  # number of branches and species: nb = 2*s - 2 -> s = nb/2 + 1
  ns <- tree$Nnode + 1
  N <- nrow(mcmc)
  
  pi <- match(tree$tip.label, rownames(M))
  if (any(is.na(pi))) stop("rownames in M must be the same as in tree$tip.label")
  M <- M[pi,]; M <- as.matrix(M)
  
  # set tip ages
  if (is.null(tip.ages)) tip.ages <- rep(0, ns)
  
  # find node ages and morphological rates
  ti <- grep(time.name, names(mcmc))
  ri <- grep(rate.name, names(mcmc))
  
  # construct ages matrix
  ages <- cbind(matrix(tip.ages, nrow=N, ncol=length(tip.ages), byrow=TRUE), mcmc[,ti])
  
  parents <- tree$edge[,1]
  daughters <- tree$edge[,2]
  
  blen <- (ages[,parents] - ages[,daughters]) * mcmc[,ri]
  blen <- apply(blen, 2, mean)
  
  daughters[daughters > ns] <- daughters[daughters > ns] - 1
  tree$edge.length <- blen[daughters]
  
  # reconstruct using postorder traversal
  recM <- .ancrec(tree, M)
  
  return(recM)
}

# Calculate transformation of trait values, x_i, and branch lengths, v_i,
# as in Felsenstein (1973)
.tin <- function(x1, x2, v1, v2, v3) {
  xp <- (v2 * x1 + v1 * x2) / (v1 + v2)
  vp <- v3 + v1 * v2 / (v1 + v2)
  # identify whether there is missing data in x1 or x3 and act accordingly
  #if (is.na(sum(x1))) { i <- is.na(x1); vp[i] <- v2 + v3; xp[i] <- x2 }
  #if (is.na(sum(x2))) { i <- is.na(x2); vp[i] <- v1 + v3; xp[i] <- x1 }
  return (list(xp=xp, vp=vp)) # x' and v'
}

# Quick reconstruction using postorder tree traversal
.ancrec <- function(tree, trait_matrix) {
  .contml <- function(node, tree, trait_matrix, vP=0) {
    # print (node)
    # tree$edge: nodes forming the edges
    # tree$edge.length: branch lengths
    desc <- which(tree$edge[,1] == node) # my descendants
    myself <- which(tree$edge[,2] == node) # myself
    
    # If I'm a tip do:
    if (!length(desc)) {
      return (list(xp=trait_matrix[node,], vp=tree$edge.length[myself] + vP))
    }
    
    # If I'm an internal node (including the root as the tree is rooted) do:
    if (length(desc) == 2) {
      l1 <- .contml (tree$edge[desc[1], 2], tree, trait_matrix, vP)
      l2 <- .contml (tree$edge[desc[2], 2], tree, trait_matrix, vP)
      # my (untransformed) branch length
      v3 <- tree$edge.length[myself]
      # get transform
      l3 <- .tin (l1$xp, l2$xp, l1$vp, l2$vp, v3) # transform
      # debugging
      #print(ll - l1$ll - l2$ll)
      recM[node - length(tree$tip.label),] <<- l3$xp 
      return (list(xp=l3$xp, vp=l3$vp))
    }
    
    # I shouldn't be here:
    stop ("Something went wrong, too many children!")
  }
  # matrix with ancestral reconstruction
  recM <- matrix(nrow = tree$Nnode, ncol = ncol(trait_matrix))
  rownames(recM) <- (tree$Nnode + 2):(2 * tree$Nnode + 1)
  
  # call recursive function for reconstruction
  .contml(length(tree$tip.label) + 1, tree, trait_matrix)
  
  return(recM)
}

# Landmark array configuration:
# p x k x n
# p: number of points
# k: point dimension
# n: number of specimens

# M %*% L %*% L' %*% M'