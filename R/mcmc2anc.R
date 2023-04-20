#' Fast ancestral character reconstruction from an MCMC sample
#' TODO: Implement covariance transform and do proper example
#' @examples 
#' tt = ape::read.tree("~/phyl/carnivores/tutorial/carnivores.tree")
#' Mdf <- read.table("~/phyl/carnivores/tutorial/M.txt")
#' M <- Mdf[,-1]; rownames(M) <- Mdf[,1]
#' mcmc = read.table("~/phyl/carnivores/tutorial/mcmc.txt", head=TRUE)
#' back.ages = as.numeric(unlist(strsplit(tt$tip.label, "\\^"))[seq(from=2, to=2*19, by=2)])
#' back.ages = max(back.ages) - back.ages
#' xx = mcmc2anc(tt, M, mcmc=mcmc, time.name="t_", rate.name="r_g1_", tip.ages=back.ages)
#' tr <-  ape::read.tree("~/phyl/carnivores/tutorial/rate_g1.tree")
#' plot(xx$rates, tr$edge.length)
#' @export
mcmc2anc <- function(tree, M, mcmc, time.name, rate.name, tip.ages) {
  # tt must be rooted and strictly bifurcating
  # number of branches and species: nb = 2*s - 2 -> s = nb/2 + 1
  ns <- tree$Nnode + 1
  #nb <- 2 * ns - 2
  #inodes <- (ns+1):(nb-1)
  
  pi <- match(tree$tip.label, rownames(M))
  M <- M[pi,]; M <- as.matrix(M)
  
  # set tip ages
  if (is.null(tip.ages)) tip.ages <- rep(0, ns)
  
  # find node ages and morphological rates
  ti <- grep(time.name, names(mcmc))
  ri <- grep(rate.name, names(mcmc))
  
  # 0%, 25%, 50%, 75% and 100% quantiles
  tm <- apply(mcmc[,ti], 2, mean)
  rm <- apply(mcmc[,ri], 2, mean)
  
  # construct ages matrix
  ages <- c(tip.ages, tm)
  
  parents <- tree$edge[,1]
  daughters <- tree$edge[,2]
  
  # Delta time (time span of branches)
  Dt <- ages[parents] - ages[daughters]
  blen <- Dt * rm
  
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