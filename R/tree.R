#' Simulated 8-species tree
#'
#' Simulated species tree with 8 species (class "phylo"). The function
#' \code{\link[ape]{read.tree}} was used to generate this object.
#'
#' @format An object of class "phylo" with a simulated species tree with
#' s = 8 species. It contains the following components:
#' \describe{
#'   \item{edge}{Two-column matrix with 14 rows, where every row
#'   is one edge in the tree, the first column is the ancestor node
#'   and the second column its daughter node}
#'   \item{edge.length}{A numeric vector with the branch lengths of
#'    the tree}
#'   \item{Nnode}{Numeric, the number of (internal) nodes}
#'   \item{tip.label}{A vector with the names of the tips, class
#'   "character"}
#' }
"tree"
