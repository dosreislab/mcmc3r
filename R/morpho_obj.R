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

#' 29 3D landmarks from the skulls of 19 carnivoran specimens before Procrustes analysis
#'
#' A 3D array containing the 29 3D landmarks collected from the skulls of 19 carnivoran
#' specimens before carrying out a Procrustes analysis (PA).
#' Please take a look at the description in morpho/data-raw/C.R to understand
#' how this object was generated.
#'
#' @format An array with k = 29 (landmarks), q = 3 (coordinates) and s = 19 (specimens):
#' \describe{
#'   \item{k}{landmark points collected from 19 carnivoran specimens, 29}
#'   \item{q}{coordinates in 3D or 2D, 3}
#'   \item{s}{number of specimens, 18}
#' }
"C.arr.unal"

#' 29 3D landmarks from the skulls of 19 carnivoran specimens before Procrustes analysis
#'
#' A matrix containing the 29 3D landmarks collected from the skulls of 19 carnivoran
#' specimens before carrying out a Procrustes analysis (PA). Please take a look at the
#' description in morpho/data-raw/C.R to understand
#' how this object was generated.
#'
#' @format A matrix with s = 19 rows and n = 87 columns (87/3 = 29 landmarks):
#'   \describe{
#'   \item{s}{Rows, specimens from which landmarks were collected, 19}
#'   \item{n}{Columns, 87 traits (29 landmarks in 3D) after the PA}
#' }
"C.mat.unal"

#' Object of class procSym output by Morpho after PA
#'
#' Object of class procSym output by Morpho, which mean shape is later
#' used to align the foxes ("Vulpes vulpes") specimens to it. Please take a
#' look at the description in morpho/data-raw/C.R to understand
#' how this object was generated.
#'
#' @format Object procSym
#'   \describe{
#'   \item{...}{Check \code{\link[Morpho]{procSym}} for more details}
#' }
"C.PS"

#' 29 3D landmarks from the skulls of 19 carnivoran specimens after Procrustes analysis
#'
#' A matrix containing the 29 3D landmarks collected from the skulls of 19 carnivoran
#' specimens after carrying out a Procrustes analysis (PA).
#' Please take a look at the description in morpho/data-raw/C.R to understand
#' how this object was generated. This is the morphological alignment.
#'
#' @format A matrix with s = 19 rows and n = 87 columns (87/3 = 29 landmarks):
#'   \describe{
#'   \item{s}{Rows, specimens from which landmarks were collected, 19}
#'   \item{n}{Columns, 87 traits (29 landmarks in 3D) after the PA}
#' }
"C"

#' 29 3D landmarks from the skulls of 19 carnivoran specimens
#'
#' A dataset containing the 29 3D landmarks collected from the skulls of 19 carnivoran specimens
#' This data.frame consists of a first column with the specimen labels used by MCMCtree,
#' a second column with the voucher names of the specimens collected, and then 87 columns
#' with the coordinates for each trait (29 landmarks x 3D coordinates).
#' Please take a look at the description in morpho/data-raw/carnivores19x29.raw.R to understand
#' how this object was generated.
#'
#' @format A data.frame with s = 19 rows and p = 89 columns ( 2 info columns + 87 traits ):
#' \describe{
#'   \item{s}{Rows, specimens from which landmarks were collected, 19}
#'   \item{p}{Columns, specimens information in 1st and 2nd column + 87 traits
#'    (29 landmarks in 3D)}
#' }
"carnivores19x29.raw"

#' Estimated shrinkage correlation matrix
#'
#' Estimated shrinkage correlation matrix obtained after using the
#' \code{corpcor::cor.shrink} package.
#'
#' @format A matrix of size n x n, where n = 87
#' (morphological traits, 29 landmarks x 3D coordinates):
#' \describe{
#'   \item{n}{Number of traits for which the correlation values
#'   have been calculated, 87}
#' }
"R.sh"

#' Correlation matrix for simulations
#'
#' True correlation matrix simulated to be used in the examples detailed in the \code{sim.pop()}
#' function. The matrix follows the constant correlation model, hence all values
#' outside the diagonal are rho = 0.50. The size is \code{p x p}, being n = 100
#' the number of characters.
#'
#' @format A matrix of size n x n, where n = 100
#' (morphological traits, the 100 simulated continuous traits):
#' \describe{
#'   \item{n}{Number of simulated continuous traits for which the correlation values
#'   have been calculated, 100}
#' }
"sim.R"

#' Simulated 8-species tree
#'
#' Simulated 8-species tree of class "phylo". The function
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
"sim.tree"

#' 29 3D landmarks from the skulls of 21 Vulpes vulpes specimens before Procrustes analysis
#'
#' A 3D array containing the 29 3D landmarks collected from the skulls of 21 "Vulpes vulpes"
#' specimens before carrying out a Procrustes analysis (PA).
#' Please take a look at the description in morpho/data-raw/V.R to understand
#' how this object was generated.
#'
#' @format An array with k = 29 (landmarks), q = 3 (coordinates) and s = 21 (specimens):
#' \describe{
#'   \item{k}{landmark points collected from 21 foxes specimens, 29}
#'   \item{q}{coordinates in 3D or 2D, 3}
#'   \item{s}{number of specimens, 21}
#' }
"V.arr.unal"

#' 29 3D landmarks from the skulls of 19 carnivoran specimens before Procrustes analysis
#'
#' A matrix containing the 29 3D landmarks collected from the skulls of 21 "Vulpes vulpes"
#' specimens before carrying out a Procrustes analysis (PA).
#' Please take a look at the description in morpho/data-raw/V.R to understand
#' how this object was generated.
#'
#' @format A matrix with s = 21 rows and n = 87 columns (87/3 = 29 landmarks):
#'   \describe{
#'   \item{s}{Rows, specimens from which landmarks were collected, 21}
#'   \item{n}{Columns, 87 traits (29 landmarks in 3D) after the PA}
#' }
"V.mat.unal"

#' Object of class array output by Morpho after PA
#'
#' Object of class array output by Morpho after being aligned to the
#' mean shape of the 19 carnivoran species previously generated (object "C.PS").
#' Please take a look at the description in morpho/data-raw/V.R to understand
#' how this object was generated.
#'
#' @format Object procSym
#'   \describe{
#'   \item{...}{Check \code{\link[Morpho]{procSym}} for more details}
#' }
"V.PS.nov1"

#' 29 3D landmarks from the skulls of 19 carnivoran specimens after Procrustes analysis
#'
#' A matrix containing the 29 3D landmarks collected from the skulls of 21 "Vulpes vulpes"
#' specimens after carrying out a Procrustes analysis (PA).
#' Please take a look at the description in morpho/data-raw/V.R to understand
#' how this object was generated.
#'
#' @format A matrix with s = 21 rows and n = 87 columns (87/3 = 29 landmarks):
#'   \describe{
#'   \item{s}{Rows, specimens from which landmarks were collected, 21}
#'   \item{n}{Columns, 87 traits (29 landmarks in 3D) after the PA}
#' }
"V"

#' Vector with the population variance of Vulpes vulpes
#'
#' Vector with 87 within-species variances calculated from the object V.
#' Please take a look at the description in morpho/data-raw/var.foxes.R to understand
#' how this object was generated.
#'
#' @format A vector with i = 87 variances regarding the "Vulpes vulpes"
#' population:
#' \describe{
#'   \item{i}{Number of variances for the foxes population, 87}
#' }
"var.foxes"

#' 29 3D landmarks from the skulls of 21 Vulpes vulpes specimens
#'
#' A dataset containing the 29 3D landmarks collected from the skulls of 19 carnivoran specimens
#' This data.drame consists of a first column with the specimens labels used by
#' MCMCtree and then 87 columns with the landmarks (29 landmarks x 3 coordinates).
#' Please take a look at the description in morpho/data-raw/vulpes21x29.raw.R to understand
#' how this object was generated.
#'
#' @format A data.frame with n = 21 rows and p = 88 columns (info column + 87 coordinates):
#' \describe{
#'   \item{n}{Rows, specimens from which landmarks were collected, 21}
#'   \item{p}{Columns, information about the taxa (1st column) and 87 coordinates (29 landmarks in 3D)}
#' }
"vulpes21x29.raw"