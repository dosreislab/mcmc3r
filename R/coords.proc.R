#' 48 3D landmarks from the skulls of 8 Triturius specimens and 1 Calotriton specimen after Procrustes analysis
#'
#' A dataset containing the 48 3D landmarks collected from the skulls of 8 Triturius specimens
#' and 1 Calotriton specimen (see Ivanovic A, Arntzen JW for more information).
#' This 3D array includes 48 landmarks collected in 3D from 9 specimens. You can follow
#' a detailed tutorial on how to perform the Procrustes analysis that generated this data set
#' here: \url{https://sabifo4.github.io/blog/Morphometrics_and_Procrustes_alignment}
#'
#' @format An array with p = 48 (landmarks), k = 3 (coordinates) and n = 9 (species):
#' \describe{
#'   \item{p}{landmark points collected from Triturius and Calotriton specimens, 48}
#'   \item{k}{coordinates, 3D}
#'   \item{n}{number of speciments, 9}
#' }
#' @source \url{https://github.com/sabifo4/Morphometrics_practical}
"coords.proc"
