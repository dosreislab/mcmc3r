#' 29 3D landmarks from the skulls of 21 Vulpes vulpes specimens
#'
#' A dataset containing the 29 3D landmarks collected from the skulls of 21 carnivoran specimens
#' This data.drame consists of a first column with the taxa and then 87 columns
#' with the landmarks (29 landmarks x 3 coordinates). You can follow
#' a detailed tutorial on how to perform a Procrustes analysis (PA) here:
#' \url{https://sabifo4.github.io/blog/Morphometrics_and_Procrustes_alignment}
#'
#' @format A data.frame with n = 21 rows and p = 87 columns (info column + 87 coordinates):
#' \describe{
#'   \item{n}{Rows, specimens from which landmarks were collected, 21}
#'   \item{p}{Columns, information about the taxa (1st column) and 87 coordinates (29 x, 29 y, and 29 z)}
#' }
#' @source \url{http://datadryad.org/resource/doi:10.5061/dryad.nr210}
"vulpes21x29.raw"
