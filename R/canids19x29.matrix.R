#' 29 3D landmarks from the skulls of 19 carnivoran specimens after Procrustes analysis
#'
#' A matrix containing the 29 3D landmarks collected from the skulls of 19 carnivoran
#' specimens after carrying out a Procrustes analysis (PA). You can follow a detailed
#' tutorial on how to perform the PA here:
#' \url{https://sabifo4.github.io/blog/Morphometrics_and_Procrustes_alignment}
#'
#' @format A matrix with n = 19 rows and p = 87 columns (87/3 = 29 landmarks):
#'   \describe{
#'   \item{n}{Rows, specimens from which landmarks were collected, 19}
#'   \item{p}{Columns, 87 coordinates (29 x, 29 y, and 29 z) after the PA}
#' }
#' @source \url{https://github.com/sabifo4/Morphometrics_practical}
"canids19x29.matrix"
