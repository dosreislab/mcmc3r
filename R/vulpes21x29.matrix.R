#' 29 3D landmarks from the skulls of 19 carnivoran specimens after Procrustes analysis
#'
#' A matrix containing the 29 3D landmarks collected from the skulls of 21 Vulpes vulpes
#' specimens after carrying out a Procrustes analysis (PA). You can follow a detailed
#' tutorial on how to perform the PA here:
#' \url{https://sabifo4.github.io/blog/Morphometrics_and_Procrustes_alignment}
#'
#' @format A matrix with s = 21 rows and n = 87 columns (87/3 = 29 landmarks):
#'   \describe{
#'   \item{s}{Rows, specimens from which landmarks were collected, 21}
#'   \item{n}{Columns, 87 coordinates (29 x, 29 y, and 29 z) after the PA}
#' }
#' @source \url{https://github.com/sabifo4/Morphometrics_practical}
"vulpes21x29.matrix"
