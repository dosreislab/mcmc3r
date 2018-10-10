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
