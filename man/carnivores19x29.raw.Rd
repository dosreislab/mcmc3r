% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/carnivores19x29.raw.R
\docType{data}
\name{carnivores19x29.raw}
\alias{carnivores19x29.raw}
\title{29 3D landmarks from the skulls of 19 carnivoran specimens}
\format{
A data.frame with s = 19 rows and p = 89 columns ( 2 info columns + 87 traits ):
\describe{
  \item{s}{Rows, specimens from which landmarks were collected, 19}
  \item{p}{Columns, specimens information in 1st and 2nd column + 87 traits
   (29 landmarks in 3D)}
}
}
\usage{
carnivores19x29.raw
}
\description{
A dataset containing the 29 3D landmarks collected from the skulls of 19 carnivoran specimens
This data.frame consists of a first column with the specimen labels used by MCMCtree,
a second column with the voucher names of the specimens collected, and then 87 columns
with the coordinates for each trait (29 landmarks x 3D coordinates).
Please take a look at the description in morpho/data-raw/carnivores19x29.raw.R to understand
how this object was generated.
}
\keyword{datasets}
