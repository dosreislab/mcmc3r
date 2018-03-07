# If this has not previously done, the R.unb is neeeded
# (see R.unb.R for more details).
# In summary, the following commented line should be run
# to generate this object. However, as it is already generated
# we will just use it.
#R.unb <- cor( canids19x29.matrix )

# Get the shrunk estimate of the correlation matrix using a
# delta parameter close to 0.
# This value is the one that has proved to be good with large
# and small variances in the simulations
delta             <- 0.01
Id                <- diag( 1, dim( canids19x29.matrix )[2] )
R.shrunk          <- delta*Id + (1 - delta)*R.unb
class( R.shrunk ) <- "matrix"

# Remember that the delta parameter controls the level of shrinkage
# between the unbiased estimate of the correlation matrix (R.unb)
# and the identity matrix (Id). Therefore, the smaller the delta
# value, the closer the estimate of the shrunk correlation matrix is
# to the unbiased estimate of correlation matrix.

# Save shrunk correlation matrix estimate
devtools::use_data( R.shrunk )
