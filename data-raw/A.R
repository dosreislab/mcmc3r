# Get the A matrix which, once multiplied by its transpose,
# yields to the inverse of the estimated shrinkage correlation
# matrix (R.sh). This matrix is later used to transform the data set as
# it corrects the morphological data set for the corresponding
# character correlation

U <- chol( R.sh )
A <- backsolve( U, diag( dim( U )[1] ) )

# Save A matrix
devtools::use_data( A )
