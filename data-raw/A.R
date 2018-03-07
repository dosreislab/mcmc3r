# Get the A matrix which, once multiplied by its transpose, 
# yields to the inverse of the estimate of the shrunk correlation 
# matrix. This matrix is later used to transform the data set as 
# it corrects the morphological data set for the corresponding
# correlation

U <- chol( R.shrunk )
A <- backsolve( U, diag( dim( U )[1] ) )

# Save A matrix 
devtools::use_data( A )
