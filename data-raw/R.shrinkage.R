# Estimate the shrinkage correlation matrix using the corpcor
# The delta parameter found by this package is supposed to be
# the optimum one. This parameter controls the level of shrinkage
# between the unbiased estimate of the correlation matrix, cor( V ),
# and the identity matrix.

# Generate shrinkage correlation matrix, R.sh
R.sh <- corpcor::cor.shrink( V )

# Convert R.sh into class matrix
R.sh <- cbind( R.sh )

# Save shrinkage correlation matrix estimate
devtools::use_data( R.sh )
