# Calculate the var-cov matrix once the PA has been done
# with the landmark points collected from the 21 foxes
# specimens. Afterwards, get the diagonal in order to obtain
# the vector of within-species variances

# Get the vector of variances for the foxes
var.foxes <- diag( cov( V ) )

# Save it in an rda object
devtools::use_data( var.foxes )
