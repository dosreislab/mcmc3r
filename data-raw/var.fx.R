# Calculate the var-cov matrix once the PA has been done
# with the landmark points collected from the 21 Vulpes vulpes
# specimens

VCV.fx <- cov( vulpes21x29.matrix )

# Get the diagonal of the var-cov matrix, i.e. the vector of
# population variances

var.fx <- diag( VCV.fx )

devtools::use_data( var.fx )
