# Get the estimate of the correlation matrix with the shrinkage
# method. Use the function corpcor::cor.shrink for
# that purpose.

library(corpcor)

R <- corpcor::cor.shrink(vulpes21x29.matrix)

# Note: The shrinkage parameter is d = 0.6352. Note that
# the closer this parameter controls the level of shrinkage
# between the unbiased estimate of the correlation matrix (R')
# and the identity matrix (I). When d = 0, R = R'; when d = 1,
# R = I.

# Convert object "shrinkage" into "matrix"

class( R ) <- "matrix"

# Save shrunk correlation matrix estimate
devtools::use_data(R)
