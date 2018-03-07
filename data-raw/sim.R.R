# Simulate a correlation matrix for the examples in
# sim.morpho()
# The matrix follows the constant correlation model, i.e. with all the
# correlations equal to a specific value rho.
# We here use rho = 0.50 and simulate a correlation matrix
# for n = 100 characters (nxn)
n             <- 100
rho           <- 0.50

Id            <- diag( 1, 100 )
sim.R         <- matrix( rep( rho ), ncol = n, nrow = n )
diag( sim.R ) <- 1

# Save the R.unb
devtools::use_data( sim.R )
