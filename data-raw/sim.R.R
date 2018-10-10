# Simulate a correlation matrix for the examples in
# sim.morpho()
# The matrix follows the constant correlation model, i.e. with all the
# correlations equal to a specific value of rho.
# We here use rho = 0.50 and simulate a correlation matrix
# for n = 100 characters ( n x n )

# Get parameters
n   <- 100
rho <- 0.50

# Get true correlation matrix used for simulation
sim.R         <- matrix( rep( rho ), ncol = n, nrow = n )
diag( sim.R ) <- 1

# Save the unbiased correlation matrix
devtools::use_data( sim.R )
