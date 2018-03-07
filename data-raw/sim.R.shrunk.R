# Estimate the shrunk correlation matrix of
# a population sample matrix.
# For this estimate, we are using the previously generated object
# "sim.population" as the population sample (20 specimens).
# We are using the function calc.pop.cor to generate the estimate
# and the default value of delta = 0.01.
# Note that this function returns a list with the estimate of the
# shrunk correlation matrix and the variance of the population.
# Therefore, we just return now the former.

vars <- calc.pop.cor( P = sim.population, delta = 0.01 )
sim.R.shrunk <- vars$R.shrunk

# Save the R.unb
devtools::use_data( sim.R.shrunk )
