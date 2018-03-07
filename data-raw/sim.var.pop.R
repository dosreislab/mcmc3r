# Estimate the population variance of a population sample.
# For this estimate, we are using the previously generated object
# "sim.population" as the population sample (20 specimens).
# We are using the function calc.pop.cor, which generates the estimate
# of the shrunk correlation matrix but additionally returns the
# variance of the population matrix. Note that these two variables
# are to be used within the write.morpho function when generating
# the alignment file to be input in MCMCTree.

vars <- calc.pop.cor( P = sim.population, delta = 0.01 )
sim.var.pop <- vars$var

# Save the R.unb
devtools::use_data( sim.var.pop )
