# Simulate a population with a sample of s = 20 specimens
# under a normal distribution with mean = 0, variance c = 0.25,
# i.e. x ~ N(0,0.25), and R = sim.R (constant correlation model with
# rho = 0.50 and n = 100 characters, thus size is n x n).

sim.population <- morpho::sim.pop( psample = 20, n = 100, c = 0.25, R = sim.R )


# Save the R.unb
devtools::use_data( sim.population )
