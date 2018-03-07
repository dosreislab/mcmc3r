context("calc.pop.cor")

test_that(
	"There is a warning if the R.shrunk matrix generated
	is not invertible",
	{
	expect_warning(
	  morpho::calc.pop.cor( P = sim.population, delta = 10^-20 )
	)
	}
)
