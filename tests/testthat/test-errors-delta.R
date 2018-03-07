context("calc.pop.cor")

test_that(
	"There is an error if delta < 0",
	{
	expect_error(
	  morpho::calc.pop.cor( P = sim.population, delta = -0.3 )
	)
	}
)


test_that(
  "There is an error if delta > 1",
  {
    expect_error(
      morpho::calc.pop.cor( P = sim.population, delta = 1.1 )
    )
  }
)
