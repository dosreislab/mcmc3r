context("sim.morpho, sim.pop")

test_that(
	"A. There is an error if n is not an integer in sim.morpho",
	{
	expect_error(
	  sim.morpho( tree = sim.tree, n = 100.56, c = 0.25 )
	)
	}
)


test_that(
  "B. There is an error if n is not an integer in sim.pop",
  {
    expect_error(
      sim.pop( psample = 20, n = 100.56, c = 0.25 )
    )
  }
)


test_that(
  "C. There is an error if n is not positive",
  {
    expect_error(
      sim.pop( psample = 20, n = -45, c = 0.25 )
    )
  }
)


test_that(
  "D. There is an error if n is not positive",
  {
    expect_error(
      sim.pop( psample = 20, n = 0, c = 0.25 )
    )
  }
)
