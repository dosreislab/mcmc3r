context("write.morpho")

test_that(
	"There is an error when tree is not of class phylo",
	{
	  tree.wrong <- c( "((((A:0.1,B:0.1):0.2,(F:0.1,C:0.2):0.1):0.5,H:0.1):0.2,(D:0.7,(G:0.2,E:0.5):0.2):0.3);" )

	expect_error(
	  sim.morpho( tree = tree.wrong, n = 100, c = 0.25 )
	)
	}
)


test_that(
  "There is an error when c is negative",
  {
    expect_error(
      sim.morpho( tree = tree.wrong, n = 100, c = -0.25 )
    )
  }
)


test_that(
  "There is an error when c is not length 1 or length n",
  {
    c <- c(0.25, 0.25, 0.30, 0.40)
    expect_error(
      sim.morpho( tree = sim.tree, n = 100, c = c )
    )
  }
)


test_that(
  "There is an error when R is not square",
  {
    R <- matrix( rep( 0.25 ),
                 ncol = dim( C )[2],
                 nrow = 20 )

    diag( R ) <- 1

    expect_error(
      sim.morpho( tree = sim.tree, n = 100, c = 0.25, R = R )
    )
  }
)


