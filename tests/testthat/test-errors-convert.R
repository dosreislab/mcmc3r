context("array2matrix, matrix2array")

test_that(
	"There is an error if not class array",
	{
	expect_error(
	  array2matrix(  X = canids19x29.matrix, coords = 3 )
	)
	}
)

test_that(
  "There is an error if not class matrix",
  {
    expect_error(
      matrix2array(  X = canids19x29.array, coords = 3 )
    )
  }
)
