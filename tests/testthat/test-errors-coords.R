context("array2matrix, matrix2array")

test_that(
	"A. There is an error if coords != 2 or 3",
	{
	expect_error(
	  array2matrix(  X = canids19x29.array, coords = 4 )
	)
	}
)

test_that(
  "B. There is an error if coords != 2 or 3",
  {
    expect_error(
      matrix2array(  X = canids19x29.matrix, coords = 4 )
    )
  }
)

test_that(
  "C. There is an error if coords != dim( array )[2] ",
  {
    expect_error(
      array2matrix(  X = canids19x29.array, coords = 2 )
    )
  }
)





