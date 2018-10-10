context("array2matrix, matrix2array")

test_that(
	"There is an error if not class array",
	{
	expect_error(
	  array2matrix(  X = C.mat.unal, coords = 3 )
	)
	}
)

test_that(
  "There is an error if not class matrix",
  {
    expect_error(
      matrix2array(  X = C.arr.unal, coords = 3 )
    )
  }
)
