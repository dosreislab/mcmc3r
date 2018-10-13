context("array2matrix, matrix2array")

test_that(
	"A. There is an error if coords != 2 or 3",
	{
	expect_error(
	  array2matrix(  X = C.arr.unal, coords = 4 )
	)
	}
)

test_that(
  "B. There is an error if coords != 2 or 3",
  {
    expect_error(
      matrix2array(  X = C.mat.unal, coords = 4 )
    )
  }
)

test_that(
  "C. There is an error if coords != dim( array )[2] ",
  {
    expect_error(
      array2matrix(  X = C.arr.unal, coords = 2 )
    )
  }
)





