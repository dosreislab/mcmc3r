context("write.morpho")

test_that(
	"There is an error when M is not a matrix",
	{
	expect_error(
	  write.morpho(  M = C.arr.unal, filename = "seqfile.aln" )
	)
	}
)


test_that(
  "There is an error when c is negative",
  {
    expect_error(
      write.morpho(  M = C, c = -0.25,
                     filename = "seqfile.aln" )
    )
  }
)


test_that(
  "There is an error when c is not length 1 or length n",
  {
    c <- c(0.25, 0.25, 0.30, 0.40)
    expect_error(
      write.morpho(  M = C, c = c,
                     filename = "seqfile.aln" )
    )
  }
)


test_that(
  "There is an error when R is not positive definite",
  {
    R <- cor( C )
    expect_error(
      write.morpho(  M = C, c = 0.25,
                     R = R, method = "chol",
                     filename = "seqfile.aln" )
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
      write.morpho(  M = C, c = 0.25,
                     R = R, method = "chol",
                     filename = "seqfile.aln" )
    )
  }
)


test_that(
  "There is an error if method provided but missing R",
  {
    expect_error(
      write.morpho(  M = C, c = 0.25,
                     method = "chol", filename = "seqfile.aln" )
    )
  }
)


test_that(
  "There is an error if method provided is not chol or eigen",
  {
    expect_error(
      write.morpho(  M = C, c = 0.25, R = R.sh,
                     method = "g", filename = "seqfile.aln" )
    )
  }
)


test_that(
  "There is an error when 'method' is not used if R is passed",
  {
    expect_error(
      write.morpho(  M = C, c = 0.25,
                     R = R.sh, filename = "seqfile.aln" )
    )
  }
)


test_that(
  "There is an error when there are more or less names than
  the rows in matrix M",
  {
    names.wrong <- list( sp1  = "Ael_sp.", sp2  = "Can_dir", sp3  = "Epi_hay", sp4  = "Hes_sp.",
                   sp5  = "Mes_cor", sp6  = "Tom_sp.", sp7  = "Enh_pah", sp8  = "Cuo_alp",
                   sp9  = "Spe_ven", sp10 = "Can_lup", sp11 = "Cer_tho", sp12 = "Oto_meg",
                   sp13 = "Vul_vul", sp14 = "Urs_ame", sp15 = "Ail_ful", sp16 = "Nan_bio",
                   sp17 = "Par_her", sp18 = "Tha_won"
                   )
    expect_error(
      write.morpho( M = C, c = 0.25, R = R.sh,
                    filename = "seqfile.aln", names = names.wrong )
    )
  }
)


test_that(
  "There is an error when there are more or less ages than
  the rows in matrix M",
  {
    ages.wrong <- list( sp1  = 15.97, sp2  =  1.80, sp3  = 13.6,  sp4  = 39.74,
                  sp5  = 30.80, sp6  = 15.97, sp7  = 30.80, sp8  =  0,
                  sp9  =  0,    sp10 =  0,    sp11 =  0,    sp12 =  0,
                  sp13 =  0,    sp14 =  0,    sp15 =  0,    sp16 =  0,
                  sp17 =  0,    sp18 =  0
                  )
    expect_error(
      write.morpho( M = C, c = 0.25, R = R.sh,
                    filename = "seqfile.aln", names = names.wrong,
                    ages = ages.wrong )
    )
  }
)
