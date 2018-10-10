#' Generate a file with a morphological alignment for MCMCtree
#'
#' @description
#' Generate an alignment file with quantitative characters in MCMCtree format.
#' The option "seqfile" in the control file used by MCMCtree
#' should read the path to the file output by this function.
#'
#' @param M Matrix, \code{s} rows (specimens) and \code{n} morphological
#' continuous characters (see details).
#'
#' @param filename Character, name for the output file.
#'
#' @param c Numeric, vector of variances within the species of a population (see details).
#' If not provided, \code{c = 0} (no population noise).
#'
#' @param R Matrix, correlation matrix. Requires \code{method}
#' (see details). If not provided, \code{R = I} (no character correlation).
#'
#' @param method (Optional) character, either \code{"eigen"} or
#' \code{"col"}, method used to decompose the inverse of the shrinkage
#' correlation matrix. Requires \code{R} (see details).
#'
#' @param A (Optional) matrix, decomposed matrix. Requires \code{R}
#' but not \code{method} (see details).
#'
#' @param names (Optional) list or character, species name included in the morphological
#' alignment (see examples B and C).
#'
#' @param ages (Optional) list or character, ages of the species included in the
#' morphological alignment (see example C).
#'
#' @details
#'
#' The matrix \code{M} has \code{s} rows, one for each specimen, and \code{n}
#' columns regarding the characters. If the data set contains landmarks,
#' they can be given in 2D or 3D.
#' For instance, if the landmarks are 3D, the first 3 columns will be the
#' coordinates x, y, and z for the first landmark; the next 3
#' columns for the second landmark; and so on.
#' \tabular{cccccccc}{
#'  specimens \tab lmk1.x \tab lmk1.y  \tab lmk1.z  \tab lmk2.x \tab lmk2.y \tab lmk2.z  \tab ... \cr
#'  Sp_1      \tab 0.143  \tab -0.028  \tab -0.044  \tab 0.129  \tab 0.028  \tab -0.043  \tab ... \cr
#'  Sp_2      \tab 0.128  \tab -0.024  \tab -0.028  \tab 0.124  \tab 0.027  \tab -0.025  \tab ... \cr
#'  ...       \tab ...    \tab ...     \tab ...     \tab ...    \tab ...    \tab ...     \tab ...
#' }
#' See descriptions in \code{data-raw/carnivores19x29.raw.R},
#' \code{data-raw/vulpes21x29.raw.R}, \code{data-raw/C.R}, \code{data-raw/V.R},
#' to know how to obtain the morphological alignment used in \code{write.morpho}.
#' Note that the explanation starts with the processing of raw data.
#'
#' If the data set contains a set of \code{n} morphological continuous characters,
#' e.g. from a simulated data set, the file should look like
#' \tabular{cccccc}{
#'  specimens \tab char.1  \tab char.2   \tab char.3   \tab char.4  \tab ... \cr
#'  Sp_1      \tab 0.143   \tab -0.028   \tab -0.044   \tab 0.129   \tab ... \cr
#'  Sp_2      \tab 0.128   \tab -0.024   \tab -0.028   \tab 0.124   \tab ... \cr
#'  ...       \tab ...     \tab ...      \tab ...      \tab ...     \tab ...
#' }
#'
#' Note that if a list with the specimens names is not passed to the parameter \code{names},
#' the name for each species will be "Species_1", "Species_2", and so on.
#'
#' The object \code{c} can be of length 1, if all characters have the same variance, or
#' a vector of length \code{n} with the variance of each of the characters. For the latter,
#' you can take a look at object \code{var.foxes}, which has been generated following
#' the steps explained in \code{data-raw/var.foxes.R}.
#'
#' The object \code{R} has to be a symmetric and positive definite object of class matrix,
#' i.e. \code{class( R ) = "matrix"}. See \code{R.sh} for an example of its format.
#' You can also read the description in \code{data-raw/R.shrinkage.R} for the details
#' about how to generate this matrix.
#'
#' The logarithm of the determinant of the correlation matrix is going to be printed
#' in the output file to later be used by MCMCtree during the likelihood calculation.
#'
#' If a correlation matrix \code{R} is provided, \code{write.morpho} can use either
#' the \code{method = "chol"} or \code{method = "eigen"} to get a matrix \code{A}
#' such that \eqn{\mathrm{R^{-1}}=\mathrm{A^{T}}\mathrm{A}}{R^{-1} = t(A) * A}.
#' This matrix \eqn{\mathrm{A^{T}}}{t(A)} is later used to transform
#' the morphological data to account for the correlation in this data set,
#' so that the transformed characters in \code{Z},
#' \eqn{\mathrm{Z}=\mathrm{M}\mathrm{A^{T}}}{Z = M * t(A)}, are independent.
#' Alternatively, this matrix \code{A} can also be provided by the user. You can
#' read the description to generate this matrix in \code{data-raw/A.R} to
#' understand how to generate this matrix, which is also available as object \code{A}.
#' If you decide to use a matrix \code{A}, it will be used to transform the data
#' and no decomposition will be performed, thus saving computational time when
#' large matrices are to be used.
#'
#' @seealso
#' \code{\link{matrix2array}}, \code{\link{array2matrix}}, \code{\link{sim.morpho}}
#'
#' @author Sandra Alvarez-Carretero and Mario dos Reis
#'
#' @examples
#' # A.1) Providing the morphological alignment (M = C) and
#' #      the name for the output file. This does not account for
#' #      correlation nor population noise
#'
#'        write.morpho(  M = C, filename = "seqfile.aln" )
#'
#' # A.2) Providing the morphological alignment (M = C), the population
#' #      noise (c = 0.25), and the name for the output file. Note that
#' #      c = 0.25 means that the population noise for all the characters
#' #      is c = 0.25, i.e. it will be considered as if
#' #      length( c ) = p characters, being all of them 0.25.
#'
#'        write.morpho(  M = C, c = 0.25,
#'                       filename = "seqfile.aln" )
#'
#' # A.3) Providing the morphological alignment (M = C), the population
#' #      noise (c = 0.25), the (estimate of the) correlation matrix (R),
#' #      the method to decompose R ("chol" in this example),
#' #      and the name for the output file. Note that the R matrix needs
#' #      to be invertible, otherwise the data will not be able to be
#' #      transformed accounting for correlation.
#'
#'        write.morpho(  M = C, c = 0.25, R = R.sh,
#'                       method = "chol", filename = "seqfile.aln" )
#'
#' # A.4) Providing the morphological alignment (M = C), a vector with
#' #      the population noise for each character (c = var.foxes),
#' #      the (shrinkage estimate of the) correlation matrix (R = R.sh),
#' #      the method to decompose R ("chol" in this example),
#' #      and the name for the output file. Note that the matrix passed
#' #      to argument R needs to be invertible, otherwise the data will
#' #      not be able to be transformed accounting for correlation.
#'
#'        write.morpho(  M = C, c = var.foxes, R = R.sh,
#'                       method = "chol", filename = "seqfile.aln" )
#'
#' # A.5) Providing the morphological alignment (C), a vector with
#' #      the population noise for each character (c = var.foxes),
#' #      the (shrinkage estimate of the) correlation matrix (R = R.sh),
#' #      the A matrix to transform the data, and the name for the
#' #      output file. Note that as the A matrix is provided, the matrix
#' #      passed to R will not be decomposed, hence the argument "method"
#' #      is not needed.
#'
#'        write.morpho(  M = C, c = var.foxes, R = R.sh,
#'                       A = A, filename = "seqfile.aln" )
#'
#' # B) Scenario A.5 but providing a list with the
#' #    names of the species
#'
#'      names <- list( sp1  = "Ael_sp.", sp2  = "Can_dir", sp3  = "Epi_hay", sp4  = "Hes_sp.",
#'                     sp5  = "Par_jos", sp6  = "Tom_sp.", sp7  = "Enh_pah", sp8  = "Cuo_alp",
#'                     sp9  = "Spe_ven", sp10 = "Can_lup", sp11 = "Cer_tho", sp12 = "Oto_meg",
#'                     sp13 = "Urs_ame", sp14 = "Ail_ful", sp15 = "Nan_bin", sp16 = "Par_her",
#'                     sp17 = "Hia_won", sp18 = "Smi_fat", sp19 = "Vul_vul"
#'                   )
#'
#'      write.morpho( M = C, c = var.foxes, R = R.sh,
#'                    A = A, filename = "seqfile.aln", names = names )
#'
#' # C) Scenario A.5 but providing a vector of type character with the names of
#' #    the specimens and a list with their corresponding ages. Please
#' #    keep the same order in both lists, so the first specimen in the
#' #    list name corresponds to the first age in the age list, and so on.
#'
#'      names <- c( "Ael_sp.", "Can_dir", "Epi_hay", "Hes_sp.",
#'                  "Par_jos", "Tom_sp.", "Enh_pah", "Cuo_alp",
#'                  "Spe_ven", "Can_lup", "Cer_tho", "Oto_meg",
#'                  "Urs_ame", "Ail_ful", "Nan_bin", "Par_her",
#'                  "Hia_won", "Smi_fat", "Vul_vul"
#'                   )
#'
#'      ages <- list( sp1  = 13.135, sp2  =  0.0285, sp3  = 11.95,  sp4  = 35.55,
#'                    sp5  = 25.615, sp6  = 14.785,  sp7  = 28.55,  sp8  =  0,
#'                    sp9  =  0,     sp10 =  0,      sp11 =  0,     sp12 =  0,
#'                    sp13 =  0,     sp14 =  0,      sp15 =  0,     sp16 =  0,
#'                    sp17 =  6.65,  sp18 =  0.0285, sp19 =  0
#'                   )
#'
#'      write.morpho( M = C, c = var.foxes, R = R.sh,
#'                    A = A, filename = "seqfile.aln",
#'                    names = names, ages = ages )
#'
#' @export

write.morpho <- function( M, filename, c = 0, R = diag( 1, dim( M )[2] ),
                          method, A = NULL, names, ages ) {

  # Check initial arguments
  .checkInArgs( X = M, filename = filename, c = c, R = R,
                method = method, A = A )

  # Create variables
  s     <- dim( M )[1]
  chars <- dim( M )[2]

  # Check optional arguments
  .checkOptArg( names, ages, M )

  # If names are provided...
  if ( ! missing( names ) ){
    vars   <- .availNames( names, ages )
    spaces <- vars$spaces
    names  <- vars$names
  }

  # If names are not provided...
  else{
    vars   <- .notNames( ages, s )
    spaces <- vars$spaces
    names  <- vars$names
  }

  # Get names with data in the array, add spaces, and put them as rownames
  names         <- paste( names, spaces, sep = "    " )
  rownames( M ) <- names

  # Scale data accounting for population noise
  # any( c != 0 ) is T if there are not 0s in c
  if ( any( c != 0 ) ){
    M        <- .scalePopVar( X = M, c = c, s = s, n = chars )
    scalevar <- 1
  }
  else{
    scalevar <- 0
  }

  # If the A matrix is not provided and R != I, calculate it and
  # then transform M
  #if ( missing( A ) & all( R == diag( 1, dim( M )[2] ) ) == FALSE ){
  if ( is.null( A ) & all( R == diag( 1, dim( M )[2] ) ) == FALSE ){
    # Match method with the argument provided
    method <- .checkMethod( method = method )
    # Calculate transformed matrix
    Z <- .CalcZ( X = M, R = R, method = method )
    # Get logarithm determinant of R
    lnd <- determinant( R )$modulus
  }
  # Otherwise, just use the A matrix provided to transform M
  #else if ( ! missing( A ) & all( R == diag( 1, dim( M )[2] ) ) == FALSE ){
  else if ( ! is.null( A )  & all( R == diag( 1, dim( M )[2] ) ) == FALSE ){
    Z <- M %*% A
    # Get logarithm determinant of R
    lnd <- determinant( R )$modulus
  }
  # Last, if there is no correlation at all, then
  else if ( all( R == diag( 1, dim( M )[2] ) ) == TRUE ){
    Z     <- M
    lnd   <- 0
  }

  # Generate output file for MCMCtree
  .outFile( X = Z , names = names, chars = chars,
            scalevar = scalevar, lndetR = lnd, filename = filename )

}

## #############################
##  SUBFUNCTIONS write.morpho ##
## #############################

# Check inp parameters
.checkInArgs <- function( X, filename, c, R, method, A ){

  # Check a name for the output file has been given
  if ( missing( filename ) ){
    stop( "\nPlease use the parameter \"filename\" to provide a name for the output file\n" )
  }

  # Check object X (morph.data) is provided and class is a matrix
  if ( missing( X ) ){
    stop( "\nPlease use an object of class \"matrix\" and dimensions \"s x n\" to convert
          into MCMCtree format\n" )
  }
  if( class( X ) != "matrix" ){
    stop( "\nPlease use an object of class \"matrix\" and dimensions \"s x n\" to convert
          into MCMCtree format\n" )
  }

  # Check population noise
  if( class( c ) != "numeric" ){
    stop( "\nThe parameter \"c\" needs a numeric value > 0 or a numeric vector
          of length \"n\" with values > 0\n" )
  }
  if( length( c ) != dim( X )[2] & length( c ) != 1 ){
    stop( "\nProvide a vector with variances of length equal to the number of
          morphological characters, \"length( c ) = n\", or to a unique numeric value,
          \"length( c ) = 1\", which will assume that all variances equal to this value\n" )
  }
  # Check any value in the vector of vars is <0
  if( length( c ) == dim( X )[2] & any( c < 0 ) ){
    stop( "\nThe parameter \"c\" needs all numeric values > 0\n" )
  }
  # Check the popvar value is >0
  if( length( c ) == 1 & any( c < 0 ) ){
    stop( "\nThe parameter \"c\" needs a numeric value\n" )
  }

  # Check that method goes always with R provided by the user
  if ( ! missing( method ) & all( R == diag( 1, dim( X )[2] ) ) == T ){
    stop( "\nPlease provide a correlation matrix so that the
             method you have selected can be used\n" )
  }
  # Check correlation matrix is size n x n and that method
  # has been provided
  if ( dim( R )[1] != dim( X )[2] & dim( R )[2] != dim( X )[2] ){
    stop( "\nPlease use a correlation matrix of size 'n x n', where
             'n' is the amount of characters present in your
             'M' matrix\n" )
  }
  if ( all( R == diag( 1, dim( X )[2] ) ) == F ){
    #if ( missing( A ) & missing( method ) ){
    if ( is.null( A ) & missing( method ) ){
      stop( "\nPlease select a method to decompose the shrinkage correlation matrix,
          either method = \"chol\" or method = \"eigen\" \n" )
    }
    .checkCorrMat( R = R, n = dim( X )[2] )
  }


  # If A provided, check it is class "matrix" and also R is given
  #if ( ! missing( A ) & all( R == diag( 1, dim( X )[2] ) ) == TRUE ){
  if ( ! is.null( A ) & all( R == diag( 1, dim( X )[2] ) ) == TRUE ){
    stop( "\nYou do not want to account for trait correlation but you
          have passed matrix \"A\" to the function. Do not use \"A\" if you
          do not want to account for trait correlation. Otherwise, please
          upload your estimated correlation matrix together with \"A\"\n")
  }
  #if ( ! missing( A ) & class( A ) != "matrix" ){
  if ( ! is.null( A ) & class( A ) != "matrix" ){
    stop( "\nObject \"A\" needs to be of class \"matrix\"\n" )
  }
    # [ DISABLED BY NOW. IT MIGHT BE TOO TIME CONSUMING WITH LARGE MATRICES ]
    # [ IF MANY PEOPLE ARE USING THE WRONG A, WE WILL ENABLE THIS CHECK AGAIN ]
    # Check that the A matrix provided is correct given their R
    #A2 <- .CalcCholesky( R = R, n = dim( X )[2] )
    #if ( A != A2 ){
    #  stop( "\nYour \"A\" matrix does not seem to be the correct
    #         one according the \"R\" matrix you have provied.
    #         If you are unsure of properly calculating A, then
    #         use write.morpho without this argument so it is
    #         internally calculated\n" )
    #}

}


# Check if optional parameters (names and ages) added and, if so,
# if class is ok
.checkOptArg <- function( names, ages, M ){

  # Check names is class "list" OR "character" and length(names) = dim(M)[1]
  if( ! missing( names ) ){
    if ( class( names ) != "list" & class( names ) != "character" ){
      stop( "\nYou need to provide an object of class list or a character vector with the species included in the alignment file.
            E.g.1 species <- list( sp1 = \"sp1\", sp2 = \"sp2\" )
            E.g.2 species <- c(\"sp1\", \"sp2\")\n" )
    }
    if( length( names ) != dim( M )[1] ){
      stop( "\nPlease, provide a list or a vector with the same amount
               of names than the rows in matrix 'M'\n")
    }
  }

  # Check ages is class "list" and length(ages) = dim(M)[1]
  if( ! missing( ages ) ){
    if ( class( ages ) != "list" & class( ages ) != "numeric" ){
      stop( "\nYou need to provide an object of class list or a numeric vector with the ages of the species included in the alignment file.
            E.g.1 ages <- list( ag1 = 30, ag2 = 15
            E.g.2 ages <- c( 30, 15 )\n" )
    }
    if( length( ages ) != dim( M )[1] ){
      stop( "\nPlease, provide a list or a vector with the same amount
            of ages than the rows in matrix 'M'\n")
    }
  }

}


# Get ages if provided
.hasAges <- function( ages, names ){

  if ( class( ages ) == "list" ){
    ages  <- unlist( ages )
  }
  names <- paste( names, - ages + max( ages ) + 0.01, sep = "^" ) # 0.01 = ct for MCMCtree
  # Return names with ages updated
  return( names )

}


# Get num spaces and add it to str
.hasSpaces <- function( names ){

  num.spaces <- max( nchar( names ) ) - nchar( names )
  str        <- sprintf( paste( "% ", num.spaces, "s", sep = "" ), c( "" ) )

  # Return str
  return( str )
}


# Get spaces and names vectors when names are not provided.
# This wraps around .hasAges and .hasSpaces
.notNames <- function( ages, s ){

  # Name species as "Species_1", "Species_2", etc.
  names <- paste( "Species_", seq( 1:s ), sep = "" )

  # If ages are not provided ...
  if ( missing(ages) ){
    # Get spaces to justify text
    spaces <- .hasSpaces( names = names )
  }

  # If ages are provided ...
  else{
    # Put together names and transformed ages and get spaces
    names <- .hasAges( ages = ages, names = names )
    # Get spaces to justify text
    spaces <- .hasSpaces( names = names )
  }

  # Return names and spaces
  return( list( names = names, spaces = spaces ) )

}


# Get spaces when names are provided.
# This wraps around .hasAges and .hasSpaces
.availNames <- function( names, ages ){

  # Get species names
  if ( class( names ) == "list" ){
    names   <- unlist( names )
  }

  # If ages are not provided ...
  if ( missing( ages ) ){
    # Get spaces to justify text
    spaces <- .hasSpaces( names = names )
  }

  # Otherwise ...
  else{
    # Put together names and transformed ages and get spaces
    names <- .hasAges( ages = ages, names = names )
    # Get spaces to justify text
    spaces <- .hasSpaces( names = names )
  }

  # Return spaces
  return( list( spaces = spaces, names = names ) )

}


# Scale matrix "X" with population noise
# Note that c has to be a vector with all the variances
# for the species of the population,
# so if length(c) == 1 it is first corrected to be of
# length "n" with the same "c" pop variance
.scalePopVar <- function( X, c, s, n ){

  if ( length( c ) == 1 ){
    c <- rep( c, n )
  }
  M.s <-  X %*% diag( 1 / sqrt( c ) )

  # Return scaled matrix
  return( M.s )

}


# Check method has been given
.checkMethod <- function( method ){

  if ( missing( method ) ){
    stop( "\nPlease select a method to decompose the shrinkage correlation matrix,
          either method = \"chol\" or method = \"eigen\" \n" )
  }

  # Match arguments and return the method
  methods <- c( "chol", "eigen" )
  method <- match.arg( arg = method, choices = methods )
  return( method )

}

# Decompose the correlation matrix using the Cholesky decomposition
.CalcCholesky <- function( R ){

  # Cholesky decomposition:
  # R = L %*% U = L %*% t( L )    = t( U ) %*% U
  # R^-1 = t( L^-1 ) %*% L^-1     = U^-1 %*% t( U^-1 )
  # R^-1 = t( A ) %*% A           = A %*% t( A )

  U <- chol( R )
  A <- backsolve( U, diag( dim( U )[1] ) )

  # Return inverse of upper triangular matrix ( U^-1 = A )
  return( A )

}


# Decompose the correlation matrix using the eigen decomposition
.CalcEigen <- function( R ){

  # "eigen" returns a list with the eigenvectors
  # and the eigenvalues
  #
  # R^-1 = t( A ) %*% A
  #
  # t( A ) = V %*% D
  # t( A ) = eigen( R^-1 )$vectors %*% diag( sqrt( eigen( R^-1 )$values ) )

  #Rinv <- solve( R )
  Rinv <- chol2inv( chol( R ) )
  tA   <- eigen( Rinv )$vectors %*% diag( sqrt( eigen( Rinv )$values ) )

  # Return t(A)
  return( tA )

}


# Calculate transformed matrix when A is not provided
# and R needs to be decomposed
.CalcZ <- function( X, R, method ){

  # Match argument
  if ( method == "chol" ){
    A <- .CalcCholesky( R = R )
  }

  else if ( method == "eigen" ){
    A <- .CalcEigen( R = R )
  }

  # Transform data
  Z <- X %*% A

  # Return Z
  return( Z )

}

# [[ NOT USED, THE CALCULATION HAS BEEN CALLED IN MAIN ]]
# Calculate transformed matrix when A is provided
.CalcZwithA <- function( X, A ){

  Z <- X %*% A

  # Return Z
  return( Z )

}


# Write output file with default parameters.
# It wraps around the function write()
.outFile <- function( X, names, chars, scalevar, lndetR, filename ){

  string <- paste( length( names ), chars, "M", scalevar, lndetR, sep = "  " )
  cat( "\n", file = filename )
  write( paste( "  ", string, sep = "" ), file = filename, append = T )
  cat( "\n", file = filename, append = T )
  write.table( X, file = filename, append = T,
               sep = " ", row.names = T, col.names = F, quote = FALSE )
  cat( "\n", file = filename, append = T )

}



# ------------------------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------------------------- #


#' Convert a matrix into an array
#'
#' @description
#' Convert a matrix with landmark points into an object of class array.
#'
#' @param X Matrix  of size \code{s x n}, \code{n} landmark points for \code{s} specimens
#' (see details).
#'
#' @param coords Numeric, 2 or 3 for 2D or 3D landmarks, respectively.
#'
#' @details
#'
#' The matrix has format \code{s x n}, with \code{s} rows, one per specimen, and \code{n}
#' columns, one for each coordinate of the landmarks.
#' Each landmark can be given in 2D or 3D. For instance,
#' if the landmarks are 3D, the first 3 columns will be the
#' coordinates x, y, and z for the first landmark, the next 3
#' columns for the second landmark; and so on:
#' \tabular{cccccccc}{
#'  specimens \tab lmk1.x \tab lmk1.y  \tab lmk1.z  \tab lmk2.x \tab lmk2.y \tab lmk2.z  \tab ... \cr
#'  Sp_1      \tab 0.143  \tab -0.028  \tab -0.044  \tab 0.129  \tab 0.028  \tab -0.043  \tab ... \cr
#'  Sp_2      \tab 0.128  \tab -0.024  \tab -0.028  \tab 0.124  \tab 0.027  \tab -0.025  \tab ... \cr
#'  ...       \tab ...    \tab ...     \tab ...     \tab ...    \tab ...    \tab ...     \tab ...
#' }
#'
#' See object \code{C} for an example of its format and \code{data-raw/C.R} to
#' see how this object is generated.
#'
#' @return
#'
#' An object of class array with format \code{k x q x s}, where \code{k} is the number of
#' landmarks, \code{q} the number of coordinates, and \code{s} the number of specimens.
#'
#' Note that if the matrix provided does not have rownames, the specimens in the
#' returned array (names for the 's' matrices accessed through the array, i.e.
#' \code{dimnames( array )[ 3 ]} will be labelled as '1', '2', and so on.
#' See object \code{C.arr.unal} for an example of the format of the object that is
#' returned and \code{data-raw/C.R} for the description of how to
#' obtain this object.
#'
#' @seealso
#' \code{\link{array2matrix}}, \code{\link{write.morpho}}
#'
#' @author Sandra Alvarez-Carretero and Mario dos Reis
#'
#' @export

matrix2array <- function( X, coords = c( 2, 3 ) ){

  # Check X and coords are provided and their classes
  .checkMat( X )
  coords <- as.numeric( match.arg( arg = as.character(coords), choices = c( "2", "3" ) ) )
  .checkCoords( coords )

  # Create variables
  chars       <- dim( X )[2]
  s           <- dim( X )[1]
  list.coords <- c("x", "y", "z")

  # Call function to create array
  ma <- .genArr( X = X, coords = coords, chars = chars,
                 s = s, list.coords = list.coords )

  # Return 3D array
  return( ma )

}

## ################
##  SUBFUNCTIONS ##
## ################

# Check if X is provided and is class "matrix"
.checkMat <- function( X ){

  # Check object with lmks is provided and it is of class "matrix"
  if ( missing( X ) ){
    stop( "\nPlease provide an object of class \"matrix\" with \"s x n\" dimensions, 's' specimens and 'n' characters\n" )
  }

  if ( class( X ) != "matrix" ){
    stop( "\nPlease use an object of class \"matrix\"\n" )
  }

}


# Check coords
.checkCoords <- function( coords ){

  # Check coords are provided and are numeric
  if ( missing( coords ) ){
    stop( "\nThe parameter \"coords\" needs a numeric value. Please use \"coords = 2\" if the coordinates
          are 2D or \"coords = 3\", if 3D\n" )
  }
  if ( as.integer( coords ) != coords ){
    stop( "\nThe parameter \"coords\" needs a numeric value. Please use \"coords = 2\" if the coordinates
          are 2D or \"coords = 3\", if 3D\n" )
  }

}

# Wraper around .coords3D.mat2arr and .coords2D.mat2arr to create the
# array with landmarks
.genArr <- function( X, coords, chars, s, list.coords ){

  # Create an empty array to store the coordinates in the format
  # k x q x s (chars x coords x specimens)
  ma             <- array( dim = c( chars / coords, coords, s ) )
  dimnames( ma ) <- list( paste( "lmk", seq( 1:( chars / coords ) ), sep = "" ),
                          list.coords[1:coords],
                          rownames( X ) )

  # Select x, y, z positions
  if ( coords == 3 ){
    ma <- .coords3D.mat2arr( X = X , ma = ma, chars = chars, s = s )
  }

  else if ( coords == 2 ){
    ma <- .coords2D.mat2arr( X = X , ma = ma, chars = chars, s = s )
  }

  # Return array with landmarks
  return( ma )

}

# Convert 3D landmarks matrix into matrix
.coords3D.mat2arr <- function( X, ma, chars, s ){

  # Select x, y, z positions
  xi <- seq( from = 1, to = chars, by = 3 )
  yi <- seq( from = 2, to = chars, by = 3 )
  zi <- seq( from = 3, to = chars, by = 3 )

  # Fill in array k x q x s
  for (i in 1:s) {
    ma[,1,i] <- unlist( X[i,xi] )
    ma[,2,i] <- unlist( X[i,yi] )
    ma[,3,i] <- unlist( X[i,zi] )
  }

  # Return array
  return( ma )

}

# Convert 2D landmarks matrix into matrix
.coords2D.mat2arr <- function( X, ma, chars, s ){

  # Select x, y positions
  xi <- seq( from = 1, to = chars, by = 2 )
  yi <- seq( from = 2, to = chars, by = 2 )

  # Fill in array k x q x s
  for (i in 1:s) {
    ma[,1,i] <- unlist( X[i,xi] )
    ma[,2,i] <- unlist( X[i,yi] )
  }

  # Return array
  return( ma )

}



# ------------------------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------------------------- #





#' Convert an array into a matrix
#'
#' @description
#' Convert an array with landmark points into an object of class matrix.
#'
#' @param X Array, \code{k} landmark points, \code{q} coordinates, and
#' \code{s} specimens.
#'
#' @param coords Integer, 2 or 3, for 2D or 3D landmarks, respectively.
#'
#' @details
#'
#' The object \code{X}, class array, has format \code{k x q x s}, where \code{k} is
#' the number of landmarks, \code{q} the number of coordinates, and \code{s} the number
#' of specimens. See \code{C.arr.unal} for an example of the format of a 3D array
#' and \code{data-raw/C.R} for the details about how to generate this object.
#'
#' @return
#'
#' An object of class matrix, with \code{s} rows, one for specimen, and \code{n} columns, one
#' for each coordinate of the landmarks.
#' Each landmark can be given in 2D or 3D. For instance,
#' if the landmarks are 3D, the first 3 columns will be the
#' coordinates x, y, and z for the first landmark, the next 3
#' columns for the second landmark, and so on.
#' \tabular{cccccccc}{
#'  specimens \tab lmk1.x \tab lmk1.y  \tab lmk1.z  \tab lmk2.x \tab lmk2.y \tab lmk2.z  \tab ... \cr
#'  Sp_1      \tab 0.143  \tab -0.028  \tab -0.044  \tab 0.129  \tab 0.028  \tab -0.043  \tab ... \cr
#'  Sp_2      \tab 0.128  \tab -0.024  \tab -0.028  \tab 0.124  \tab 0.027  \tab -0.025  \tab ... \cr
#'  ...       \tab ...    \tab ...     \tab ...     \tab ...    \tab ...    \tab ...     \tab ...
#' }
#'
#' See object \code{C} for a more detailed example of the format of the
#' object that is returned and \code{data-raw/C.R} for the explanation
#' about how to generate this object.
#'
#' Note that if the names for the \code{s} matrices, one per specimen, in the array
#' are not provided, i.e. the names for specimens are not given,
#' the specimens in the returned matrix, \code{row.names(matrix)}, will be
#' labelled as '1', '2', and so on.
#'
#' @seealso
#' \code{\link{matrix2array}}, \code{\link{write.morpho}}
#'
#' @author Sandra Alvarez-Carretero and Mario dos Reis
#'
#' @export

array2matrix <- function( X, coords = c( 2, 3 ) ){

  # Check X and coords are provided and their classes
  coords <- as.numeric( match.arg( arg = as.character(coords), choices = c( "2", "3" ) ) )
  .checkCoords( coords )
  .checkArr( X, coords )

  # Create variables
  lmks  <- NULL
  lmks  <- dim( X )[1]
  coords <- coords
  s      <- dim( X )[3]

  # Call function to create array
  faln <- .genMat( X = X, coords = coords, lmks = lmks, s = s )

  # Return matrix
  return( faln )

}


## ################
##  SUBFUNCTIONS ##
## ################

# Check if X is provided and is class "array"
.checkArr <- function( X, coords ){

  # Check object with lmks is provided and it is of class "array"
  if ( missing( X ) ){
    stop( "\nPlease provide an object of class \"array\" with \"k x q x s\" dimensions,
          'k' landmarks, 'q' coordinates, and 's' specimens\n" )
  }
  if ( class( X ) != "array" ){
    stop( "\nPlease use an object of class \"array\"\n" )
  }

  # Check argument coords is alright compared to array X
  if ( coords != dim( X )[2] ){
    stop( "\nThe number of coords does not match the number of
             coordinates in the array you have provided" )
  }

}

# Wraper around .coords3D.arr2mat and .coords2D.arr2mat to create the
# array with landmarks
.genMat <- function( X, coords, lmks, s ){

  # Create an empty matrix to store the coordinates in the format
  # s x n (specimens x characters )
  faln <- matrix( 0, nrow = s, ncol = lmks * coords )

  # Select x, y, z positions
  if ( coords == 3 ){
    faln <- .coords3D.arr2mat( X = X , faln = faln, lmks = lmks, s = s )
  }

  else if ( coords == 2 ){
    faln <- .coords2D.arr2mat( X = X , faln = faln, lmks = lmks, s = s )
  }

  # Return array with positions x, y, and z
  return( faln )

}

# Convert 3D landmarks into matrix
.coords3D.arr2mat <- function( X, faln, lmks, s ){

  # Select x, y, z positions
  xi <- seq( from = 1, to = lmks * 3, by = 3 )
  yi <- seq( from = 2, to = lmks * 3, by = 3 )
  zi <- seq( from = 3, to = lmks * 3, by = 3 )

  # Fill in matrix and put rownames
  for ( i in 1:s ) {
    faln[i,xi] <- X[,1,i]
    faln[i,yi] <- X[,2,i]
    faln[i,zi] <- X[,3,i]
  }
  rownames( faln ) <- dimnames( X )[[ 3 ]]

  # Put colnames
  colnames( faln ) <- paste( "lmk", seq( 1:(lmks*3) ), sep = "" )
  colnames( faln )[xi] <- paste( "lmk", seq( 1:lmks ), ".x", sep = "" )
  colnames( faln )[yi] <- paste( "lmk", seq( 1:lmks ), ".y", sep = "" )
  colnames( faln )[zi] <- paste( "lmk", seq( 1:lmks ), ".z", sep = "" )

  # Return matrix
  return( faln )

}

# Convert 2D landmarks into matrix
.coords2D.arr2mat <- function( X, faln, lmks, s ){

  # Select x and y positions
  xi <- seq( from = 1, to = lmks * 2, by = 2 )
  yi <- seq( from = 2, to = lmks * 2, by = 2 )

  # Fill in matrix and put rownames
  for ( i in 1:s ) {
    faln[i,xi] <- X[,1,i]
    faln[i,yi] <- X[,2,i]
  }
  rownames( faln ) <- dimnames( X )[[ 3 ]]

  # Put colnames
  colnames( faln ) <- paste( "lmk", seq( 1:(lmks*2) ), sep = "" )
  colnames( faln )[xi] <- paste( "lmk", seq( 1:lmks ), ".x", sep = "" )
  colnames( faln )[yi] <- paste( "lmk", seq( 1:lmks ), ".y", sep = "" )

  # Return matrix
  return( faln )

}


# ------------------------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------------------------- #



#' Simulate a continuous morphological alignment
#'
#' @description
#' Simulate a continuous morphological alignment using \code{\link[ape]{rTraitCont}}
#' and later allowing to account for population noise and correlation among characters.
#'
#' @param tree Phylo, object with a phylogenetic tree
#' (see \code{\link[ape]{rTraitCont}}).
#'
#' @param n Numeric, number of morphological traits to be simulated.
#'
#' @param c (Optional) numeric, vector with variances for the speccies within
#' a population to add as population noise to the simulated morphological
#' traits (see details).
#'
#' @param R (Optional) matrix, correlation matrix (see details).
#'
#' @param ... Further arguments passed to \code{\link[ape]{rTraitCont}}
#' (see details).
#'
#' @details
#'
#' The function \code{\link[ape]{rTraitCont}} simulates continuous traits and
#' can take different parameters to adjust the simulation
#' (e.g. the model, the rate drift, etc.).
#' These parameters are the ones the user can pass to the argument \code{...} in
#' \code{sim.morpho}.
#' The default values that \code{sim.morpho} uses are
#' \code{model = "BM"}, \code{sigma = 1}, \code{ancestor = F},
#' and \code{root.value = 0}. For this kind of simulation,
#' \code{sim.morpho} allows only \code{ancestor = F}, so please
#' do not change this parameter.
#' In the \code{\link[ape]{rTraitCont}} package, the parameter
#' \code{model} can be \code{model = BM}, \code{model = OU}, or
#' a function \code{model = FUN} provided by the user. Currently,
#' \code{sim.morpho} supports only the first two.
#'
#' The parameter \code{c} contains the population noise, which is used to simulate
#' the noise matrix. Each parameter follows a normal distribution, \code{x ~ N(0,c)}.
#' If the variances are assumed to be the same for all characters within the species,
#' then the length of \code{c} is 1 and equals to the value of this variance.
#' If it differs, then a vector of length \code{n} has to
#' be provided specifying the variance for each of the characters.
#'
#' The simulated noise is later added to the morphological data previously
#' generated, so we obtain the noisy matrix.
#' If a correlation matrix, \code{R}, is provided, then it is added to the
#' noisy matrix. See object \code{sim.R} for an example of its format and
#' \code{data-raw/sim.R} to understand how it can be generated.
#' Note that the correlation matrix needs to be of class "matrix"
#' and symmetric.
#'
#' @return
#'
#' \item{M}{Matrix with the simulated morphological continuous data accounting
#' for noise and, if provided, for population variance and trait correlation too.}
#'
#' @seealso
#' \code{\link{write.morpho}}
#'
#' @author Sandra Alvarez-Carretero and Mario dos Reis
#'
#' @examples
#' # A) Simulation setup: Simulate a morphological alignment
#' #    with n = 100 continuous characters for a phylogeny
#' #    defined in object 'tree', with the default parameters in
#' #    'sim.morpho' to run 'rTraitCont'.
#' #    Population noise and character correlation are not considered,
#' #    i.e. c = 0 and R not provided.
#'
#'      sim.morpho( tree = sim.tree, n = 100 )
#'
#'
#' # B) Simulation setup: Simulate a morphological alignment
#' #    with n = 100 continuous characters for a phylogeny
#' #    defined in object 'tree', but with different parameters
#' #    than the default ones in 'sim.morpho' to run
#' #    'rTraitCont'. Population noise and trait correlation are not
#' #    considered, i.e. c = 0 and R not provided.
#'
#'      sim.morpho( tree  = sim.tree, n = 100,
#'                          model = "OU", sigma = 0.2, alpha = 2 )
#'
#'
#' # C) Simulation setup: Simulate a morphological alignment
#' #    with n = 100 continuous characters for a phylogeny
#' #    defined in object 'tree', with the default parameters in
#' #    'sim.morpho' to run 'rTraitCont'.
#' #    Population noise is low, c = 0.25, but trait correlation is not
#' #    considered, i.e. R not provided.
#'
#'      sim.morpho( tree = sim.tree, n = 100, c = 0.25 )
#'
#'
#' # D) Simulation setup: Simulate a morphological alignment
#' #    with n = 100 continuous characters for a phylogeny
#' #    defined in object 'tree', with the default parameters in
#' #    'sim.morpho' to run 'rTraitCont'.
#' #    Population noise is low, c = 0.25, and a correlation
#' #    matrix to simulate trait correlation (rho = 0.50) is provided.
#'
#'      sim.morpho( tree = sim.tree, n = 100, c = 0.25, R = sim.R )
#'
#' @export

sim.morpho <- function( tree, n , c = 0, R, ... ) {

  # Check objects tree and n are provided and that c !< 0
  .checkFiles( tree = tree, n = n, c = c )

  # Get ... parameters in a list in case the user has provided them.
  # Else, assign default values
  pars <- .checkRtraitCont( rtraitcont.pars = list( ... ) )

  # Get amount of specimens
  s <- length(tree$tip.label)

  # Simulate continous data with ape::rTraitCont
  M <- .simMorphData( tree = tree,  s = s, n = n, pars = pars )

  # Simulate noise matrix, which follows a normal distribution
  # with mu = 0 and sd = c
  N <- .simNoise( s = s, n = n, c = c )

  # Add character correlation if R matrix is provided
  if ( missing( R ) ){
    # Generate noisy matrix
    M.n <- M + N
  }
  else if ( ! missing( R ) ){
    .checkCorrMat( R = R, n = n )
    # Generate noisy matrix with correlation
    M.n <- .addCorrMat( X = M, Y = N, R = R )
  }

  # Return the M.n matrix (with or without corr), which
  # will be transformed when
  # written in the output file by write.morpho()
  return( M.n )

}


## ################
##  SUBFUNCTIONS ##
## ################

# Checking tree is provided and is class "phylo".
# Also check n are c are provided and that "c" is either
# a unique value or a vector length "n"
.checkFiles <- function( tree, n, c ){

  if ( missing( tree ) ){
    stop( "\nPlease provide an object of class \"phylo\" with a phylogenetic tree\n" )
  }

  if ( class( tree ) != "phylo" ){
    stop( "\nPlease use an object of class \"phylo\" with a phylogenetic tree\n" )
  }

  if( missing( n ) ){
    stop( "\nPlease, provide a numerical value with the amount of
          morphological continuous characters to be simulated for
          the population\n" )
  }
  if( as.integer( n ) != n | n < 0 ){
    stop( "\nPlease, provide a numerical value with the amount of
          morphological continuous characters to be simulated for
          the population, i.e. n > 0 and not a decimal\n" )
  }

  # Check population noise
  if ( missing( c ) ){
    stop( "\nPlease, a vector of variances of length equal to the number of
          morphological characters, n, or to a unique numeric value,
          which will assume that all variances equal to this value\n" )
  }
  if( class( c ) != "numeric" ){
    stop( "\nThe parameter \"c\" needs a numeric value > 0 or a numeric vector
          of length \"n\" with values > 0\n" )
  }
  if( length( c ) != n & length( c ) != 1 ){
    stop( "\nProvide a vector with variances of length equal to the number of
          morphological characters, \"length( c ) = n\", or to a unique numeric value,
          \"length( c ) = 1\", which will assume that all variances equal to this value\n" )
  }
  # Check any value in the vector of vars is <0
  if( length( c ) == n & any( c < 0 ) ){
    stop( "\nThe parameter \"c\" needs all numeric values > 0\n" )
  }
  # Check the popvar value is >0
  if( length( c ) == 1 & any( c < 0 ) ){
    stop( "\nThe parameter \"c\" needs a numeric value\n" )
  }


}


# Checking if optional values to pass to ape::rTraitCont are provided
.checkRtraitCont <- function( rtraitcont.pars ){

  pars <- list( model = "BM", sigma = 1, alpha = 1,
                theta = 0, ancestor = FALSE,
                root  = 0 )

  if ( "model" %in% names( rtraitcont.pars ) ){
    pars$model <- rtraitcont.pars$model
  }

  if ( "sigma" %in% names( rtraitcont.pars ) ){
    pars$sigma <- rtraitcont.pars$sigma
  }

  if ( "alpha" %in% names( rtraitcont.pars ) ){
    pars$alpha <- rtraitcont.pars$alpha
  }

  if ( "theta" %in% names( rtraitcont.pars ) ){
    pars$theta <- rtraitcont.pars$theta
  }

  if ( "ancestor" %in% names( rtraitcont.pars ) ){
    pars$ancestor <- rtraitcont.pars$ancestor
  }

  if ( "root.value" %in% names( rtraitcont.pars ) ){
    pars$root <- rtraitcont.pars$root.value
  }

  # Return list with updated values
  return( pars )

}


# Simulate continous data with ape::rTraitCont and parameters in pars
.simMorphData <- function( tree,  s, n, pars ){

  # Simulate matrix s x n (species x characters)
  # with ape:rTraitCont
  M <- replicate( n,
                  ape::rTraitCont( phy      = tree,          model      = pars$model,
                                   sigma    = pars$sigma,    root.value = pars$root,
                                   ancestor = pars$ancestor, alpha      = pars$alpha,
                                   theta    = pars$theta ) )
  #Return matrix M
  return( M )

}


# Simulate noise matrix
# Sample num.species*traits samples from a normal distribution
# distrib with mu = 0 and sd = c
.simNoise <- function( s, n, c ){

  N <- t( replicate( s, rnorm( n, mean = 0, sd = sqrt( c ) ) ) )

  # Return noise matrix
  return( N )

}


# Checking correlation matrix, if provided, is class = "matrix",
# and symmetric
.checkCorrMat <- function( R, n ){

  if ( class( R ) != "matrix" ){
    stop( "\nPlease use a correlation matrix of class \"matrix\"\n" )
  }

  if ( dim( R )[1] != n & dim( R )[2] != n ){
    stop( "\nPlease use a symmetric correlation matrix of size \"n x n\"\n" )
  }

  if ( isSymmetric( R ) != TRUE ){
    stop( "\nPlease provide a symmetric correlation matrix of class \"matrix\"\n" )
  }

}

# [[ NOT USED NOW AS IT MIGHT BE TOO TIME CONSUMING IF MATRICES
#                         ARE TOO BIG                           ]]
# Checking correlation matrix, if provided, is class = "matrix",
# positive definite, and symmetric
# .checkCorrMat2 <- function( R, n ){
#
#   if ( class( R ) != "matrix" ){
#     stop( "\nPlease use a correlation matrix of class \"matrix\"\n" )
#   }
#
#   if ( dim( R )[1] != n & dim( R )[2] != n ){
#     stop( "\nPlease use a correlation matrix of size \"n x n\"\n" )
#   }
#
#   if ( isSymmetric( R ) != TRUE ){
#     stop( "\nPlease provide a symmetric correlation matrix of class \"matrix\"\n" )
#   }
#
#   # Check R is a positive definite matrix -- This might take a while to check
#   # depending on R size
#   eigenvals     <- eigen( R )$values
#   neg.eigenvals <- which ( eigenvals < 0 )
#   if ( length( neg.eigenvals ) != 0 ){
#     stop( "\nPlease provide a positive definite correlation matrix of class \"matrix\"\n" )
#   }
#
# }


# Add correlation to noisy matrix
.addCorrMat <- function( X, Y, R ){

  # Check the class of R is "matrix"
  # [ALREADY CHECKED IN MAIN FUNC]
  #.checkCorrMat( R )

  # Code: X = M, Y = N, R = R

  # Get correlated data
  MR <- X %*% chol( R )

  # Add correlation to noise matrix
  NR <- Y %*% chol( R )

  # Get matrix with correlated data and pop noise
  M.n.R <- MR + NR

  # Return M.n.R
  return( M.n.R )

}


# ------------------------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------------------------- #

#' Simulate a population matrix
#'
#' @description
#' Simulate a population sample and return a list with (i) a matrix of
#' size \code{s x n}, \code{s} specimens and \code{n} characters, (ii) a
#' vector with the estimated population variances for each character,
#' and (iii) the estimated shrinkage correlation matrix if the true
#' correlation matrix is provided.
#'
#' @param psample Numeric, number of specimens the simulated
#' population sample should include.
#'
#' @param n Numeric, number of morphological traits to be simulated.
#'
#' @param c Numeric, vector with the variances for the species 
#' within a population (see details).
#'
#' @param R (Optional) matrix, correlation matrix.
#' (see details).
#'
#' @return
#'
#' \item{$P}{Matrix with the simulated population sample}
#' \item{$var}{Vector with the estimated variances}
#' \item{$Rsh}{Estimated shrinkage correlation matrix, only returned if \code{R} is provided}
#'
#' @details
#'
#' The parameter \code{c} is the population noise and it is used to sample
#' \code{n} characters for each of the \code{psample} specimens from a
#' normal distribution \code{x ~ N(0,c)}.
#' If the population noise is assumed to be the same for all the characters
#' within the species, then the length of \code{c} is 1 and equals to the value
#' of this variance.
#' If it differs, then a vector of length \code{n} has to
#' be provided specifying the variance for each of the characters
#'
#' If a correlation matrix, \code{R}, is provided, then it is added to the
#' population matrix. Note that the correlation matrix needs to be of class "matrix"
#' and symmetric. You can take a look at \code{data-raw/sim.R.R} to follow
#' the commands used to generate this matrix, object \code{R.sim}, which is used
#' in the examples.
#'
#' @seealso
#' \code{\link{sim.morpho}}, \code{\link{write.morpho}}
#'
#' @author Sandra Alvarez-Carretero and Mario dos Reis
#'
#' @examples
#'
#' # A) Simulation setup: Simulate a population with
#' #    psample = 20 specimens, and sample n = 100 characters with
#' #    a low population noise, c = 0.25.
#'
#'      sim.pop( psample = 20, n = 100, c = 0.25 )
#'
#' # B) Simulation setup: Simulate a population with
#' #    psample = 20 specimens, and sample n = 100 characters with
#' #    a low population noise, c = 0.25, and a low trait correlation
#' #    rho = 0.50 (correlation matrix that follows
#' #    the constant correlation model, i.e. all non-diagonal values
#' #    equal to rho).
#'
#'      sim.pop( psample = 20, n = 100, c = 0.25, R = sim.R )
#'
#' @export

sim.pop <- function( psample, n, c, R ){

  # Check input parameters
  .checkInSimPop( psample = psample, n = n, c = c )

  # Sample from a normal distribution the morph. cont. traits
  P <- t( replicate( psample, rnorm( n, mean = 0, sd = sqrt( c ) ) ) )

  # Add correlation if provided
  if ( ! missing( R ) ){
    .checkCorrMat( R = R, n = n )
    P    <- P %*% chol( R )
    R.sh <- as.matrix( corpcor::cor.shrink( P ) )
  }

  # Get pop variance
  var.P <- diag( cov( P ) )

  # Get names for the specimens as "Species_1", "Species_2", and so on
  vars          <- .notNames( s = psample )
  rownames( P ) <- vars$names

  # Return P
  if ( missing( R ) ){
    return( list( P = P, var = var.P, Rsh = R.sh ) )
  }else
    return( list( P = P, var = var.P ) )

}

## ################
##  SUBFUNCTIONS ##
## ################

# Check input values
.checkInSimPop <- function( psample, n, c ){

  if( missing( psample ) ){
    stop( "\nPlease, provide a numerical value with the amount of
          specimens to be included in the population\n" )
  }
  if( class( psample ) != "numeric" ){
    stop( "\nPlease, provide a numerical value with the amount of
          specimens to be included in the population\n" )
  }

  if( missing( n ) ){
    stop( "\nPlease, provide a numerical value with the amount of
          morphological continuous characters to be simulated for
          the population\n" )
  }
  if( as.integer( n ) != n | n < 0 ){
    stop( "\nPlease, provide a numerical value with the amount of
          morphological continuous characters to be simulated for
          the population, i.e. n > 0 and without decimal part\n" )
  }

  # Check population variance
  if ( missing( c ) ){
    stop( "\nPlease, a vector of variances of length equal to the number of
          morphological characters, n, or to a unique numeric value,
          which will assume that all variances equal to this value\n" )
  }
  if( class( c ) != "numeric" ){
    stop( "\nThe parameter \"c\" needs a numeric value > 0 or a numeric vector
          of length \"n\" with values > 0\n" )
  }
  if( length( c ) != n & length( c ) != 1 ){
    stop( "\nProvide a vector with variances of length equal to the number of
          morphological characters, \"length( c ) = n\", or to a unique numeric value,
          \"length( c ) = 1\", which will assume that all variances equal to this value\n" )
  }
  # Check any value in the vector of vars is <0
  if( length( c ) == n & any( c < 0 ) ){
    stop( "\nThe parameter \"c\" needs all numeric values > 0\n" )
  }
  # Check the popvar value is >0
  if( length( c ) == 1 & any( c < 0 ) ){
    stop( "\nThe parameter \"c\" needs a numeric value\n" )
  }

}

# ------------------------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------------------------- #

#' Procrustes alignment output in MCMCtree format for Bayesian inference
#'
#' @description
#' Perform a Procrustes alignment with the \code{\link[Morpho]{procSym}}
#' and \code{\link[Morpho]{align2procSym}} and generate a morphological
#' alignment in MCMCtree format.
#'
#' @param data Matrix or array, object with the data set with one specimen per
#' species. See details.
#'
#' @param popdata Matrix or array, containing the specimens for the population of
#' species with a larger variation (more than one specimen for this species).
#' See details.
#'
#' @param sp.data Numeric, position of the specimen in the \code{data} object
#' also present in the \code{popdata} object (see details).
#' 
#' @param sp.popdata Numeric, position of the specimen in the \code{popdata} object
#' also present in the \code{data} object (see details).
#'
#' @param filename Character, name for the output file.
#'
#' @param coords Numeric, 2 or 3 for 2D or 3D landmarks, respectively.
#'
#' @param method (Optional) character, either \code{"eigen"} or
#' \code{"col"}, method used to decompose the inverse of the shrinkage
#' correlation matrix. See details.
#'
#' @param ages (Optional) list or vector, ages of the species included in the
#' morpholical alignment.
#'
#' @param ... (Optional) Further arguments passed to \code{\link[Morpho]{procSym}}
#' (see details)
#'
#' @details
#'
#' If the objects \code{data} and \code{popdata} are of class "matrix", then they have
#' \code{s} or \code{ps} rows, respectively, and \code{n} columns as of the amount of characters.
#' These data sets are supposed to contain landmarks, which can be given in 2D or 3D.
#' For instance, if the landmarks are 3D, the first 3 columns will be the
#' coordinates x, y, and z for the first landmark, the next 3
#' columns for the second landmark; and so on.
#' \tabular{cccccccc}{
#'  specimens \tab lmk1.x \tab lmk1.y  \tab lmk1.z  \tab lmk2.x \tab lmk2.y \tab lmk2.z  \tab ... \cr
#'  Sp_1      \tab 0.143  \tab -0.028  \tab -0.044  \tab 0.129  \tab 0.028  \tab -0.043  \tab ... \cr
#'  Sp_2      \tab 0.128  \tab -0.024  \tab -0.028  \tab 0.124  \tab 0.027  \tab -0.025  \tab ... \cr
#'  ...       \tab ...    \tab ...     \tab ...     \tab ...    \tab ...    \tab ...     \tab ...
#' }
#' You can load the object \code{C.mat.unal} to see an example of this format and also
#' take a look at the \code{data-raw/C.R} file for details about how to
#' generate it. Otherwise, if the objects \code{data} and \code{popdata} are of
#' class "array", they have format \code{k x q x s} or \code{k x q x ps}, respectively,
#' where \code{k} is the number of landmarks, \code{q} the number of coordinates,
#' \code{s} the number of specimens for object \code{data}, and \code{ps} the number of
#' specimens in a sampled population of one species for object \code{popdata}.
#' You can load the object \code{C.arr.unal} to see an example of this array format and
#' also take a look at the \code{data-raw/C.R} file for details about how to
#' generate it.
#'
#' The names of the specimens will be taken from either the row names of the data
#' matrices, if \code{data} and \code{popdata} objects are of class "matrix", or from
#' the names of the third dimension of the array, if these objects are of class
#' "array". If no names are found, the specimens will be labelled as "Specimen_1",
#' "Specimen_2", and so on.
#'
#' Note that you are providing a population matrix with the landmarks collected from
#' more than on specimen belonging to the same species. Therefore, it is assumed that
#' population noise is accounted for.
#' Furthermore, the landmarks of one of these specimens are also present in the
#' morphological alignment. First, a Procrusts alignment is performed with the data set
#' with one specimen per species (\code{data}), and then all the specimens in \code{popdata}
#' except for the specimen common in \code{data}, which position is specified in the
#' argument \code{sp.popdata}, are aligned to the mean shape of the
#' PA generated with \code{data}. Take a look at \code{data.raw/V.R} for more details.
#'
#' The logarithm of the determinant of the correlation matrix is going to be printed
#' in the output file to later be used by MCMCtree during the likelihood calculation.
#'
#' The function \code{\link[Morpho]{procSym}} can perform a Procrustes superimposition
#' alignment given a dataset of landmarks in array format and allows the user to pass
#' different options to the arguments defined.
#' The current function runs with the default parameters
#' in \code{\link[Morpho]{procSym}} but allows the user to pass three arguments to
#' \code{\link{proc2MCMCtree}} through \code{...}: scale, reflect, and pairedLM.
#' If you are thinking of using these arguments, read the documentation in \code{\link[Morpho]{procSym}},
#' otherwise do not pass further arguments to \code{\link{proc2MCMCtree}} and let it run
#' \code{\link[Morpho]{procSym}} in default mode.
#'
#' @return
#'
#' \item{$dataPS}{Object of class \code{symproc} output by the function
#' \code{Morpho::procSym}. The user can access the rest of the items in this object
#' by checking <your_object>$dataPS$<available_Morpho_objects>. It contains the
#' array with the morphological alignment generated by this function, which
#' corresponds to the specimens in the object passed to \code{data}. This can
#' be accessed by loading <your_object>$dataPS$rotated.}
#' \item{$popdataPS}{Array with the morphological alignment with the object
#' passed to \code{popdata}, i.e. the alignment with the population sample}
#' \item{$M}{Matrix with the morphological alignment with the specimens of all
#' species. Note that this alignment is not corrected for character correlation
#' nor population noise. The corrected alignment is only printed in the output alignment file.}
#' \item{$Rsh}{Estimated shrinkage correlation matrix.}
#' \item{$c}{Estimated population variance.}
#'
#' @seealso
#' \code{\link{write.morpho}}, \code{\link{matrix2array}}, \code{\link{array2matrix}}
#'
#' @author Sandra Alvarez-Carretero and Mario dos Reis
#'
#' @examples
#'
#' # A. Use the unaligned, but processed, carnivoran data (data = C.mat.unal) and
#' #    vulpes data (popdata = V.mat.unal) to obtain a morphological alignment.
#' #    The fox specimen that is common in V.mat.unal and C.mat.unal is the
#' #    one in the first row of V.mat.unal (sp.popdata = 1) and the one in position
#' #    13 in C.mat.unal (sp.data = 13). The method to
#' #    decompose the estimate shrinkage correlation matrix (internally calculated
#' #    in this function) is the Cholesky decomposition (method = c("chol")).
#' #    We do not add ages.
#'      right    <- c( 11, 22, 13, 19, 15, 20, 24, 5, 7, 2, 9, 26, 29 )
#'      left     <- c( 10, 21, 12, 17, 14, 18, 23, 4, 6, 1, 8, 25, 28 )
#'      pairedLM <- cbind( left, right )
#'      obj.aln <- proc2MCMCtree( data = C.mat.unal, popdata = V.mat.unal, sp.data = 13,
#'                                sp.popdata = 1, filename = "./seqfile.aln", coords = 3,
#'                                method = c("chol"), pairedLM = pairedLM )
#'
#' # B. Use the unaligned, but processed, carnivoran data (data = C.mat.unal) and
#' #    vulpes data (popdata = V.mat.unal) to obtain a morphological alignment.
#' #    The fox specimen that is common in V.mat.unal and C.mat.unal is the
#' #    one in the first row of V.mat.unal (sp.popdata = 1) and the one in position
#' #    13 in C.mat.unal (sp.data = 13). The method to
#' #    decompose the estimate shrinkage correlation matrix (internally calculated
#' #    in this function) is the Cholesky decomposition (method = c("chol")).
#' #    We add ages.
#'      ages <- list( sp1  = 13.135, sp2  =  0.0285, sp3  = 11.95,  sp4  = 35.55,
#'                    sp5  = 25.615, sp6  = 14.785,  sp7  = 28.55,  sp8  =  0,
#'                    sp9  =  0,     sp10 =  0,      sp11 =  0,     sp12 =  0,
#'                    sp13 =  0,     sp14 =  0,      sp15 =  0,     sp16 =  0,
#'                    sp17 =  6.65,  sp18 =  0.0285, sp19 =  0 )
#'      right    <- c( 11, 22, 13, 19, 15, 20, 24, 5, 7, 2, 9, 26, 29 )
#'      left     <- c( 10, 21, 12, 17, 14, 18, 23, 4, 6, 1, 8, 25, 28 )
#'      pairedLM <- cbind( left, right )
#'      obj.aln <- proc2MCMCtree( data = C.mat.unal, popdata = V.mat.unal, sp.data = 13,
#'                                sp.popdata = 1, filename = "./seqfile.aln", coords = 3,
#'                                method = c("chol"), pairedLM = pairedLM, ages = ages )
#'
#' @export
proc2MCMCtree <- function( data, popdata, sp.data, sp.popdata, filename, coords, method = c( "eigen", "chol" ), ages, ... ){

  # Check data and popdata are either array or matrix class
  .checkData( data = data, popdata = popdata, sp.data = sp.data, sp.popdata = sp.popdata )

  # Get ... parameters in a list in case the user has passed
  # arguments for Morpho::procSym function
  pars <- .checkProcSym( procsym.pars = list( ... ) )

  # Create vars for specimens names
  dat.names <- NULL
  pop.names <- NULL

  # If data and/or popdata are a matrix, convert them into an array
  if ( class( data ) == "matrix" ){
    data    <- mcmc3r::matrix2array( X = data, coords = coords )
    if( ! is.null( row.names( data ) ) ){
      dat.names <- row.names( data )
    }else{
      dat.names <- paste( "Specimen_", seq(1:dim(data)[1]), sep = "" )
    }
  }else{
    if( ! is.null( dimnames( data )[[3]] ) ){
      dat.names <- dimnames( data )[[3]]
    }else{
      dat.names <- paste( "Specimen_", seq(1:length(dimnames( data )[[3]])), sep = "" )
    }
  }
  if( class( popdata ) == "matrix" ){
    popdata <- mcmc3r::matrix2array( X = popdata, coords = coords )
    if( ! is.null( row.names( popdata ) ) ){
      pop.names <- row.names( popdata )
    }else{
      pop.names <- paste( "P.Specimen_", seq(1:dim(popdata)[1]), sep = "" )
    }
  }else{
    if( ! is.null( dimnames( popdata )[[3]] ) ){
      pop.names <- dimnames( popdata )[[3]]
    }else{
      pop.names <- paste( "Specimen_", seq(1:length(dimnames( popdata )[[3]])), sep = "" )
    }
  }

  # Create variables
  k  <- length( dimnames( data )[[1]] )
  q  <- length( dimnames( data )[[2]] )
  s  <- length( dimnames( data )[[3]] )
  ps <- length( dimnames( popdata )[[3]] )
  method <- as.character( match.arg( arg = as.character( method ),
                                     choices = c( "eigen", "chol" ) ) )

  # Perform Procrustes analysis with Morpho package
  # 1. Align specimens in data
  data.PS    <- Morpho::procSym( dataarray = data,
                                 scale = pars$scale, reflect = pars$reflect,
                                 pairedLM = pars$pairedLM )
  # 2. Align all the specimens in popdata except for the common specimen in
  #    data to the mean shape of the "data.PS" alignment so they have the same
  #    orientation
  #    Get also popdata.PS.filt in matrix format
  popdata.filt           <- popdata[,,-sp.popdata]
  popdata.PS.filt        <- Morpho::align2procSym( x = data.PS,
                                                   newdata = popdata.filt )
  dimnames( popdata.PS.filt ) <- dimnames( popdata.filt )

  # 3. Get the common specimen in popdata.filt and data.PS and add it to
  #    "P" object.
  P <- .alnPopDat( k = k, q = q, s = ps, pos.sp.dat = sp.data,
                   name.sp = rownames( popdata )[ sp.popdata ],
                   dataPS = data.PS, popdataPS = popdata.PS.filt )
  M <- mcmc3r::array2matrix( X = data.PS$rotated, coords = 3)
  
  # Estimate Rshrinkage with corpcor::cor.shrink
  R.sh <- corpcor::cor.shrink( P )
  R.sh <- cbind( R.sh )

  # Transform data set accounting for popvar and corr using
  # write.morpho() and get the alignment file in MCMCtree format

  mcmc3r::write.morpho( M = M, filename = filename,
                        c = diag( cov( P ) ), R = R.sh,
                        method = method, names = row.names( M ), ages = ages )

  # Return list object with Procrustes alignments, morpho matrices,
  # pop variances, R shrinkage matrix, and number of specimen from popdata
  # included in final morpho alignment
  return( list( dataPS = data.PS, popdataPS = P,
          M = M, Rsh = R.sh, c = diag( cov( P ) ) ) )

}


## ################
##  SUBFUNCTIONS ##
## ################

# Check data and popdata are either of class array or matrix
# If specimen is provided, check it is of class numeric
.checkData <- function( data, popdata, sp.data, sp.popdata ){

  # Check object with lmks is provided and it is of class "matrix"
  if ( missing( data ) | missing( popdata ) ){
    stop( "\nPlease provide an object of class \"matrix\" or \"array\" for arguments \"data\" and \"popdata\"\n" )
  }

  if ( class( data ) != "matrix" & class( data ) != "array" ){
    stop( "\nPlease provide an object of class \"matrix\" or \"array\" for argument \"data\"\n" )
  }
  if ( class( popdata ) != "matrix" & class( popdata ) != "array" ){
    stop( "\nPlease provide an object of class \"matrix\" or \"array\" for argument \"popdata\"\n" )
  }

  # Check a name for the output file has been given
  #if ( missing( filename ) ){
  #  stop( "\nPlease use the parameter \"filename\" to provide a name for the output file\n" )
  #}

  # Check specimen is numeric
  if( ! missing( sp.data ) & as.integer( sp.data ) != sp.data ){
    stop( "\nPlease provide a numeric value for object \"sp.data\"\n" )
  }
  if( ! missing( sp.popdata ) & as.integer( sp.popdata ) != sp.popdata ){
    stop( "\nPlease provide a numeric value for object \"sp.popdata\"\n" )
  }

  # Option "coords" is checked by array2matrix and matrix2array
  # Options "filename", "method", and "ages" are checked by write.morpho

}

# Check if user has passed arguments for ProcSym
.checkProcSym <- function( procsym.pars ){

  pars <- list( scale = TRUE, reflect = TRUE, pairedLM = NULL )

  if( "scale" %in% names( procsym.pars ) ){
    pars$scale <- procsym.pars$scale
  }
  if( "reflect" %in% names( procsym.pars ) ){
    pars$reflect <- procsym.pars$reflect
  }
  if( "pairedLM" %in% names( procsym.pars ) ){
    pars$pairedLM <- procsym.pars$pairedLM
  }

  # Return list with updated values
  return( pars )
}


# Get alignment with "data" and one randomly sampled specimen from "popdata"
# as morphological alignment "M" and then the alignment for the population sample
# within "popdata"
.alnPopDat <- function( k, q, s, pos.sp.dat, name.sp, dataPS, popdataPS ){

  # Create a new array to save the morpho alignment
  combPS <- array( dim = c( k, q, s ) )

  # Add data depending on how many coords are
  if ( q == 2 ){
    dimnames( combPS ) <- list( paste( "lmk", seq( 1:k ), sep="" ),
                                c( "x", "y" ),
                                c( name.sp,
                                   dimnames( popdataPS )[[3]] ) )
  }
  else if( q == 3 ){
    dimnames( combPS ) <- list( paste( "lmk", seq( 1:k ), sep="" ),
                                c( "x", "y", "z" ),
                                c( name.sp,
                                   dimnames( popdataPS )[[3]] ) )
  }
  combPS[,,1]     <- dataPS$rotated[,,pos.sp.dat]
  combPS[,,(2:s)] <- popdataPS

  # Convert array into matrix
  combPS.mat <- mcmc3r::array2matrix( X = combPS, coords = q )

  # Return morpho aln in matrix format
  return( combPS.mat )

}



# ------------------------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------------------------- #


#' Generating a tree file in MCMCtree format when using
#' morphological data
#'
#' @description
#' Outputs a tree file to be used in MCMCtree. The path to this file
#' needs to be written next to the option "treefile = " in the
#' control file for MCMCtree.
#' @param tree Character, path to the file with the tree topology
#' in Newick format, without branch lengths (see details).
#' @param aln Character, path to the alignment file output by
#' \code{\link{write.morpho}} or \code{\link{proc2MCMCtree}}.
#' @param filename Character, name for the output file.
#'
#' @details
#'
#' If you used \code{\link{write.morpho}} or \code{\link{proc2MCMCtree}}
#' and passed a vector or a list with the ages of the specimens in your
#' alignment, you will see that the names in the output file with this
#' alignment are followed by a "^" and a value. This value is transformed
#' to the age you input within MCMCtree, as if it was a tip date, and used
#' to estimate the divergence times of the species in your phylogeny.
#' This function generates a tree file with the same species names than in
#' the morphological alignment, i.e. with the ages to be used by MCMCtree
#' next to the names of each specimen.
#' Remember to write the path to the tree file next to the "treefile = "
#' option in the control file to run MCMCtree.
#' Remember that the names without "^" followed by the ages have to be the
#' same in the file you pass to parameter "tree" than the ones you have
#' in the alignment file passed to "aln".
#'
#' @examples
#'
#' # Use the file with the tree topology (no branch lengths) and the
#' # file with the morphological alignment that are saved in the
#' # inst/extdata directory to generate a tree file with the
#' # ages used by MCMCtree included.
#' # We call the output tree file "treefile.txt".
#'    tree <- system.file( "extdata", "19s.trees", package = "mcmc3r")
#'    aln  <- system.file( "extdata", "seqfile.aln", package = "mcmc3r")
#'    treeMCMCtree( tree = tree, aln = aln, filename = "treefile.txt" )
#'
#' @export

treeMCMCtree <- function(tree, aln, filename){

  # Check parameters
  .checkTree( tree = tree, aln = aln , filename = filename )

  # Read tree and alignment files
  tt  <- readLines( tree )
  inp <- read.table( file = aln, skip = 3, dec = " ", header = F, stringsAsFactors = F )

  # Find the pattern (species names), and replace it with the species name together with the
  # age in MCMCtree format
  tt <- .matchNames( inp = inp, tt = tt )

  # Write header info for MCMCtree (species nºtrees)
  .writeTree( inp = inp, filename = filename, tt = tt )

}


# ###############
# SUBFUNCTIONS ##
#################

.checkTree <- function( tree, aln, filename ){

  if( missing( tree ) ){
    stop( "\nPlease, type the path to the tree file next to the
          parameter \"tree\"\n" )
  }
  if( missing( aln ) ){
    stop( "\nPlease, type the path to the alignment file next to the
          parameter \"aln\"\n")
  }
  if( missing( filename ) ){
    stop( "\nPlease, introduce a name for your output file next
          to the parameter \"filename\"\n" )
  }

}


.matchNames <- function( inp, tt ){

  # Replace names with names+ages
  for ( i in seq( 1:dim( inp )[1] ) ){
    name.inp = inp[i,1]
    #cat("This is name.inp= ", name.inp, "\n")
    name.sp = gsub( "\\^.*", "", name.inp )
    #print(name.sp)
    tt <- gsub( name.sp, name.inp, tt )
    #cat("This is name.sp =", name.sp, "\n")
    #cat("This is tt=, ", tt, "\n")
  }

  # Return tree with new species labels with ages
  return( tt )

}


.writeTree <- function( inp, filename, tt ){

  # Generate and append header for tree file
  info <- paste( " ", length(inp[,1]), "1", sep = " " )
  write( info, file = filename, append = T )
  #final.tree <- paste( " ", tt, sep = " " )
  # Append tree to tree file
  write( tt, file = filename, append = T )

}


# ------------------------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------------------------- #

#' Generating the control file to run MCMCtree when
#' using morphological data
#'
#' @description
#' Function that outputs the control file to run MCMCtree. It can be used to run
#' only with morphological data or morphological+molecular data
#' (more than one partition in the alignment file) together.
#' @param filename  Character, name for the output control file.
#' @param mol Boolean, TRUE if you include molecular data, FALSE otherwise.
#' Default = FALSE.
#' @param seed Numeric, seed value. Default = -1 (assigns random seed using
#' computer's current time).
#' @param seqfile Character, path to the alignment file.
#' @param treefile Character, path to the tree file.
#' @param mcmcfile Character, path to the file with the report of MCMC runs. By
#' default it generates a file called "mcmc.txt" in the
#' directory where MCMCtree is run.
#' @param outfile Character, path to the summary results file. By default it
#' generates a file called "out.txt" in the directory where
#' MCMCtree is run.
#' @param ndata Numeric, number of partitions in the alignment file.
#' @param seqtype Numeric, 0 for nucleotide sequences, 1 for codon sequences,
#' and 2 for amino acid sequences. Default = 0.
#' @param usedata Numeric, 1: Calculate the likelihood function in the normal
#' way, 0: Likelihood is not calculated (likelihood = 1), 2 and
#' 3: approximate likelihood calculation and ML estimation of
#' branch lengths (see details). Default = 1.
#' @param clock Numeric, 1: strict clock model, 2: independent rates model,
#' 3: autocorrelated rates model (see details).
#' @param RootAge Character, calibration for the root.
#' @param TipDate Numeric, time unit to scale the estimated divergence times.
#' See details.
#' @param alpha Numeric, alpha value for the discrete-gamma model of rate
#' variation. Only used if molecular data is available.
#' @param ncatG Numeric, number of categories for the discrete-gamma model of
#' rate variation. Only used if molecular data is available.
#' @param cleandata Numeric, 0: alignment gaps and ambiguity characters are
#' treated as missing data, 1: any site where at least one
#' sequences has an alignment gap or ambiguity character is
#' deleted (see details). Default = 0.
#' @param BDparas Numeric, vector with the parameters controlling the
#' birth-death-sequential-sampling (BDSS) process (see details).
#' @param kappa_gamma Numeric, vector with the parameters for the substitution
#' model parameter kappa (transition/transversion rate ratio).
#' Only used if molecular data is available.
#' @param alpha_gamma Numeric, vector with the parameters for the substitution
#' model parameter gamma (gamma shape parameter for variable rates among sites).
#' Only used if molecular data is available.
#' @param rgene_gamma Numeric, vector with the parameters for the Dirichlet-gamma
#' prior for the mean substitution rate (see details).
#' @param sigma2_gamma Numeric, vector with the parameters for the Dirichlet-gamma
#' prior for the rate drift parameter (see details).
#' @param print Numeric, 0: results are printed to screen only, 1: MCMC is
#' written to "mcmcfile" and the summary to the "outfile", 2: same as 1 but
#' rates for branches for each partitions are appended to "outfile". Default = 2.
#' @param burnin Numeric, number of iterations to be discarded (burn-in).
#' @param sampfreq Numeric, number of iterations after which a sample will be collected.
#' @param nsample Numeric, number of samples to be gathered.
#' @param model Numeric, substitution model to be used (see details). 0:JC69, 1:K80,
#' 2:F81, 3:F84, 4:HKY85, 5:T92, 6:TN93, 7:REV, 8:UNREST, 9:REVu; 10:UNRESTu.
#' Only used if molecular data is available.
#'
#' @details
#'
#' For more information, please check the MCMCTree tutorial and
#' the PAML documentation.
#'
#' @examples
#'
#' # First create objects with the path to alignment and tree files and then
#' # call the function to generate the control file. The parameters not passed
#' # to the function are used as the default values
#' tree <- system.file( "extdata", "19sMCMCTree.trees", package = "morpho")
#' aln  <- system.file( "extdata", "seqfile.aln", package = "morpho")
#' ctlMCMCTree( filename = "../mcmctree.ctl", mol = FALSE, seqfile = aln, treefile = tree,
#' ndata = 1, clock = 2, TipDate = 1, RootAge = c("B(37.3, 66.0, 0.025, 0.025)"),
#' BDparas = c( 1, 1, 0, 0.001 ), rgene_gamma = c( 2, 5 ),
#' sigma2_gamma = c( 2, 2 ), burnin = 50000, sampfreq = 50, nsample = 20000 )
#'
#' @export
ctlMCMCTree <- function( filename, mol = FALSE, seed = -1, seqfile, treefile, mcmcfile = "mcmc.txt",
                         outfile = "out.txt", ndata, seqtype = 0, usedata = 1, clock, RootAge,
                         TipDate, alpha, ncatG, cleandata = 0, BDparas, kappa_gamma, alpha_gamma,
                         rgene_gamma, sigma2_gamma, print = 2, burnin, sampfreq, nsample, model ){
  
  # Check parameters
  .checkCTL( filename = filename, mol = mol, seed = seed, seqfile = seqfile, treefile = treefile,
             mcmcfile = mcmcfile, outfile = outfile, seqtype = seqtype, usedata = usedata,
             ndata = ndata, clock = clock, TipDate = TipDate, RootAge = RootAge,
             alpha = alpha, ncatG = ncatG, cleandata = cleandata, BDparas = BDparas,
             rgene_gamma = rgene_gamma, sigma2_gamma = sigma2_gamma, print = print,
             burnin = burnin, sampfreq = sampfreq, nsample = nsample,
             kappa_gamma = kappa_gamma, alpha_gamma = alpha_gamma, model = model )
  
  # Write control file
  
  .writeCTL( filename = filename, mol = mol, seed = seed, seqfile = seqfile, treefile = treefile,
             mcmcfile = mcmcfile, outfile = outfile, seqtype = seqtype, usedata = usedata,
             ndata = ndata, clock = clock, TipDate = TipDate, RootAge = RootAge,
             alpha = alpha, ncatG = ncatG, cleandata = cleandata, BDparas = BDparas,
             rgene_gamma = rgene_gamma, sigma2_gamma = sigma2_gamma, print = print,
             burnin = burnin, sampfreq = sampfreq, nsample = nsample,
             kappa_gamma = kappa_gamma, alpha_gamma = alpha_gamma, model = model )
  
}

# ###############
# SUBFUNCTIONS ##
#################

.checkCTL <- function( filename, mol, seed, seqfile, treefile, mcmcfile, outfile, seqtype,
                       usedata, ndata, clock, TipDate, RootAge, alpha, ncatG,
                       cleandata, BDparas, rgene_gamma, sigma2_gamma, print, burnin,
                       sampfreq, nsample, kappa_gamma, alpha_gamma, model ){
  
  if( missing( filename ) ){
    stop( "\nPlease, provide a name for the output control file\n" )
  }else{
    if( class( filename ) != "character" ){
      stop( "\nPlease, provide a name for the output control file, class\"character\"\n" )
    }
  }
  # If a molecular partition is included, check kappa_gamma, alpha_gamma, alpha,
  # and ncatG are provided
  if( mol == TRUE ){
    if( missing( kappa_gamma ) ){
      stop( "\nPlease, provide a vector to parameter \"kappa_gamma\" with the alpha
            and beta values for the gamma prior for the substitution model
            parameter kappa (transition/transversion ratio)\n." )
    }else{
      if( length( kappa_gamma ) != 2 | class( kappa_gamma ) != "numeric" ){
        stop( "\nPlease, provide alpha and beta values (vector length 2) to parameter \"kappa_gamma\"
              for the gamma prior for the substitution model parameter kappa\n" )
      }
      }
    
    if( missing( alpha_gamma ) ){
      stop( "\nPlease, provide a vector to parameter \"alpha_gamma\" with the alpha
            and beta values for the gamma prior for the substitution model
            parameter alpha (gamma shape parameter for variable rates among sites)\n." )
      
    }else{
      if( length( alpha_gamma ) != 2 | class( alpha_gamma ) != "numeric" ){
        stop( "\nPlease, provide alpha and beta values (vector length 2) to paramter \"alpha_gamma\"
              for the gamma prior for the substitution model parameter alpha\n" )
      }
      }
    
    if( missing( alpha ) ){
      stop( "\nPlease, provide a value for the parameter \"alpha\" to set the gamma
            model of rate variation\n" )
    }else{
      if( class( alpha ) != "numeric" | length( alpha ) != 1 ){
        stop( "\nPlease, provide a numeric value for the parameter \"alpha\" to set
              the discrete-gamma model of rate variation\n" )
      }
      }
    
    if( missing( ncatG ) ){
      stop( "\nPlease, provide a numeric value for the number of categories for the
            discrete-gamma model of rate variation\n" )
    }else{
      if( as.integer( ncatG ) != ncatG | length( ncatG ) != 1 ){
        stop( "\nPlease, provide a numeric value for the paramter \"ncatG\" to set
              the number of categories for the discrete-gamma model of rate variation\n" )
      }
      }
    
    if( model != 0 & any( model == c( 1:10 ) ) ){
      stop( "\nPlease, provide a numeric value from 0 to 10 to set the substitution model\n" )
    }
    
  }
  
  if( seed != -1 & as.integer( seed ) != seed ){
    stop( "\nPlase, provide and odd or an even number for the seed.\n" )
  }
  
  if( mcmcfile != "mcmc.txt" & class( mcmcfile ) != "character" ){
    stop( "\nPlease, provide the path to the output \"mcmc.txt\" file,
          class \"character\".\n" )
  }
  
  if( outfile != "out.txt" & class( mcmcfile ) != "character" ){
    stop( "\nPlease, provide the path to the output \"out.txt\" file,
          class \"character\".\n" )
  }
  
  if( missing( seqfile ) ){
    stop( "\nPlease, provide the path to the alignment file\n.")
  }else{
    if( class( seqfile ) != "character" ){
      stop( "\nPlease, provide the path to the alignment file, class \"character\"\n.")
    }
  }
  
  if( missing( treefile ) ){
    stop( "\nPlease, provide the path to the tree file\n.")
  }else{
    if( class( treefile ) != "character" ){
      stop( "\nPlease, provide the path to the tree file, class \"character\"\n.")
    }
  }
  
  if( seqtype != 0 & any( seqtype != 1 & seqtype != 2 ) ){
    stop( "\nPlease, provide either 0, 1, or 2 to specify the type of molecular
          sequences in the alignment\n" )
  }
  
  if( usedata != 1 & any( usedata != 0 & usedata != 2 ) ){
    stop( "\nPlease, provide either 0, 1, or 2 to set the parameter \"usedata\"\n" )
  }
  
  if( missing( ndata ) ){
    stop( "\nPlease, provide the number of partitions in the alignment file\n." )
  }else{
    if( length( ndata ) != 1 | as.integer( ndata ) != ndata ){
      stop( "\nPlease, provide the number of partitions in the alignment file\n." )
    }
  }
  
  if( missing( clock ) ){
    stop( "\nPlease, provide a numerical value to set the clock model\n.")
  }else{
    if( length( clock ) != 1 | as.integer( clock ) != clock | any( clock!=1 & clock!=2 & clock!=3 ) ){
      stop( "\nPlease, provide a numerical value to set the clock model.
            Values allowed: 1, 2, or 3.\n.")
    }
    }
  
  if( missing( TipDate ) ){
    stop( "\nPlease, provide a numeric value to set the parameter \"TipDate\" with
          a time unit\n" )
  }else{
    if( as.integer( TipDate ) != TipDate | length( TipDate ) != 1 ){
      stop( "\nPlease, provide a numeric value to set the parameter \"TipDate\" with
            a time unit\n" )
    }
    }
  
  if( missing( RootAge ) ){
    stop( "\nPlease, provide a calibration for the root age.\n" )
  }else{
    if( class( RootAge ) != "character" ){
      stop( "\nPlease, provide a calibration for the root age, class = \"character\".\n" )
    }
  }
  
  if( cleandata != 0 & cleandata != 1 ){
    stop( "\nPlease, provide either 0 or 1 to specify parameter \"cleandata\"\n" )
  }
  
  if( missing( BDparas) ){
    stop( "\nPlease, provide a vector with the values of the
          birth-death-sequential-sampling (BDSS) process\n")
  }else{
    if( length( BDparas ) != 4 | class( BDparas ) != "numeric" ){
      stop( "\nPlease, provide a numeric vector with the values of the
            birth-death-sequential-sampling (BDSS) process. It needs
            a vector with four parameters: birth (lambda), death (mu),
            frequency to sample extant species (rho), and
            frequency to sample extinct species (psi).\n")
    }
    }
  
  if( missing( rgene_gamma ) ){
    stop( "\nPlease, provide a vector with the alpha and beta parameters
          for the Dirichlet-gamma prior for the mean substitution rate.\n" )
  }else{
    if( length( rgene_gamma ) != 2 | class( rgene_gamma ) != "numeric" ){
      stop( "\nPlease, provide a numeric vector with the alpha and beta parameters
            for the Dirichlet-gamma prior for the mean substitution rate.\n" )
    }
    }
  
  if( missing( sigma2_gamma ) ){
    stop( "\nPlease, provide a vector with the alpha and beta parameters
          for the Dirichlet-gamma prior for the rate drift parameter.\n" )
  }else{
    if( length( sigma2_gamma ) != 2 | class( sigma2_gamma ) != "numeric" ){
      stop( "\nPlease, provide a numeric vector with the alpha and beta parameters
            for the Dirichlet-gamma prior for the rate drift parameter.\n" )
    }
    }
  
  if( print != 0 & print != 1 & print != 2 ){
    stop( "\nPlease, provide either 0, 1, or 2 to specify parameter \"print\"\n" )
  }
  
  if( missing( burnin ) ){
    stop( "\nPlease, provide a numeric value for the burn-in\n" )
  }else{
    if( as.integer( burnin ) != burnin | length( burnin ) != 1 ){
      stop( "\nPlease, provide a numeric value for the burn-in\n" )
    }
  }
  
  if( missing( sampfreq ) ){
    stop( "\nPlease, provide a numeric value for the sampling frequency parameter\n" )
  }else{
    if( as.integer( sampfreq ) != sampfreq | length( sampfreq ) != 1 ){
      stop( "\nPlease, provide a numeric value for the sampling frequency parameter\n" )
    }
  }
  
  if( missing( nsample ) ){
    stop( "\nPlease, provide a numeric value for the amount of samples to collect\n" )
  }else{
    if( as.integer( nsample ) != nsample | length( nsample ) != 1 ){
      stop( "\nPlease, provide a numeric value for the amount of samples to collect\n" )
    }
  }
  
  }


.writeCTL <- function( filename, mol, seed, seqfile, treefile, mcmcfile, outfile, seqtype,
                       usedata, ndata, clock, TipDate, RootAge, alpha, ncatG,
                       cleandata, BDparas, rgene_gamma, sigma2_gamma, print, burnin,
                       sampfreq, nsample, kappa_gamma, alpha_gamma, model ){
  
  write( paste( "        seed   = ", seed, sep = "" ), file = filename )
  write( paste( "       seqfile = ", seqfile, sep = "" ), file = filename, append = T )
  write( paste( "      treefile = ", treefile, sep = "" ), file = filename, append = T )
  write( paste( "       outfile = ", outfile, sep = "" ), file = filename, append = T )
  write( paste( "      mcmcfile = ", mcmcfile, "\n", sep = "" ), file = filename, append = T )
  write( paste( "       seqtype = ", seqtype, "  * 0: nucleotides; 1:codons; 2:AAs", sep = "" ), file = filename, append = T )
  write( paste( "       usedata = ", usedata, "  * 0: no data; 1:seq like; 2:use in.BV; 3: out.BV", sep = "" ), file = filename, append = T )
  write( paste( "         ndata = ", ndata, "   * data partitions\n", sep = "" ), file = filename, append = T )
  write( paste( "         clock = ", clock, "   * 1: global clock; 2: independent rates; 3: correlated rates", sep = "" ), file = filename, append = T )
  write( paste( "       TipDate = 1 ", TipDate, "   * TipDate (1) & time unit", sep = "" ), file = filename, append = T )
  write( paste( "       RootAge = ", RootAge, "\n", sep = "" ), file = filename, append = T )
  if ( mol == TRUE ){
    write( paste( "         model = ", model, "\n", sep = "" ), file = filename, append = T )
  }else{
    write( "         model = 0    * 0:Brownian model (only model implemented with morpho data)", file = filename, append = T )
  }
  if( mol == TRUE ){
    write( paste( "         alpha = ", alpha, "   * alpha for gamma rates at sites", sep = ""), file = filename, append = T )
    write( paste( "         ncatG = ", ncatG, "   * No. categories in discrete gamma", sep = ""), file = filename, append = T )
  }
  write( paste( "\n     cleandata = ", cleandata, "    * remove sites with ambiguity data (1:yes, 0:no)?\n", sep = "" ), file = filename, append = T )
  write( paste( "       BDparas = ", paste( BDparas, collapse = " " ), "   * lambda, mu, rho, psi for birth-death-sequential-sampling model\n", sep = "" ), file = filename, append = T )
  if( mol == TRUE ){
    write( paste( "   kappa_gamma = ", paste( kappa_gamma, collapse = " " ), "   * gamma prior for kappa", sep = "" ), file = filename, append = T )
    write( paste( "   alpha_gamma = ", paste( alpha_gamma, collapse = " " ), "   * gamma prior for alpha\n", sep = "" ), file = filename, append = T )
  }
  write( paste( "   rgene_gamma = ", paste( rgene_gamma, collapse = " " ), "   * gamma prior for rate for genes", sep = "" ), file = filename, append = T )
  write( paste( "  sigma2_gamma = ", paste( sigma2_gamma, collapse = " " ), "   * gamma prior for sigma^2  (for clock=2)\n", sep = "" ), file = filename, append = T )
  write( "      finetune = 1: .1 .1 .1 .1 .1 .1  * auto (0 or 1) : times, musigma2, rates, mixing, paras, FossilErr\n", file = filename, append = T )
  write( paste( "         print = ", print, sep = "" ), file = filename, append = T )
  write( paste( "        burnin = ", burnin, sep = "" ), file = filename, append = T )
  write( paste( "      sampfreq = ", sampfreq, sep = "" ), file = filename, append = T )
  write( paste( "       nsample = ", nsample, sep = "" ), file = filename, append = T )
  
}



# ------------------------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------------------------- #

# --- CURRENTLY NOT USED --- #

# Estimate shrinkage correlation matrix from population matrix
#
# @description
# Estimate the shrinkage correlation matrix from a population that could
# have been either previously simulated with \code{\link{write.morpho}}
# or could refer to a population sample from a real data set.
#
# @param P Matrix, population matrix of size \code{s x n}, 's' specimens
# and 'n' characters.
#
# @param delta Numeric, shrinkage value used to generate the estimated
# shrinkage correlation matrix (see details).
#
# @return
#
# \item{P}{Matrix with the simulated population sample of
# \code{psample} specimens and \code{n} morphological continuous traits
# sampled from a normal distribution with mean = 0 and sd = sqrt(c)}
#
# @details
#
# Thre is a problem commonly faced with morphological data because, usually,
# the matrices with these data have more columns (more traits) than rows
# (specimens from which the traits were sampled). As the unbiased correlation
# matrices that can be calculated from these morphological matrices tend to
# be singular, they are not invertible. Consequently, they cannot be used
# to transform the data sets when accounting for trait correlation nor used to
# calculate the corresponding determinants, needed during the likelihood
# calculation in MCMCtree.
#
# The \code{delta} parameter is a shrinkage value used to find the estimate
# of the shrinkage correlation matrix, which is to be symmetric and invertible,
# thus it can be used to transform the morphological continuous data set
# and during the likelihood calculation in MCMCtree. This estimate can be found
# by using the following equation,
# \eqn{{R}^{*} = \delta I + (1-\delta){R'}}{R.sh = delta I + (1-delta) R.unb},
# as detailed in Schaffer & Strimmer, 2005.
#
# @seealso
# \code{\link{sim.morpho}}, \code{\link{write.morpho}}
#
# @author Sandra Alvarez-Carretero and Mario dos Reis
#
# @references
# \insertRef{Schafer2005}{morpho}
#
# @examples
#
# # Simulation setup: Estimate the shrinkage correlation matrix
# # of a sample population (20 specimens). The shrinkage value
# # is the one also used as a default.
#
#   calc.pop.cor( P = sim.population, delta = 0.01 )
#
# @export

# calc.pop.cor <- function( P, delta = 0.01 ){
#
#   # Check input parameters
#   .checkPopSamp( P, delta )
#
#   # Calculate the population variance
#   var.P <- diag( cov( P ) )
#
#   # Calculate estimated shrinkage correlation matrix
#   R.sh <- .CorrPopSamp( P = P, delta = delta )
#
#   # Test if R.sh is invertible. Otherwise,
#   # throw a warning saying that the value of delta
#   # should be changed
#   .testInvRsh( X = R.sh )
#
#   # Return a list with R.sh and the pop.var
#   return( list( R.sh = R.sh, var = var.P ) )
#
# }


## ################
##  SUBFUNCTIONS ##
## ################

# Check P is matrix
# .checkPopSamp <- function( P, delta ){
#
#   if( missing( P ) ){
#     stop( "\nPlease, provide an object of class matrix with the sampled population\n" )
#   }
#   if ( class( P ) != "matrix" ){
#     stop( "\nPlease, provide an object of class matrix with the sampled population\n" )
#   }
#
#   if ( class( delta ) != "numeric" ){
#     stop( "\nPlease, provide a numeric value for delta\n" )
#   }
#   if ( delta < 0 | delta > 1 ){
#     stop( "\n", "The shrinkage value delta used to generate the estimate of the
#           shrinkage correlation matrix should be positive and not larger than 1,
#           i.e. 1 >= delta > 0 ", "\n" )
#   }
#
# }


# Generate R.sh
# .CorrPopSamp <- function( P, delta ){
#
#   # Generate the identity matrix and the unbiased estimated
#   # correlation matrix
#   Id    <- diag( 1, dim( P )[2] )
#   R.unb <- cor( P )
#
#   # Get estimated shrinkage matrix
#   R.sh          <- delta*Id + ( 1 - delta )*R.unb
#   class( R.sh ) <- "matrix"
#
#   # Return R.sh
#   return( R.sh )
#
# }


# Test R.sh is invertible.
# Otherwise, throw a warning suggesting to change the delta value
# .testInvRsh <- function( X ){
#
#   test <- tryCatch( chol( X ), error = function( e ) e, warning = function( w ) w )
#   warn <- gsub( ".*Error..*", "Warning", test )
#   if( warn == "Warning" ){
#     warning( "The generated shrinkage correlation matrix is not invertible\nYou might like to rerun calc.pop.cor with\nanother value of delta")
#   }
#
# }



# ------------------------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------------------------- #



# ------------------------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------------------------- #


#####################
## EXTRA FUNCTIONS ##
#####################


#' Import various landmark files for different specimens at once
#'
#' @description
#' Import more than one csv file with landmark points and return an array
#' object with dimensions p x k x n, being \code{p} the number of
#' landmarks, \code{k} the dimension of the coordinates (2D or 3D), and
#' \code{n} the number of specimens (the number of files, as each file
#' contains the landmarks for one specimen).
#'
#' @param path Character, absolute path to the directory with the
#' csv files with the landmark points are.
#'
#' @param lmk.names Logical, TRUE if there is an extra column for
#' landmark names, FALSE otherwise (see details).
#'
#' @details
#' Note that all files need to be comma separated files (csv).
#'
#' If \code{lmk.names = TRUE}, the format expected for 3D landmarks is
#' the following:
#' \tabular{cccccccc}{
#' landmarks \tab x      \tab y      \tab z      \cr
#' lmk_1     \tab 0.143  \tab -0.028 \tab -0.030 \cr
#' lmk_2     \tab 0.128  \tab -0.024 \tab -0.035 \cr
#' ...       \tab ...    \tab ...
#' }
#'
#' Otherwise, if \code{lmk.names = TRUE}, then the format is:
#' \tabular{cccccccc}{
#' x      \tab y      \tab z      \cr
#' 0.143  \tab -0.028 \tab -0.030 \cr
#' 0.128  \tab -0.024 \tab -0.035 \cr
#' ...    \tab ...    \tab ...
#' }
#'
#' Note that you can always have 2D landmarks, so the format is the same
#' but the column with the \code{z} landmarks will not appear in the files.
#'
#' @author Sandra Alvarez-Carretero
#'
#' @export

lmk_imp <- function( path = NULL, lmk.names = FALSE ){

  # Set working directory

  setwd( path )

  # Get a list with the files

  f.names <- list.files( ".", pattern = "*.csv", full.names = TRUE )

  f.list  <- lapply( f.names, read.csv )

  if ( lmk.names == TRUE ){

    # Delete row names

    f.list <- lapply( f.list, function( x ) x[ -1 ] )

  }

  # Get specimen names and add them as names in the list

  sp.names <- rep( "", length( f.list ) )
  sp.names <- list.files( ".", pattern = ".csv" )
  sp.names <- gsub( ".csv", "", sp.names )

  names(f.list) <- sp.names

  # Get number of landmarks, specimens, and coordinates

  s     <- length( f.list )
  lmk   <- dim( f.list[[ 1 ]] )[ 1 ]
  coord <- dim( f.list[[ 1 ]] )[ 2 ]

  # Create empty 3D array ( p (lmk) x k (coord) x n (s) )

  arr             <- array( dim = c( lmk, coord, s ) )
  dimnames( arr ) <- list( paste( "lmk", seq( 1:( lmk ) ), sep="" ),
                           c( "x", "y", "z" ),
                           sp.names )

  # Fill in 3D array and return 3D array

  for ( i in 1:length( f.list ) ){

    arr[ , , i ] <- unlist( f.list[[ i ]] )

  }

  return( arr )

}







