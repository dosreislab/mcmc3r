#' Generate a phylip file for MCMCTree
#'
#' @description
#' Generate an alignment file in phylip format for MCMCTree.
#' The option "seqfile" in the control file used by MCMCTree
#' should read the path to the file output by this function.
#'
#' @param M Matrix, 's' rows (specimens) and 'n' morphological
#' continuous characters(see details).
#'
#' @param filename Character, name for the output file.
#'
#' @param c Numeric, vector of population variances (see details).
#' If not provided, c = 0 (no population variance).
#'
#' @param R Matrix, correlation matrix. Requires \code{method}
#' (see details). If not provided, R = I (no trait correlation).
#'
#' @param method (Optional) character, either \code{"eigen"} or
#' \code{"col"}, method used to decompose the inverse of the shrunk
#' correlation matrix. Requires \code{R} (see details).
#'
#' @param A (Optional) matrix, decomposed matrix. Requires \code{R}
#' but not \code{method} (see details).
#'
#' @param names (Optional) list, species name included in the morphological
#' alignment (see examples B and C).
#'
#' @param ages (Optional) list, ages of the species included in the
#' morpholical alignment (see example C).
#'
#' @details
#'
#' The matrix \code{X} has 's' rows, one for each specimen, and 'n' columns regarding
#' the characters. If the data set contains landmarks, they can be given in 2D or 3D.
#' For instance, if the landmarks are 3D, the first 3 columns will be the
#' coordinates x, y, and z for the first landmark, the next 3
#' columns for the second landmark, and so on.
#' \tabular{cccccccc}{
#'  specimens \tab lmk1.x \tab lmk1.y  \tab lmk1.z  \tab lmk2.x \tab lmk2.y \tab lmk2.z  \tab ... \cr
#'  Sp_1      \tab 0.143  \tab -0.028  \tab -0.044  \tab 0.129  \tab 0.028  \tab -0.043  \tab ... \cr
#'  Sp_2      \tab 0.128  \tab -0.024  \tab -0.028  \tab 0.124  \tab 0.027  \tab -0.025  \tab ... \cr
#'  ...       \tab ...    \tab ...     \tab ...     \tab ...    \tab ...    \tab ...     \tab ...
#' }
#' See \code{canids19x29.matrix} in the \code{data} directory for a
#' more detailed example and \code{morpho/data-raw/canids19x29.matrix.R} for
#' details on how to generate this object from the raw data.
#' If the data set contains a set of 'n' morphological continuous characters,
#' e.g. from a simulated data set, the file should look like
#' \tabular{cccccc}{
#'  specimens \tab trait.1 \tab trait.2  \tab trait.2  \tab trait.2 \tab ... \cr
#'  Sp_1      \tab 0.143   \tab -0.028   \tab -0.044   \tab 0.129   \tab ... \cr
#'  Sp_2      \tab 0.128   \tab -0.024   \tab -0.028   \tab 0.124   \tab ... \cr
#'  ...       \tab ...     \tab ...      \tab ...      \tab ...     \tab ...
#' }
#'
#' Note that if a list with the specimens names is not passed to the parameter \code{names},
#' the name for each species will be "Species_1", "Species_2", and so on.
#'
#' The object \code{c} can be of length 1, if all characters have the same variance, or
#' a vector of length 'n' with the variance of each of the characters.
#'
#' The object \code{R} has to be a symmetric and positive definite object of class matrix
#' (\code{class( R ) = "matrix"}).
#' The logarithm of the determinant of the correlation matrix is going to be printed
#' in the output file to later be used by MCMCTree during the likelihood calculation.
#'
#' If a correlation matrix \code{R} is provided, \code{write.morpho()} can use either
#' the \code{method = "chol"} or \code{method = "eigen"} to get a matrix \code{A}
#' such that \eqn{\mathrm{R^{-1}}=\mathrm{A^{T}}\mathrm{A}}{R^{-1} = t(A) * A}.
#' This matrix \eqn{\mathrm{A^{T}}}{t(A)} is later used to transform
#' the morphological data while to account for the correlation in this data set,
#' so that the transformed characters in \code{Z},
#' \eqn{\mathrm{Z}=\mathrm{M}\mathrm{A^{T}}}{Z = M * t(A)}, are independent.
#' Alternatively, this matrix \code{A} can also be provided by the user. If so,
#' this matrix will be used to transform the data and no decomposition will be
#' performed, saving computational time for big matrices.
#'
#' @seealso
#' \code{\link{matrix2array}}, \code{\link{array2matrix}}, \code{\link{sim.morpho}}
#'
#' @author Sandra Alvarez-Carretero and Mario dos Reis
#'
#' @examples
#' # A.1) Providing the morphological alignment (M) and
#' #      the name for the output file. This does not account for
#' #      correlation nor population variance.
#'
#'        write.morpho(  M = canids19x29.matrix, filename = "seqfile.aln" )
#'
#' # A.2) Providing the morphological alignment (M), the population
#' #      variance (c), and the name for the output file. Note that
#' #      c = 0.25 means that the population variance for all the traits
#' #      is c = 0.25, i.e. it will be considered as if
#' #      length( c ) = p characters, being all of them 0.25.
#'
#'        write.morpho(  M = canids19x29.matrix, c = 0.25,
#'                       filename = "seqfile.aln" )
#'
#' # A.3) Providing the morphological alignment (M), the population
#' #      variance (c), the estimate of the correlation matrix (R),
#' #      the method to decompose R ("chol" in this example),
#' #      and the name for the output file. Note that the R matrix needs
#' #      to be invertible, otherwise the data will not be able to be
#' #      transformed accountign for correlation.
#'
#'        write.morpho(  M = canids19x29.matrix, c = 0.25, R = R.shrunk,
#'                       method = "chol", filename = "seqfile.aln" )
#'
#' # A.4) Providing the morphological alignment (M), the population
#' #      variance (c), the estimate of the correlation matrix (R),
#' #      the A matrix to transform the data, and the name for the
#' #      output file. Note that as the A matrix is provided, the R matrix
#' #      will not be decomposed, hence the argument "method" is no needed.
#'
#'        write.morpho(  M = canids19x29.matrix, c = 0.25, R = R.shrunk,
#'                       A = A, filename = "seqfile.aln" )
#'
#' # B) Scenario A.3 but providing a list with the
#' #    names of the species
#'
#'      names <- list( sp1  = "Ael_sp.", sp2  = "Can_dir", sp3  = "Epi_hay", sp4  = "Hes_sp.",
#'                     sp5  = "Mes_cor", sp6  = "Tom_sp.", sp7  = "Enh_pah", sp8  = "Cuo_alp",
#'                     sp9  = "Spe_ven", sp10 = "Can_lup", sp11 = "Cer_tho", sp12 = "Oto_meg",
#'                     sp13 = "Vul_vul", sp14 = "Urs_ame", sp15 = "Ail_ful", sp16 = "Nan_bio",
#'                     sp17 = "Par_her", sp18 = "Tha_won", sp19 = "Smi_fat"
#'                   )
#'
#'      write.morpho( M = canids19x29.matrix, c = 0.25, R = R.shrunk,
#'                    A = A, filename = "seqfile.aln", names = names )
#'
#' # C) Scenario A.3 but providing a list with the names of
#' #    the specimens and a list with their corresponding ages. Please
#' #    keep the same order in both lists, so the first specimen in the
#' #    list name corresponds to the first age in the age list, and so on.
#'
#'      names <- list( sp1  = "Ael_sp.", sp2  = "Can_dir", sp3  = "Epi_hay", sp4  = "Hes_sp.",
#'                     sp5  = "Mes_cor", sp6  = "Tom_sp.", sp7  = "Enh_pah", sp8  = "Cuo_alp",
#'                     sp9  = "Spe_ven", sp10 = "Can_lup", sp11 = "Cer_tho", sp12 = "Oto_meg",
#'                     sp13 = "Vul_vul", sp14 = "Urs_ame", sp15 = "Ail_ful", sp16 = "Nan_bio",
#'                     sp17 = "Par_her", sp18 = "Tha_won", sp19 = "Smi_fat"
#'                   )
#'
#'      ages <- list( sp1  = 15.97, sp2  =  1.80, sp3  = 13.6,  sp4  = 39.74,
#'                    sp5  = 30.80, sp6  = 15.97, sp7  = 30.80, sp8  =  0,
#'                    sp9  =  0,    sp10 =  0,    sp11 =  0,    sp12 =  0,
#'                    sp13 =  0,    sp14 =  0,    sp15 =  0,    sp16 =  0,
#'                    sp17 =  0,    sp18 =  0,    sp19 =  0
#'                   )
#'
#'      write.morpho( M = canids19x29.matrix, c = 0.25, R = R.shrunk,
#'                    A = A, filename = "seqfile.aln",
#'                    names = names, ages = ages )
#'
#' @export

write.morpho <- function( M, filename, c = 0, R = diag( 1, dim( M )[2] ),
                          method, A, names, ages ) {

  # Check initial arguments
  .checkInArgs( X = M, filename = filename, c = c, R = R, method = method, A = A )

  # Create variables
  s     <- dim( M )[1]
  chars <- dim( M )[2]

  # Check optional arguments
  .checkOptArg( names, ages, M )

  # If names are provided...
  if ( ! missing( names ) ){
    spaces <- .availNames( names, ages )
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

  # Scale data accounting for population variance
  if ( c != 0 ){
    M        <- .scalePopVar( X = M, c = c, s = s, n = chars )
    scalevar <- 1
  }
  else{
    scalevar <- 0
  }

  # If the A matrix is not provided, calculate it and
  # then transform M
  if ( missing( A ) & all( R == diag( 1, dim( M )[2] ) ) == F ){
    # Match method with the argument provided
    method <- .checkMethod( method = method )
    # Calculate transformed matrix
    Z <- .CalcZ( X = M, R = R, method = method, n = chars,  s = s )
    # Get logarithm determinant of R
    lnd <- determinant( R )$modulus
  }
  # Otherwise, just use the A matrix provided to transform M
  else if ( ! missing( A ) & all( R == diag( 1, dim( M )[2] ) ) == F ){
    Z <- .CalcZwithA( X = M, A = A, n = chars, s = s )
    # Get logarithm determinant of R
    lnd <- determinant( R )$modulus
  }
  # Last, if there is no correlation at all, then
  else if ( all( R == diag( 1, dim( M )[2] ) ) == T ){
    Z     <- M
    lnd   <- 0
  }

  # Generate output file for MCMCTree
  .outFile( X = Z , names = names, chars = chars,
            scalevar = scalevar, lndetR = lnd, filename = filename )

}

## ################
##  SUBFUNCTIONS ##
## ################

# Check inp parameters
.checkInArgs <- function( X, filename, c, R, method, A ){

  # Check a name for the output file has been given
  if ( missing( filename ) ){
    stop( "\nPlease use the parameter \"filename\" to provide a name for the output file\n" )
  }

  # Check object X (morph.data) is provided and class is a matrix
  if ( missing( X ) ){
    stop( "\nPlease use an object of class \"matrix\" and dimensions \"s x n\" to convert
          into MCMCTree format\n" )
  }
  if( class( X ) != "matrix" ){
    stop( "\nPlease use an object of class \"matrix\" and dimensions \"s x n\" to convert
          into MCMCTree format\n" )
  }

  # Check population variance
  if( c < 0 | class( c ) != "numeric" ){
    stop( "\nThe parameter \"c\" needs a numeric value > 0 or a numeric vector
          of length \"n\" with values > 0\n" )
  }
  if( length( c ) != dim( X )[2] & length( c ) != 1 ){
    stop( "\nProvide a vector with variances of length equal to the number of
          morphological characters, \"length( c ) = n\", or to a unique numeric value,
          \"length( c ) = 1\", which will assume that all variances equal to this value\n" )
  }

  # Check that method goes always with R
  if ( ! missing( method) & all( R == diag( 1, dim( X )[2] ) ) == T ){
    stop( "\nPlease provide a correlation matrix so that the
             method you have selected can be used\n" )
  }
  # Check correlation matrix is invertible and that method
  # has been providd
  if ( dim( R )[1] != dim( X )[2] & dim( R )[2] != dim( X )[2] ){
    stop( "\nPlease use a correlation matrix of size 'n x n', where
             'n' is the amount of characters present in your
             'M' matrix\n" )
  }
  if ( all( R == diag( 1, dim( X )[2] ) ) == F ){
    if ( missing( A ) & missing( method ) ){
      stop( "\nPlease select a method to decompose the shrunk correlation matrix,
          either method = \"chol\" or method = \"eigen\" \n" )
    }
    .checkCorrMat( R = R, n = dim( X )[2] )
  }

  # If A provided, check it is class matrix
  if ( ! missing( A ) ){

    if ( missing( R ) ){
      stop( "\nIf you provide A, you need to specify also R\n")
    }
    if ( class( A ) != "matrix" ){
      stop( "\nObject \"A\" needs to be of class \"matrix\"\n" )
    }
    # Check that the A matrix provided is correct given their R
    A2 <- .CalcCholesky( R = R, n = dim( X )[2] )
    if ( A != A2 ){
      stop( "\nYour \"A\" matrix does not seem to be the correct
             one according the \"R\" matrix you have provied.
             If you are unsure of properly calculating A, then
             use write.morpho without this argument so it is
             internally calculated\n" )
    }

  }

}


# Check if optional parameters (names and ages) added and, if so,
# if class is ok
.checkOptArg <- function( names, ages, M ){

  # Check names is class "list" and length(names) = dim(M)[1]
  if( ! missing( names ) ){
    if ( class( names ) != "list" & class( names ) != "character" ){
      stop( "\nYou need to provide an object of class list or character with the species included in the alignment file.
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
      stop( "\nYou need to provide an object of class list or numeric with the ages of the species included in the alignment file.
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

  ages       <- unlist( ages )
  names      <- paste( names, - ages + max( ages ) + 0.01, sep = "^" ) # 0.01 = ct for MCMCTree

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
    spaces <- .hasSpaces( names )
  }

  # If ages are provided ...
  else{
    # Put together names and transformed ages and get spaces
    names <- .hasAges( ages )
    # Get spaces to justify text
    spaces <- .hasSpaces( names )
  }

  # Return names and spaces
  return( list( names = names, spaces = spaces ) )

}


# Get spaces when names are provided.
# This wraps around .hasAges and .hasSpaces
.availNames <- function( names, ages ){

  # Get species names
  names   <- unlist( names )

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
  return( spaces )

}


# Scale matrix "X" with population variance
# Note that c has to be a vector with all the pop variances,
# so if length(c) == 1 it is first corrected to be of
# length "n" with the same "c" pop variance
.scalePopVar <- function( X, c, s, n ){

  M.s <- matrix( rep( 0, s*n ), nrow = s, ncol = n )
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
    stop( "\nPlease select a method to decompose the shrunk correlation matrix,
          either method = \"chol\" or method = \"eigen\" \n" )
  }

  # Match arguments and return the method
  methods <- c( "chol", "eigen" )
  method <- match.arg( arg = method, choices = methods )
  return( method )

}

# Decompose the correlation matrix using the Cholesky decomposition
.CalcCholesky <- function( R, n ){

  # Cholesky decomposition:
  # R = L %*% U = L %*% t( L )    = t( U ) %*% U
  # R^-1 = t( L^-1 ) %*% L^-1     = U^-1 %*% t( U^-1 )
  # R^-1 = t( A ) %*% A           = A %*% t( A )

  U <- matrix( rep( 0 ), ncol = n, nrow = n )
  U <- chol( R )
  A <- backsolve( U, diag( dim( U )[1] ) )

  # Return inverse of upper triangular matrix ( U^-1 = A )
  return( A )

}


# Decompose the correlation matrix using the eigen decomposition
.CalcEigen <- function( R, n ){

  # "eigen" returns a list with the eigenvectors
  # and the eigenvalues
  #
  # R^-1 = t( A ) %*% A
  #
  # t( A ) = V %*% D
  # t( A ) = eigen( R^-1 )$vectors %*% diag( sqrt( eigen( R^-1 )$values ) )

  Rinv <- tA <- matrix( rep( 0 ), ncol = n, nrow = n )
  #Rinv <- solve( R )
  Rinv <- chol2inv( chol( R ) )
  tA   <- eigen( Rinv )$vectors %*% diag( sqrt( eigen( Rinv )$values ) )

  # Return t(A)
  return( tA )

}


# Calculate transformed matrix when A is not provided
# and R needs to be decomposed
.CalcZ <- function( X, R, method, n, s ){

  # Create empty matrices
  Z <- matrix( rep( 0, s*n ), nrow = s, ncol = n )
  # Match argument
  if ( method == "chol" ){
    A <- .CalcCholesky( R = R, n = n )
  }

  else if ( method == "eigen" ){
    A <- .CalcEigen( R = R, n = n )
  }

  # Transform data
  Z <- X %*% A

  # Return Z
  return( Z )

}


# Calculate transformed matrix when A is provided
.CalcZwithA <- function( X, A, n, s ){

  Z <- matrix( rep( 0, s*n ), nrow = s, ncol = n )
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


#' Convert a matrix into a 3D array
#'
#' @description
#' Convert a matrix with landmark points into a 3D array object.
#'
#' @param X Matrix (s x n), 'n' landmark points for 's' specimens
#' (see details).
#'
#' @param coords Integer, 2 or 3 for 2D or 3D landmarks, respectively.
#'
#' @details
#'
#' The matrix has format s x n, with 's' rows, one for specimen, and 'n'
#' columns, one for each coordinate of the landmarks.
#' Each landmark can be given in 2D or 3D. For instance,
#' if the landmarks are 3D, the first 3 columns will be the
#' coordinates x, y, and z for the first landmark, the next 3
#' columns for the second landmark, and so on:
#' \tabular{cccccccc}{
#'  specimens \tab lmk1.x \tab lmk1.y  \tab lmk1.z  \tab lmk2.x \tab lmk2.y \tab lmk2.z  \tab ... \cr
#'  Sp_1      \tab 0.143  \tab -0.028  \tab -0.044  \tab 0.129  \tab 0.028  \tab -0.043  \tab ... \cr
#'  Sp_2      \tab 0.128  \tab -0.024  \tab -0.028  \tab 0.124  \tab 0.027  \tab -0.025  \tab ... \cr
#'  ...       \tab ...    \tab ...     \tab ...     \tab ...    \tab ...    \tab ...     \tab ...
#' }
#'
#' See \code{canids19x29.matrix} in the \code{data} directory for an example
#' and \code{morpho/data-raw/canids19x29.matrix.R} to see how this object
#' was generated from the raw file.
#'
#' @return
#'
#' An object of class array with format k x n x s, where 'k' is the number of
#' landmarks, 'n' the number of coordinates, and 's' the number of specimens.
#'
#' Note that if the matrix provided does not have rownames, the specimens in the
#' returned array (dimension 's') will be labelled as '1', '2', and so on.
#' See \code{canids19x29.array} in the \code{data} directory for an example of
#' the format of the object that is returned and
#' \code{/morpho/data-raw/canids19x29.array.R} for the description of how to
#' obtain this object from the raw data.
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
    stop( "\nPlease provide an object of class \"matrix\" with \"n x p\" dimensions, 'n' specimens and 'p' characters\n" )
  }
  if ( class( X ) != "matrix" ){
    stop( "\nPlease use an object of class \"matrix\"\n" )
  }

}


# Check coords
.checkCoords <- function( coords ){

  # Check coords are provided and are numeric
  if ( missing( coords ) ){
    stop( "\nThe parameter \"coords\" needs a numeric value. Please use \"coords = 2\" if the coordinates]
          are 2D or \"coords = 3\" if 3D\n" )
  }
  if ( class( coords ) != "numeric" ){
    stop( "\nThe parameter \"coords\" needs a numeric value. Please use \"coords = 2\" if the coordinates]
          are 2D or \"coords = 3\" if 3D\n" )
  }

}

# Wraper around .coords3D.mat2arr and .coords2D.mat2arr to create the
# array with landmarks
.genArr <- function( X, coords, chars, s, list.coords ){

  # Create an empty array to store the coordinates in the format
  # p x k x n (chars x coord x s)
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

  # Fill in array k x n x s
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

  # Fill in array k x n x s
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





#' Convert a 3D array into a matrix
#'
#' @description
#' Convert a 3D array with landmark points into an object of class matrix.
#'
#' @param X 3D array, 'k' landmark points, 'n' coordinates, and
#' 's' specimens.
#'
#' @param coords Integer, 2 or 3, for 2D or 3D landmarks, respectively.
#'
#' @details
#'
#' The object \code{X}, class array, has format k x n x s, where 'k' is
#' the number of landmarks, 'n' the number of coordinates, and 's' the number
#' of specimens. See \code{canids19x29.array} in the \code{data} directory for an
#' example of the format of a 3D array and \code{morpho/data-raw/canids19x29.array.R}
#' for the details about how to generate this object from the raw data.
#'
#' @return
#'
#' An object of class matrix, with 's' rows, one for specimen, and 'n' columns, one
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
#' See \code{canids19x29.matrix} in the \code{data} directory for a
#' more detailed example of the format of the object that is returned and
#' \code{morpho/data-raw/canids19x29.matrix.R} for the explanation about
#' how to generate this object from the raw data.
#'
#' Note that if the 's' dimension (specimens) of the array provided does
#' not have names, the specimens in the returned matrix will be
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
  .checkCoords <- function( coords )
  .checkArr    <- function( X, coords )

  # Create variables
  chars  <- NULL
  chars  <- dim( X )[1]
  coords <- coords
  s      <- dim( X )[3]

  # Call function to create array
  faln <- .genMat( X = X, coords = coords,chars = chars, s = s )

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
    stop( "\nPlease provide an object of class \"array\" with \"p x k x n\" dimensions,
          'p' landmarks, 'k' coordinates, and 'n' specimens\n" )
  }
  if ( class( X ) != "array" ){
    stop( "\nPlease use an object of class \"array\"\n" )
  }

  # Check argument coords is alright compared to array X
  if ( coords != dim( X )[2] ){
    stop( "\nThe number of coords does not match the number of
             coordinates in the array you have provided")
  }

}

# Wraper around .coords3D.arr2mat and .coords2D.arr2mat to create the
# array with landmarks
.genMat <- function( X, coords, chars, s ){

  # Create an empty matrix to store the coordinates in the format
  # s x n (specimens x characters )
  faln <- matrix( 0, nrow = s, ncol = chars * coords )

  # Select x, y, z positions
  if ( coords == 3 ){
    faln <- .coords3D.arr2mat( X = X , faln = faln, chars = chars, coords = coords, s = s )
  }

  else if ( coords == 2 ){
    faln <- .coords2D.arr2mat( X = X , faln = faln, chars = chars, s = s )
  }

  # Return array with positions x, y, and z
  return( faln )

}

# Convert 3D landmarks into matrix
.coords3D.arr2mat <- function( X, faln, chars, coords, s ){

  # Select x, y, z positions
  xi <- seq( from = 1, to = chars * coords, by = 3 )
  yi <- seq( from = 2, to = chars * coords, by = 3 )
  zi <- seq( from = 3, to = chars * coords, by = 3 )

  # Fill in matrix and put rownames
  for ( i in 1:s ) {
    faln[i,xi] <- X[,1,i]
    faln[i,yi] <- X[,2,i]
    faln[i,zi] <- X[,3,i]
  }
  rownames( faln ) <- dimnames( X )[[ 3 ]]

  # Return matrix
  return( faln )

}

# Convert 2D landmarks into matrix
.coords2D.arr2mat <- function( X, faln, chars, s ){

  # Select x and y positions
  xi <- seq( from = 1, to = chars, by = 2 )
  yi <- seq( from = 2, to = chars, by = 2 )

  # Fill in matrix and put rownames
  for ( i in 1:s ) {
    faln[i,xi] <- X[,1,i]
    faln[i,yi] <- X[,2,i]
  }
  rownames( faln ) <- dimnames( X )[[ 3 ]]

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
#' and later allowing to account for population variance and trait correlation.
#'
#' @param tree Phylo, object with a phylogenetic tree
#' (see \code{\link[ape]{rTraitCont}}).
#'
#' @param n Numeric, number of morphological traits to be simulated.
#'
#' @param c (Optional) numeric, vector with population variances to
#' add to the simulated morphological traits (see details).
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
#' These parameters are the ones the user can pass through
#' \code{sim.morpho()}.
#' The default values that \code{sim.morpho()} uses are
#' \code{model = "BM"}, \code{sigma = 1}, \code{ancestor = F},
#' and \code{root.value = 0}. For this kind of simulation,
#' \code{sim.morpho} allows only \code{ancestor = F}, so please
#' do not change this parameter.
#' In the \code{\link[ape]{rTraitCont}} package, the parameter
#' \code{model} can be \code{model = BM}, \code{model = OU}, or
#' a function \code{model = FUN} provided by the user. Currently,
#' \code{sim.morpho()} supports only the first two.
#'
#' The parameter \code{c} is the population variance and used to simulate
#' the noise matrix, which follows a normal distribution \code{x ~ N(0,c)}.
#' If the population variance is assumed to be the same for all traits,
#' then the length of \code{c} is 1 and equals to the value of this variance.
#' If it differs from trait to trait, then a vector of length \code{n} has to
#' be provided specifying the variance for each of the traits.
#'
#' The simulated noise is later added to the morphological data previously
#' generated, so we obtain the noisy matrix.
#' If a correlation matrix, \code{R}, is provided, then it is added to the
#' noisy matrix. Note that the correlation matrix needs to be of class "matrix"
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
#' #    Population variance and correlation are not considered,
#' #    i.e. c = 0 and R not provided.
#'
#'      morpho::sim.morpho( tree = sim.tree, n = 100 )
#'
#'
#' # B) Simulation setup: Simulate a morphological alignment
#' #    with n = 100 continuous characters for a phylogeny
#' #    defined in object 'tree', but with different parameters
#' #    than the default ones in 'sim.morpho' to run
#' #    'rTraitCont'. Population variance and correlation are not
#' #    considered, i.e. c = 0 and R not provided.
#'
#'      morpho::sim.morpho( tree  = sim.tree, n = 100,
#'                          model = "OU", sigma = 0.2, alpha = 2 )
#'
#'
#' # C) Simulation setup: Simulate a morphological alignment
#' #    with n = 100 continuous characters for a phylogeny
#' #    defined in object 'tree', with the default parameters in
#' #    'sim.morpho' to run 'rTraitCont'.
#' #    Population variance is low, c = 0.25, but correlation is not
#' #    considered, i.e. R not provided.
#'
#'      morpho::sim.morpho( tree = sim.tree, n = 100, c = 0.25 )
#'
#'
#' # D) Simulation setup: Simulate a morphological alignment
#' #    with n = 100 continuous characters for a phylogeny
#' #    defined in object 'tree', with the default parameters in
#' #    'sim.morpho' to run 'rTraitCont'.
#' #    Population variance is low, c = 0.25, and a correlation
#' #    matrix is provided.
#'
#'      morpho::sim.morpho( tree = sim.tree, n = 100, c = 0.25, R = sim.R )
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

  # Sample num.species*traits samples from a normal distribution
  # with mu = 0 and sd = c
  N <- .simNoise( s = s, n = n, c = c )

  # Generate noisy matrix
  M.n <- M + N

  # Add trait correlation if R matrix is provided
  if ( ! missing( R ) ){
    .checkCorrMat2( R = R, n = n )
    M.n <- .addCorrMat( X = M.n, R = R )
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
  if( class( n ) != "numeric" ){
    stop( "\nPlease, provide a numerical value with the amount of
          morphological continuous characters to be simulated for
          the population\n" )
  }
  if( ! n > 0 ){
    stop( "\nPlease, provide a numerical value with the amount of
          morphological continuous characters to be simulated for
          the population ( n > 0 )\n" )
  }
  if ( as.integer( n ) != n ){
    stop( "\nPlease provide an integer value for the amount of characters to
          be simulated\n")
  }

  if ( c < 0 ){
    stop( "\n", "The population variance should be positive, i.e. c > 0 ", "\n" )
  }
  if( length( c ) != n & length( c ) != 1 ){
    stop( "\nProvide a vector of variances of length equal to the number of
          morphological characters, n, or to a unique numeric value,
          which will assume that all variances equal to this value\n" )
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

  # Create matrix s x n (species x characters)
  # and then fill it in
  M <- matrix( rep( 0, s*n ), nrow = s, ncol = n )
  M <- replicate( n,
                  ape::rTraitCont( phy      = tree,          model      = pars$model,
                                   sigma    = pars$sigma,    root.value = pars$root,
                                   ancestor = pars$ancestor, alpha      = pars$alpha,
                                   theta    = pars$theta ) )

  #Return matrix M
  return( M )

}


# Checking correlation matrix, if provided, is class = "matrix",
# positive definite, and symmetric
.checkCorrMat <- function( R, n ){

  if ( class( R ) != "matrix" ){
    stop( "\nPlease use a correlation matrix of class \"matrix\"\n" )
  }

  if ( dim( R )[1] != n & dim( R )[2] != n ){
    stop( "\nPlease use a correlation matrix of size \"n x n\"\n" )
  }

  if ( isSymmetric( R ) != TRUE ){
    stop( "\nPlease provide a symmetric correlation matrix of class \"matrix\"\n" )
  }

  # Check R is a positive definite matrix -- This might take a while to check
  # depending on R size
  eigenvals     <- eigen( R )$values
  neg.eigenvals <- which ( eigenvals < 0 )
  if ( length( neg.eigenvals ) != 0 ){
    stop( "\nPlease provide a positive definite correlation matrix of class \"matrix\"\n" )
  }

}


# Simulate noise matrix
# Sample num.species*traits samples from a normal distribution
# distrib with mu = 0 and sd = c
.simNoise <- function( s, n, c ){

  N <- matrix( rep( 0, s*n ), nrow = s, ncol = n )
  N <- t( replicate( s, rnorm( n, mean = 0, sd = sqrt( c ) ) ) )

  # Return noise matrix
  return( N )

}


# Checking correlation matrix, if provided, is class = "matrix",
# and symmetric
.checkCorrMat2 <- function( R, n ){

  if ( class( R ) != "matrix" ){
    stop( "\nPlease use a correlation matrix of class \"matrix\"\n" )
  }

  if ( dim( R )[1] != n & dim( R )[2] != n ){
    stop( "\nPlease use a correlation matrix of size \"n x n\"\n" )
  }

  if ( isSymmetric( R ) != TRUE ){
    stop( "\nPlease provide a symmetric correlation matrix of class \"matrix\"\n" )
  }


}


# Add correlation to noisy matrix
.addCorrMat <- function( X, R ){

  # Check the class of R is "matrix"
  # [ALREADY CHECKED IN MAIN FUNC]
  #.checkCorrMat2( R )

  # Add correlation to noisy matrix M.n as an independent variable
  M.n.R <- t( chol( R ) ) %*% t( X )
  M.n.R <- t( M.n.R )

  # Return M.n.R
  return( M.n.R )

}


# ------------------------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------------------------- #

#' Simulate a population matrix
#'
#' @description
#' Simulate a population sample and return in a matrix of
#' size \code{s x n}, 's' specimens and 'n' characters.
#'
#' @param psample Numeric, number of specimens the simulated
#' population sample should include.
#'
#' @param n Numeric, number of morphological traits to be simulated.
#'
#' @param c Numeric, vector with the population variances with
#' which the population should be simulated (see details).
#'
#' @param R (Optional) matrix, correlation matrix.
#' (see details).
#'
#' @return
#'
#' \item{P}{Matrix with the simulated population sample of
#' \code{psample} specimens and \code{n} morphological continuous traits
#' sampled from a normal distribution with mean = 0 and sd = sqrt(c)}
#
#' @details
#'
#' The parameter \code{c} is the population variance and it is used to sample
#' \code{n} characters for each of the \code{psample} specimens from a
#' normal distribution \code{x ~ N(0,c)}.
#' If the population variance is assumed to be the same for all traits,
#' then the length of \code{c} is 1 and equals to the value of this variance.
#' If it differs from trait to trait, then a vector of length \code{n} has to
#' be provided specifying the variance for each of the traits.
#'
#' If a correlation matrix, \code{R}, is provided, then it is added to the
#' noisy matrix. Note that the correlation matrix needs to be of class "matrix"
#' and symmetric.
#'
#' @seealso
#' \code{\link{sim.morpho}}, \code{\link{write.morpho}}
#'
#' @author Sandra Alvarez-Carretero and Mario dos Reis
#'
#' @examples
#'
#' # A) Simulation setup: Simulate a population with a sample of
#' #    psample = 20 specimens, with n = 100 characters, and
#' #    a low population variance, c = 0.25.
#'
#'      morpho::sim.pop( psample = 20, n = 100, c = 0.25 )
#'
#' # B) Simulation setup: Simulate a population with a sample of
#' #    psample = 20 specimens, with n = 100 characters,
#' #    a low population variance, c = 0.25, and a low correlation
#' #    rho = 0.50 used to generate a correlation matrix that follows
#' #    the constant correlation model (all non-diagonal values
#' #    equal to rho). This matrix already exists as an object
#' #    within this package, sim.R, so take a look to
#' #    morpho/data-raw/sim.R.R to see how it was generated
#'
#'      morpho::sim.pop( psample = 20, n = 100, c = 0.25, R = sim.R )
#'
#' @export

sim.pop <- function( psample, n, c, R ){

  # Check input parameters
  .checkInSimPop( psample = psample, n = n, c = c )

  # Sample from a normal distribution the morph. cont. traits
  P <- matrix( rep( 0, psample*n ), nrow = psample, ncol = n )
  P <- t( replicate( psample, rnorm( n, mean = 0, sd = sqrt( c ) ) ) )

  # Add correlation if provided
  if ( ! missing( R ) ){
    .checkCorrMat2( R = R, n = n )
    P <- .addCorrMat( X = P, R = R )
  }
  #var <- diag( cov( P ) )

  # Get names for the specimens as "Species_1", "Species_2", and so on
  vars          <- .notNames( s = psample )
  rownames( P ) <- vars$names

  # Return P
  return( P )

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
  if( class( n ) != "numeric" ){
    stop( "\nPlease, provide a numerical value with the amount of
          morphological continuous characters to be simulated for
          the population\n" )
  }
  if( ! n > 0 ){
    stop( "\nPlease, provide a numerical value with the amount of
          morphological continuous characters to be simulated for
          the population ( n > 0 )\n" )
  }
  if ( as.integer( n ) != n ){
    stop( "\nPlease provide an integer value for the amount of characters to
          be simulated\n")
  }

  if ( missing( c ) ){
    stop( "\nPlease, a vector of variances of length equal to the number of
          morphological characters, n, or to a unique numeric value,
          which will assume that all variances equal to this value\n" )
  }
  if ( c < 0 ){
    stop( "\n", "The population variance should be positive, i.e. c > 0 ", "\n" )
  }
  if( length( c ) != n & length( c ) != 1 ){
    stop( "\nProvide a vector of variances of length equal to the number of
          morphological characters, n, or to a unique numeric value,
          which will assume that all variances equal to this value\n" )
  }

}


# ------------------------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------------------------- #

#' Estimate shrunk correlation matrix from population matrix
#'
#' @description
#' Estimate the shrunk correlation matrix from a population that could
#' have been either previously simulated with \code{\link{write.morpho}}
#' or could refer to a population sample from a real data set.
#'
#' @param P Matrix, population matrix of size \code{s x n}, 's' specimens
#' and 'n' characters.
#'
#' @param delta Numeric, shrinkage value used to generate the estimated
#' shrunk correlation matrix (see details).
#'
#' @return
#'
#' \item{P}{Matrix with the simulated population sample of
#' \code{psample} specimens and \code{n} morphological continuous traits
#' sampled from a normal distribution with mean = 0 and sd = sqrt(c)}
#
#' @details
#'
#' Thre is a problem commonly faced with morphological data because, usually,
#' the matrices with these data have more columns (more traits) than rows
#' (specimens from which the traits were sampled). As the unbiased correlation
#' matrices that can be calculated from these morphological matrices tend to
#' be singular, they are not invertible. Consequently, they cannot be used
#' to transform the data sets when accounting for trait correlation nor used to
#' calculate the corresponding determinants, needed during the likelihood
#' calculation in MCMCTree.
#'
#' The \code{delta} parameter is a shrinkage value used to find the estimate
#' of the shrunk correlation matrix, which is to be symmetric and invertible,
#' thus it can be used to transform the morphological continuous data set
#' and during the likelihood calculation in MCMCTree. This estimate can be found
#' by using the following equation,
#' \eqn{{R}^{*} = \delta I + (1-\delta){R'}}{R.sh = delta I + (1-delta)R.unb},
#' as detailed in Schaffer & Strimmer, 2005.
#'
#' @seealso
#' \code{\link{sim.morpho}}, \code{\link{write.morpho}}
#'
#' @author Sandra Alvarez-Carretero and Mario dos Reis
#'
#' @references
#' \insertRef{Schafer2005}{morpho}
#'
#' @examples
#'
#' # Simulation setup: Estimate the shrunk correlation matrix
#' # of a sample population (20 specimens). The shrinkage value
#' # is the one also used as a default.
#'
#'   morpho::calc.pop.cor( P = sim.population, delta = 0.01 )
#'
#' @export

calc.pop.cor <- function( P, delta = 0.01 ){

  # Check input parameters
  .checkPopSamp( P, delta )

  # Calculate the population variance
  var.P <- diag( cov( P ) )

  # Calculate estimated shrunk correlation matrix
  R.shrunk <- .CorrPopSamp( P = P, delta = delta )

  # Test if R.shrunk is invertible. Otherwise,
  # throw a warning saying that the value of delta
  # should be changed
  .testInvRsh( X = R.shrunk )

  # Return a list with R.shrunk and the pop.var
  return( list( R.shrunk = R.shrunk, var = var.P ) )

}


## ################
##  SUBFUNCTIONS ##
## ################

# Check P is matrix
.checkPopSamp <- function( P, delta ){

  if( missing( P ) ){
    stop( "\nPlease, provide an object of class matrix with the sampled population\n" )
  }
  if ( class( P ) != "matrix" ){
    stop( "\nPlease, provide an object of class matrix with the sampled population\n" )
  }

  if ( class( delta ) != "numeric" ){
    stop( "\nPlease, provide a numeric value for delta\n" )
  }
  if ( delta < 0 | delta > 1 ){
    stop( "\n", "The shrinkage value delta used to generate the estimate of the
          shrunk correlation matrix should be positive and not larger than 1,
          i.e. 1 >= delta > 0 ", "\n" )
  }

}

# Generate R.shrunk
.CorrPopSamp <- function( P, delta ){

  # Generate the identity matrix and the unbiased estimated
  # correlation matrix
  Id    <- diag( 1, dim( P )[2] )
  R.unb <- cor( P )

  # Get estimated shrunk matrix
  R.shrunk          <- delta*Id + ( 1 - delta )*R.unb
  class( R.shrunk ) <- "matrix"

  # Return R.shrunk
  return( R.shrunk )

}


# Test R.shrunk is invertible.
# Otherwise, throw a warning suggesting to change the delta value
.testInvRsh <- function( X ){

  test <- tryCatch( chol( X ), error = function( e ) e, warning = function( w ) w )
  warn <- gsub( ".*Error..*", "Warning", test )
  if( warn == "Warning" ){
    warning( "The generated R.shrunk is not invertible\nYou might like to rerun calc.pop.cor with\nanother value of delta")
  }

}



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







