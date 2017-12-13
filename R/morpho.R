#' Generate a phylip file for MCMCTree
#'
#' @description
#' Generate an alignment file in phylip format for MCMCTree.
#' The option "seqfile" in the control file used by MCMCTree
#' should read the path to the file output by this function.
#'
#' @param filename Character, name for the output file.
#'
#' @param proc 3D array (p x k x n) or matrix (n x p), landmark coordinates
#' after a Procrustes analysis (see details).
#'
#' @param coords Integer, only needed if \code{class(proc) == "array"}
#' (see details).
#'
#' @param names (Optional) list, species name included in the morphological
#' alignment (see examples B and C).
#'
#' @param ages (Optional) list, ages of the species included in the
#' morpholical alignment (see example C).
#'
#' @param popvar (Optional) vector of integers, population variance (see details).
#'
#' @param R (Optional) matrix, correlation matrix. Requires \code{popvar} and
#' \code{method}.
#'
#' @param method (Optional) character, method to decompose the inverse of the
#' correlation matrix R, \code{eigen} or \code{chol} (see details).
#' Requires \code{popvar} and \code{R}.
#'
#' @param ... Further arguments passed to \code{write_morpho()} (see details).
#'
#' @details
#'
#' If \code{class(proc) == "array"}, the array has format p x k x n, where
#' 'p' is the number of landmarks, 'k' the number of coordinates,
#' and 'n' the number of specimens. See \code{canids19x29.array} in the
#' \code{data} directory for an example.
#' You also need to provide \code{coords = 2} (2D landmark points) or
#' \code{coords = 3} (3D landmark points) if the \code{proc} object is an array.
#'
#' If \code{class(proc) == "matrix"}, the matrix has 'n' rows, one for specimen,
#' and 'p' columns.Each landmark can be given in 2D or 3D. For instance,
#' if the landmarks are 3D, the first 3 columns will be the
#' coordinates x, y, and z for the first landmark, the next 3
#' columns for the second landmark, and so on.
#' \tabular{cccccccc}{
#'  specimens \tab lmk1.x \tab lmk1.y  \tab lmk1.z  \tab lmk2.x \tab lmk2.y \tab lmk2.z  \tab ... \cr
#'  Sp_1      \tab 0.143  \tab -0.028  \tab -0.044  \tab 0.129  \tab 0.028  \tab -0.043  \tab ... \cr
#'  Sp_2      \tab 0.128  \tab -0.024  \tab -0.028  \tab 0.124  \tab 0.027  \tab -0.025  \tab ... \cr
#'  ...       \tab ...    \tab ...     \tab ...     \tab ...    \tab ...    \tab ...     \tab ...
#' }
#' See \code{canids19x29.matrix} in the \code{data} directory for a
#' more detailed example.
#'
#' If \code{names} is not provided, the name for each species will be
#' "Species_1", "Species_2", and so on.
#'
#' The parameter \code{popvar} accounts for the population variance. First of all,
#' a PA regarding the landmark points collected from different specimens
#' belonging to different species should be carried out, which will result into
#' a matrix of superimposed coordinates, say \code{M}. If there are many specimens of
#' one of the species included in \code{M}, then another PA is carried out with
#' the landmark points collected from these specimens belonging to the same species.
#' The resulting matrix, say \code{V}, is used to obtain the population variance.
#' The variance-covariance matrix of \code{V}, say \code{VCV}, is first calculated
#' and, afterwards, the diagonal of \code{VCV}, the vector of variances, is used
#' to scale the matrix \code{M} according to the standard deviation.
#' Therefore, \code{popvar = diag(VCV)}. If this parameter is not provided, then
#' the population variance is assumed to be 0.
#'
#' If \code{method = "eigen"}, the inverse of the correlation matrix provided,
#'  \code{R}, will be decomposed following the spectral eigendecomposition
#' (see \code{eigen}).
#' If \code{method = "chol"}, the inverse of \code{R} will be decomposed
#' following the Cholesky decomposition (see \code{chol}).
#'
#' There are three further parameters you can pass to \code{write_morpho} through
#' \code{ ... }: \code{R.sh}, a shrunk matrix, \code{scaled}, a logical value
#' indicating if the object \code{proc} has been already scaled, and a numeric value
#' for the maximum age, which is always 1.
#' These arguments are used within the function \code{\link{simulate_morpho}}
#' and are not thought to be called by the user.
#'
#' @seealso
#' \code{\link{matrix2array}}, \code{\link{array2matrix}}, \code{\link{simulate_morpho}}
#'
#' @author Sandra Alvarez-Carretero
#'
#' @examples
#' # A.1) Providing only the morphological alignment (proc) after the
#' #      Procrustes analysis (PA) in an object of class \"array\".
#'
#'        write_morpho( filename = "seqfile.aln", proc = canids19x29.array,
#'                      coords = 3 )
#'
#' # A.2) Providing only the morphological alignment (proc) after the
#' #      Procrustes analysis (PA) in an object of class \"matrix\".
#'
#'        write_morpho( filename = "seqfile.aln", proc = canids19x29.matrix )
#'
#' # B) Providing the morphological alignment (proc) after the
#' #    PA in an object of class \"array\" and a list with the
#' #    names of the species
#'
#'      names <- list( sp1  = "Ael_sp.", sp2  = "Can_dir", sp3  = "Epi_hay", sp4  = "Hes_sp.",
#'                     sp5  = "Mes_cor", sp6  = "Tom_sp.", sp7  = "Enh_pah", sp8  = "Cuo_alp",
#'                     sp9  = "Spe_ven", sp10 = "Can_lup", sp11 = "Cer_tho", sp12 = "Oto_meg",
#'                     sp13 = "Vul_vul", sp14 = "Urs_ame", sp15 = "Ail_ful", sp16 = "Nan_bio",
#'                     sp17 = "Par_her", sp18 = "Tha_won", sp19 = "Smi_fat"
#'                   )
#'
#'      write_morpho( filename = "seqfile.aln", proc = canids19x29.array,
#'                    coords = 3, names = names )
#'
#' # C) Providing the morphological alignment (proc) after the
#' #    PA in an object of class \"array\", a list with the names
#' #    of the species, and a list with the ages of the species
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
#'      write_morpho( filename = "seqfile.aln", proc = canids19x29.array,
#'                    coords = 3, names = names, ages = ages )
#'
#' # D) Providing an object of class \"array\" after having
#' #    carried out a PA. As an example, we use the function
#' #    geomorph::gpagen, but you can use your preferred
#' #    function meanwhile the format of the \"proc\" object
#' #    (class array) for write_morpho is
#' #    "p landmarks x k coordinates x n species"
#'
#'      df <- canids19x29.raw
#'
#'      # Process raw data
#'
#'      mm <- df[,2:dim( df )[2]] # 19 sp x 87 coords (87/3=29 lmks)
#'      rownames( mm ) <- df[,1]
#'      mm <- as.matrix( mm )
#'
#'      # Convert matrix into 3D array
#'
#'      ma <- matrix2array( proc = mm, coords = 3 )
#'
#'      # Get procrustes analysis done
#'
#'      ma.paln <- geomorph::gpagen( ma )
#'
#'      # Check object is class \"array\"
#'
#'      class( ma.paln$coords )
#'
#'      # Run write_morpho
#'
#'      write_morpho( filename = "seqfile.aln", proc = ma.paln$coords, coords = 3 )
#'
#' # E) Providing the morphological alignment (proc) after the
#' #    PA in an object of class \"matrix\" and population variance
#'
#'      write_morpho( filename = "seqfile.aln", proc = canids19x29.array,
#'                    coords = 3, popvar = var.fx )
#'
#' # F) Providing the morphological alignment (proc) after the
#' #    PA in an object of class \"matrix\", and the population variance
#' #    and the correlation matrix (choosing the Cholesky
#' #    decomposition)
#'
#'      write_morpho( filename = "seqfile.aln", proc = canids19x29.array,
#'                    coords = 3, popvar = var.fx, R = R, method = "chol" )
#'
#'
#'
#'
#' @export

write_morpho <- function( filename, proc, coords = c( 2, 3 ), names = NULL, ages = NULL,
                          popvar = NULL, R = NULL, method = c( "eigen", "chol" ), ... ) {

  #\\ Check a name for the output file has been given

  if ( missing( filename ) ){
    stop( "Please add a name for the output file" )
  }

  #\\ Check object with lmks is provided

  if ( missing( proc ) ){
         stop( "Please use an object of class \"matrix\" or \"array\" with the result of a Procrustes analysis" )
  }

  if ( class( proc ) == "array" | class( proc ) == "matrix" ){

    #\\ Check class of object with lmks and
    #\\ get variables

    if( class( proc ) == "array" ){

      # The coords info is needed

      if ( is.null( coords ) | length( coords ) != 1 | class( coords ) != "numeric" ){
        stop( "Please enter coords = 2 or coords = 3 depending on using 2D or 3D landmarks, respectively" )
      }

      chars  <- dim( proc )[1] * dim( proc )[2]
      num.sp <- dim( proc )[3]

    }

    else if( class( proc ) == "matrix" ){

      num.sp <- dim( proc )[1]
      chars  <- dim( proc )[2]

    }

    #\\ Get passed arguments through ...

    simul.pars <- list( ... )

    #\\ If names are provided...

    if ( !is.null( names ) ){

      # Check class for names

      if ( class( names ) != "list" ){

        stop( "You need to provide an object of class list with the species included in the alignment file.
               E.g. species <- list( sp1 = \"sp1\", sp2 = \"sp2\"" )
      }

      # Get species names

      names   <- unlist( names )

      #-- If ages are not provided ...

      if ( is.null( ages ) ){

        # Get spaces to justify text

        num.spaces <- max( nchar( names ) ) - nchar( names )
        x          <- sprintf( paste( "% ", num.spaces, "s", sep = ""), c( "" ) )

      }

      #-- If ages are provided ...

      else{

        # Check class for ages

        if ( class( ages ) != "list" ){

          stop( "You need to provide an object of class list with the ages of the species included in the alignment file.
                 E.g. ages <- list( ag1=30, ag2=15" )
        }

        # Put together names and transformed ages and get spaces

        ages       <- unlist( ages )
        if ( "max.age" %in% names( simul.pars ) ){
          names      <- paste( names, - ages + simul.pars$max.age, sep = "^" )
        }
        else{
          names      <- paste( names, - ages + max( ages ) + 0.01, sep = "^" ) # 0.01 = ct for MCMCTree
        }
        num.spaces <- max( nchar( names ) ) - nchar( names )
        x          <- sprintf( paste( "% ", num.spaces, "s", sep = "" ), c( "" ) )

      }

    }


    #\\ If names are not provided...

    else{

      # Name species as "Species_1", "Species_2", etc.

      names <- paste( "Species_", seq( 1:num.sp ), sep = "" )

      #-- If ages are not provided ...

      if ( is.null(ages) ){

        # Get spaces to justify text

        num.spaces <- max( nchar( names ) ) - nchar( names )
        x          <- sprintf( paste( "% ", num.spaces, "s", sep = "" ), c( "" ) )

      }

      #-- If ages are provided ...

      else{

        # Check class for ages

        if ( class(ages) != "list" ){

          stop( "You need to provide an object of class list with the ages of the species included in the alignment file.
              E.g. ages <- list(ag1=30, ag2=15" )
        }

        # Put together names and transformed ages and get spaces

        ages       <- unlist( ages )
        if ( "max.age" %in% names( simul.pars ) ){
          names      <- paste( names, - ages + simul.pars$max.age, sep = "^" )
        }
        else{
          names      <- paste( names, - ages + max( ages ) + 0.01, sep = "^" ) # 0.01 = ct for MCMCTree
        }
        num.spaces <- max( nchar(names) ) - nchar( names )
        x          <- sprintf( paste( "% ", num.spaces, "s", sep = "" ), c( "" ) )

      }

    }

    #\\ Get names with data

    names <- paste( names, x, sep = "    " )

    if ( class( proc ) != "array" ){

      rownames( proc ) <- names

    }

    else{

      faln             <- array2matrix( proc = proc, coords = coords )
      rownames( faln ) <- names
      proc             <- faln

    }

    #\\ Get popvar

    if ( !is.null( popvar ) ){

      # If the variances are different...
      if ( length( unique( popvar ) ) != 1 ){

        proc <- proc %*% diag( 1 / sqrt( popvar ) )
        pvar <- 1 # The matrix has been scaled, so now

      }

      # If the variances are all equal...
      if ( length( unique( popvar ) ) == 1 ){

        pvar <- 1 # No need for scale

      }

    }

    else{

      pvar <- 0 # No variance at all

    }

    #\\ Decompose R and get lnd

    if ( !is.null( R ) ){

      # Check method is not missing

      if ( missing(method) ){
        stop( "Please select a method to decompose the correlation matrix")
      }

      # If R is not class "matrix", stop

      if ( class( R ) != "matrix" ){

        stop( "Please provide an object of class \"matrix\" for the correlation matrix" )
      }

      # Match argument

      method <- match.arg(method)

      if ( method == "chol" ){

        # chol returns upper triangular matrix

        U <- chol( R )
        L <- t( U )

        # R = L %*% U = L %*% t( L )
        # solve( R ) = t( solve(L) ) %*% solve( L )
        # solve( R ) = t( A )        %*% A
        #
        # all.equal( t( solve( L ) ) %*% solve( L ), solve( R ) )

        Linv <- solve( L )

        # Transform data "M", object proc, so
        # Z = M %*% t( A ) = proc %*% t( Linv )

        Z <- proc %*% t( Linv )

      }

      else if ( method == "eigen" ){

        # eigen returns a list with the eigenvectors
        # and the eigenvalues
        #
        # solve( R ) = Rinv = t( A ) %*% A
        #
        # t( A ) = V %*% D
        # t( A ) = eigen( Rinv )$vectors %*% diag( sqrt( eigen( Rinv )$values ) )
        #
        # all.equal( tA %*% t( tA ), solve( R ) )

        tA <- eigen( solve(R) )$vectors %*% diag( sqrt( eigen( solve(R) )$values ) )

        # Transform data "M", object proc, so
        # Z = M %*% t( A ) = proc %*% t( A )

        Z <- proc %*% tA

      }

      # Get logarithm of the determinant of the correlation matrix
      # This variable is used to correct the likelihood value when the data have
      # been transformed to account for the correlation.
      # This is done to allow the correct calculation of bayes factors.

      lnd <- determinant( R )$modulus # lnd = lod(det(R))

    }

    else{

      # Check now if further arguments have been passed
      # just in case the matrix has simulated data

      if ( length( simul.pars ) != 0 ){

        if ( "R.sh" %in% names( simul.pars ) ){
          # cat( "You provided the shrunk matrix. The logarithm
          #        is being calcualted", "\n")
          R.sh <- simul.pars$R.sh
          lnd <- determinant( R.sh )$modulus

        }
        else{
          lnd <- 0
        }

        if ( "scaled" %in% names( simul.pars ) ){

          scaled.log <- simul.pars$scaled

          if ( scaled.log == TRUE ){
            # cat( "The matrix you provided has been previously
            #       scaled.", "\n")
            pvar <- 1
          }

          else{
            pvar <- 0
          }

        }
        else{
          pvar <- 0
        }

        Z <- proc

      }

      # There is no correlation at all and we are not
      # dealing with simulated data

      else{

        lnd <- 0
        Z   <- proc

      }

    }

    #\\ Generate file

    string <- paste( length( names ), chars, "M", pvar, lnd, sep = "  " )
    cat( "\n", file = filename )
    write( paste( "  ", string, sep = "" ), file = filename, append = T )
    #cat(paste("  ", string, sep=""))
    cat( "\n", file = filename, append = T )
    #cat("\n")
    write.table( Z, file = filename, append = T,
                 sep = " ", row.names = T, col.names = F, quote = FALSE )
    cat( "\n", file = filename, append = T )
    #Z

  }

  #\\ If object proc is not a matrix nor an array...

  else{

    stop( "You need to use an object of class \"matrix\" or \"array\" with the result of a Procrustes analysis" )

  }

}






#' Convert a matrix into a 3D array
#'
#' @description
#' Convert a matrix with landmark points into a 3D array object.
#'
#' @param proc Matrix (n x p), 'p' landmark points for 'n' specimens
#' (see details).
#'
#' @param coords Integer, 2 or 3 for 2D or 3D landmarks, respectively.
#'
#' @details
#'
#' The matrix has format n x p, with 'n' rows, one for specimen, and 'p' columns.
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
#' See \code{canids19x29.matrix} in the
#' \code{data} directory for an example.
#'
#' @return
#'
#' An object of class array with format p x k x n, where 'p' is the number of
#' landmarks, 'k' the number of coordinates, and 'n' the number of specimens.
#'
#' Note that if the matrix provided does not have rownames, the specimens in the
#' returned array (dimension 'n') will be labelled as '1', '2', and so on.
#' See \code{canids19x29.array} in the \code{data} directory for an example of
#' the format of the object that is returned.
#'
#' @seealso
#' \code{\link{array2matrix}}, \code{\link{write_morpho}}
#'
#' @author Sandra Alvarez-Carretero
#'
#' @export

matrix2array <- function( proc, coords = c( 2, 3 ) ){

  #\\ Check object with lmks is provided

  if ( missing( proc ) ){
    stop( "Please provide an object of class \"matrix\"" )
  }

  #\\ Only proceed it the object is a matrix

  if ( class( proc ) == "matrix" ){

    if ( length( dim( proc ) ) != 2 ){
      stop( "Please provide a matrix \"n x p\", with 'n' specimens and 'p' characters" )
    }

    if ( is.null( coords ) | length( coords ) != 1 | class( coords ) != "numeric" ){
      stop( "Please enter coords = 2 or coords = 3 depending on using 2D or 3D landmarks, respectively" )
    }

    chars       <- dim( proc )[2]
    num.sp      <- dim( proc )[1]
    list.coords <- c("x", "y", "z")

    # Create an empty array to store the coordinates in the format
    # p x k x n (chars x coord x num.sp)

    ma             <- array( dim = c( chars / coords, coords, num.sp ) )
    dimnames( ma ) <- list( paste( "lmk", seq( 1:( chars/coords ) ), sep="" ),
                            list.coords[1:coords],
                            rownames( proc ) )

    # Select x, y, z positions

    if ( coords == 3 ){

      xi <- seq( from = 1, to = chars, by = 3 )
      yi <- seq( from = 2, to = chars, by = 3 )
      zi <- seq( from = 3, to = chars, by = 3 )

      # Fill in array p x k x n

      for (i in 1:num.sp) {
        ma[,1,i] <- unlist( proc[i,xi] )
        ma[,2,i] <- unlist( proc[i,yi] )
        ma[,3,i] <- unlist( proc[i,zi] )
      }

    }

    else if ( coords == 2 ){

      xi <- seq( from = 1, to = chars, by = 2 )
      yi <- seq( from = 2, to = chars, by = 2 )

      # Fill in array p x k x n

      for (i in 1:num.sp) {
        ma[,1,i] <- unlist( proc[i,xi] )
        ma[,2,i] <- unlist( proc[i,yi] )
      }

    }


  #\\ Return 3D array

    return(ma)

  }

  else{
    stop( "Please use an object of class \"matrix\"" )
  }

}







#' Convert a 3D array into a matrix
#'
#' @description
#' Convert a 3D array with landmark points into an object of class matrix.
#'
#' @param proc 3D array, 'p' landmark points, 'k' coordinates, and
#' 'n' specimens.
#'
#' @param coords Integer, 2 or 3, for 2D or 3D landmarks, respectively.
#'
#' @details
#'
#' The object \code{proc}, class array, has format p x k x n, where 'p' is
#' the number of landmarks, 'k' the number of coordinates, and 'n' the number
#' of specimens. See \code{canids19x29.array} in the \code{data} directory for an
#' example of the format of a 3D array.
#'
#' @return
#'
#' An object of class matrix, with 'n' rows, one for specimen, and 'p' columns.
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
#' more detailed example of the format of the object that is returned.
#'
#' Note that if the 'n' dimension (specimens) of the array provided does
#' not have names, the specimens in the returned matrix will be
#' labelled as '1', '2', and so on.
#'
#' @seealso
#' \code{\link{matrix2array}}, \code{\link{write_morpho}}
#'
#' @author Sandra Alvarez-Carretero
#'
#' @export

array2matrix <- function( proc, coords = c( 2, 3 ) ){

  #\\ Check object with lmks is provided

  if ( missing( proc ) ){
    stop( "Please provide an object of class \"array\"" )
  }

  #\\ Only proceed if the object is an array

  if ( class( proc ) == "array" ){

    if ( length( dim( proc ) ) != 3 ){
      stop( "Please provide a 3D array \"p x k x n\", with 'p' landmarks, 'k' coordinates, and 'n' specimens" )
    }

    if ( is.null( coords ) | length( coords ) != 1 | class( coords ) != "numeric" ){
      stop( "Please enter coords = 2 or coords = 3 depending on using 2D or 3D landmarks, respectively" )
    }

    chars  <- dim( proc )[1]
    coords <- coords
    num.sp <- dim( proc )[3]

    faln <- matrix( 0, nrow = num.sp, ncol = chars * coords )

    if ( coords == 2 ){

      xi <- seq( from = 1, to = chars, by = 2 ); yi <- seq( from = 2, to = chars, by = 2 )

      for ( i in 1:num.sp ) {
        faln[i,xi] <- proc[,1,i]
        faln[i,yi] <- proc[,2,i]
      }

      rownames(faln) <- dimnames( proc )[[ 3 ]]

    }

    else if ( coords == 3 ){
      xi <- seq( from = 1, to = chars * coords, by = 3 )
      yi <- seq( from = 2, to = chars * coords, by = 3 )
      zi <- seq( from = 3, to = chars * coords, by = 3 )

      for ( i in 1:num.sp ) {
        faln[i,xi] <- proc[,1,i]
        faln[i,yi] <- proc[,2,i]
        faln[i,zi] <- proc[,3,i]

      }

      rownames(faln) <- dimnames( proc )[[ 3 ]]

    }

  #\\ Return matrix

    return(faln)

  }

  else{
    stop( "Please use an object of class \"array\"" )
  }

}


#' Simulate a morphological alignment
#'
#' @description
#' Simulate a continuous morphological alignment using \code{\link[ape]{rTraitCont}}
#' and later allowing to account for population variance and trait correlation.
#'
#' @param tree Phylo, object with a phylogenetic tree
#' (see \code{\link[ape]{rTraitCont}}).
#'
#' @param mtraits Numeric, number of morphological traits to be simulated.
#'
#' @param reps Numeric, number of replicates to be simulated.
#'
#' @param c (Optional) numeric, population variance to add to the simulated
#' morphological traits. Requires \code{psample} (see details).
#'
#' @param rho (Optional) numeric, correlation value used to generate a
#' correlation matrix. Requires \code{c}, \code{psample}, and \code{method}
#' (see details).
#'
#' @param R (Optional) matrix, correlation matrix. Requires \code{c},
#' \code{psample}, and \code{method} (see details).
#'
#' @param psample (Optional) numeric, number of individuals to sample
#' from the simulated population.
#'
#' @param method (Optional) character, either \code{"eigen"} or
#' \code{"col"}, method used to decompose the inverse of the shrunk
#' correlation matrix. Requires \code{c}, \code{psample}, and either \code{R}
#' or \code{rho} (see details).
#'
#' @param out (Optional) character, name for the output file with the
#' simulated data in phylip format (see details and \code{\link{write_morpho}}).
#'
#' @param names (Optional) list, species name included in the
#' morphological alignment (see examples in \code{\link{write_morpho}}).
#'
#' @param ages (Optional) list, ages of the species included in
#' the morpholical alignment (see examples in \code{\link{write_morpho}}).
#'
#' @param ... Further arguments passed to \code{\link[ape]{rTraitCont}}.
#'
#' @details
#'
#' The function \code{\link[ape]{rTraitCont}} simulates continuous traits and
#' can take different parameters to adjust the simulation
#' (e.g. the model, the rate drift, etc.).
#' These parameters are the ones the user can pass through
#' \code{simulate_morpho()}.
#' The default values that \code{simulate_morpho()} uses are
#' \code{model = "BM"}, \code{sigma = 1}, \code{ancestor = F},
#' and \code{root.value = 0}. Currently, \code{simulate_morpho}
#' supports only \code{ancestor = F}, so do not change this logical
#' value or the function will not work. See \code{\link[ape]{rTraitCont}}
#' for more details on this function.
#'
#' The parameter \code{c} is the population variance. If a correlation
#' matrix, \code{R}, or the parameter \code{rho} needed to generate
#' a correlation matrix following the constant correlation model are
#' not provided, a total of \code{s x p} samples per replicate are
#' generated following a normal distribution with mean 0
#' and variance \code{c}, where \code{s} is the number of specimens and
#' \code{p} the number of morphological traits, \code{mtraits}.
#' Otherwise, if either \code{R} or \code{rho} are provided, then the
#' \code{s x p} samples follow a multivariate
#' normal distribution with mean 0 and variance \code{c R}.
#' The resulting randomly sampled variables are used to obtain one noise
#' matrix per replicate \eqn{i}, \code{N_i}, which is added to the
#' corresponding simulated matrix, \code{M_i}, such as
#' \eqn{\mathrm{M.n_{i}}=\mathrm{M_{i}}+\mathrm{N_{i}}}{M.n_i = M_i + N_i}.
#' For each replicate, this results into the noisy matrix \eqn{i},
#' \code{M.n_i}, which accounts for population noise.
#'
#' The parameter \code{rho} is used to generate one correlation matrix per
#' replicate, \code{R_i}, with dimensions \code{p x p}, where \code{p} is
#' the number of continuous traits, \code{mtraits}. All elements in each
#' \code{R_i} generated have value \code{rho} (constant correlation model).
#' If you want to use another correlation matrix, please do not provide
#' any numeric value to the parameter \code{rho}. Just provide your
#' preferred correlation matrix as the argument of the parameter \code{R},
#' class "matrix".
#'
#' In order to account for population variance, the parameter \code{psample}
#' is required. This parameter indicates the number of individuals \code{n}
#' needed to generated one population matrix per replicate, \code{P_i},
#' with dimensions \code{n x p}, where \code{p} is the number of
#' continuous traits to be simulated, \code{mtraits}. The simulated
#' continuous traits for each sampled population per replicate follow a
#' normal distribution with mean 0 and variance \code{c} (required parameter).
#' The variance of each \code{P_i} is then calculated and used to scale
#' the corresponding noisy matrix, \code{M.n_i}, such as
#' \eqn{\mathrm{M.s_{i}}=\mathrm{M.n_{i}}\times diag\left(\frac{1}{
#' \sqrt{diag(cov(\mathrm{P_{i}}))}}\right)}{
#' M.s_i = M.n_i diag( 1 / \sqrt( diag(cov(P_i)) ) )}.
#'
#' If either the parameter \code{rho} or \code{R} are provided, trait
#' correlation is considered. Therefore, each scaled simulated matrix with
#' continuous traits, \code{M.s_i}, needs to be transformed in order
#' to account for this correlation.
#' Specifically, \code{simulate_morpho} estimates the shrunk correlation
#' matrix for each replciate, \code{R.sh_i}, with the function
#' \code{\link[corpcor]{cor.shrink}}. Later, it decomposes the inverse
#' of \code{R.sh_i} either using the Cholesky decomposition
#' (if \code{method = "chol"}, see default usage at \code{\link[base]{chol}})
#' or the eigendecomposition (if \code{method = "eigen"}, see default
#' usage at \code{\link[base]{eigen}}), such as
#' \eqn{\mathrm{R.sh_{i}}=\mathrm{A_{i}^{T}}\mathrm{A_{i}}}{R.sh_i = t(A_i) A_i}.
#' Each matrix \code{t(A_i)} is used to transform the corresponding \code{M.s_i}
#' such as \eqn{\mathrm{Z_{i}}=\mathrm{M.s_{i}}\times\mathrm{A{i}^{T}}}{Z_i =
#' M.s_i t(A_i)},
#' where each matrix \code{Z_i} is the transformed data set, a matrix that
#' has been scaled and accounts for trait correlation.
#'
#' If the user wants to output the resulting simulated alignment in
#' phylip format readable by MCMCTree, then a name for the output file
#' should be provided as the argument of \code{out}.
#'
#' @return
#'
#' A) If \code{c} nor either \code{rho} or \code{R} are provided, then
#' a list with a list of each simulated matrix with continuous traits
#' per replicate is returned:
#'   \item{M}{  List with \eqn{i} matrices \code{s x p}, with \code{p} simulated continuous
#' traits for \code{s} specimens}
#'
#' B) If \code{c} is provided, then a list with a list of the following
#' matrices, one per replicate, is returned: \eqn{i} simulated matrices
#' with continuous traits, \eqn{i} noise matrices, \eqn{i}
#' population matrices, \eqn{i} noisy matrices, and \eqn{i} scaled matrices.
#'   \item{M}{  List with \eqn{i} matrices \code{s x p}, with \code{p}
#'   simulated continuous traits for \code{s} specimens}
#'   \item{N}{  List with \eqn{i} matrices \code{s x p}, with \code{s x p}
#'   randomly sampled values from a normal distribution with mean 0
#'   and variance \code{c}}
#'   \item{P}{  List with \eqn{i} matrices \code{n x p}, where \code{n}
#'   is the amount of individuals sampled from a population with \code{p}
#'   continuous traits simulated under a normal distribution
#'   with mean 0 and variance \code{c}}
#'   \item{M.n}{  List with \eqn{i} matrices \code{s x p}, result of
#'   \eqn{\mathrm{M.n_{i}}=\mathrm{M_{i}}+\mathrm{N_{i}}}{M.n_i = M_i + N_i}}
#'   \item{M.s}{  List with \eqn{i} matrices \code{s x p}, result of
#'   \eqn{\mathrm{M.s_{i}}=\mathrm{M.n_{i}}\times diag\left(\frac{1}{\sqrt{
#'   diag(cov(\mathrm{P_{i}}))}}\right)}{
#'   M.s_i = M.n_i diag( 1 / \sqrt( diag(cov(P_i)) ) )}}
#'
#' C) If \code{c} and either \code{rho} or \code{R} are provided,
#' then a list with a list of the following matrices, one per replicate,
#' is returned: \eqn{i} simulated matrices with continuous traits,
#' \eqn{i} noise matrices, \eqn{i} population matrices, \eqn{i} noisy
#' matrices, \eqn{i} scaled matrices, \eqn{i} shrunk correlation matrices,
#' and \eqn{i} transformed matrices.
#'   \item{M}{  List with \eqn{i} matrices \code{s x p}, with \code{p}
#'   simulated continuous traits for \code{s} specimens}
#'   \item{N}{  List with \eqn{i} matrices \code{s x p}, with \code{s x p}
#'   randomly sampled values from a normal distribution with mean 0
#'   and variance \code{c}}
#'   \item{P}{  List with \eqn{i} matrices \code{n x p}, where \code{n}
#'   is the amount of individuals sampled from a population with \code{p}
#'   continuous traits simulated under a normal distribution
#'   with mean 0 and variance \code{c}}
#'   \item{M.n}{  List with \eqn{i} matrices \code{s x p}, result of
#'   \eqn{\mathrm{M.n_{i}}=\mathrm{M_{i}}+\mathrm{N_{i}}}{M.n_i = M_i + N_i}}
#'   \item{M.s}{  List with \eqn{i} matrices \code{s x p}, result of
#'   \eqn{\mathrm{M.s_{i}}=\mathrm{M.n_{i}}\times diag\left(\frac{1}{\sqrt{
#'   diag(cov(\mathrm{P_{i}}))}}\right)}{
#'   M.s_i = M.n_i diag( 1 / \sqrt( diag(cov(P_i)) ) )}}
#'   \item{R.sh}{  List with \eqn{i} matrices \code{p x p}, where \code{p}
#'   is the amount of simulated continuous traits. These matrices are the
#'   estimated shrunk correlation matrices computed for each replicate
#'   with the function \code{\link[corpcor]{cor.shrink}}}
#'   \item{Z}{  List of \eqn{i} matrices \code{s x p}, matrices with the
#'   transformed data calculated as \eqn{\mathrm{Z_{i}}=\mathrm{M.s_{i}}\times
#'   \mathrm{A_{i}^{T}}}{Z_i = M.s_i t(A_i)}}
#'
#' @seealso
#' \code{\link{write_morpho}}
#'
#' @author Sandra Alvarez-Carretero
#'
#' @examples
#' # A) Simulation setup: Simulate a morphological alignment
#' #    with 'mtraits' = 87 continuous characters that follows a
#' #    fixed tree, object 'tree', with the default parameters in
#' #    'simulate_morpho' to run 'rTraitCont'.
#' #    Population variance and correlation are not considered.
#' #
#' #    Number of replicates: 2.
#'
#'      morpho::simulate_morpho( tree = tree, mtraits = 87, reps = 2 )
#'
#'
#' # B) Simulation setup: Simulate a morphological alignment
#' #    with 'mtraits' = 87 continuous characters that follows a
#' #    fixed tree, object 'tree', but with different parameters
#' #    than the default ones in 'simulate_morpho' to run
#' #    'rTraitCont'. Population variance and correlation are not
#' #    considered.
#' #
#' #    Number of replicates: 2.
#'
#'      morpho::simulate_morpho( tree  = tree, mtraits = 87, reps  = 2,
#'                               model = "OU", root    = 1,  sigma = 0.2,
#'                               alpha = 2
#'                              )
#'
#'
#' # C) Simulation setup: Simulate a morphological alignment
#' #    with 'mtraits' = 87 continuous characters that follows a
#' #    fixed tree, object 'tree', with the default parameters in
#' #    'simulate_morpho' to run 'rTraitCont'.
#' #    Population variance is c = 0.2 and the within population from
#' #    where to sample has 'psample' = 5 individuals.
#' #    Correlation is not considered.
#' #
#' #    Number of replicates: 2.
#'
#'      morpho::simulate_morpho( tree = tree, mtraits = 87, reps = 2,
#'                               c    = 0.2,  psample = 5
#'                              )
#'
#'
#' # D) Simulation setup: Simulate a morphological alignment
#' #    with 'mtraits' = 87 continuous characters that follows a
#' #    fixed tree, object 'tree', with the default parameters in
#' #    'simulate_morpho' to run 'rTraitCont'.
#' #    Population variance is c = 0.2 and the within population from
#' #    where to sample has 'psample' = 5 individuals.
#' #    A correlation matrix is provided, so once it is shrunk, it can
#' #    be decomposed using the Cholesky decomposition method.
#' #
#' #    Number of replicates: 2.
#'
#'
#'      morpho::simulate_morpho( tree   = tree, mtraits = 87, reps = 2,
#'                               c      = 0.2,  psample = 5,  R    = R,
#'                               method = "chol"
#'                              )
#'
#'
#' # E) Simulation setup: Simulate a morphological alignment
#' #    with 'mtraits' = 87 continuous characters that follows a
#' #    fixed tree, object 'tree', with the default parameters in
#' #    'simulate_morpho' to run 'rTraitCont'.
#' #    Population variance is c = 0.2 and the within population from
#' #    where to sample has 'psample' = 5 individuals.
#' #    A correlation matrix is provided, so once it is shrunk, it can
#' #    be decomposed using the eigen decomposition method.
#' #
#' #    Number of replicates: 2.
#'
#'      morpho::simulate_morpho( tree   = tree, mtraits = 87, reps = 2,
#'                               c      = 0.2,  psample = 5,  R    = R,
#'                               method = "eigen"
#'                              )
#'
#'
#' # F) Simulation setup: Simulate a morphological alignment
#' #    with 'mtraits' = 87 continuous characters that follows a
#' #    fixed tree, object 'tree', with the default parameters in
#' #    'simulate_morpho' to run 'rTraitCont'.
#' #    Population variance is c = 0.2 and the within population from
#' #    where to sample has 'psample' = 5 individuals.
#' #    The rho parameter, 'rho' = 0.5, is provided so a correlation
#' #    matrix following the constant correlation model is generated, which
#' #    later will be shrunk. The shrunk matrix will be decomposed
#' #    using the shrunk correlation matrix.
#' #
#' #    Number of replicates: 2.
#'
#'      morpho::simulate_morpho( tree   = tree, mtraits = 87, reps = 2,
#'                               c      = 0.2,  psample = 5,  rho  = 0.5,
#'                               method = "chol"
#'                              )
#'
#'
#' # G) Simulation setup: Simulate a morphological alignment
#' #    with 'mtraits' = 87 continuous characters that follows a
#' #    fixed tree, object 'tree', with the default parameters in
#' #    'simulate_morpho' to run 'rTraitCont'.
#' #    Population variance is c = 0.2 and the within population from
#' #    where to sample has 'psample' = 5 individuals.
#' #    The rho parameter, 'rho' = 0.5, is provided so a correlation
#' #    matrix following the constant correlation model is generated, which
#' #    later will be shrunk. The shrunk matrix will be decomposed
#' #    using the shrunk correlation matrix.
#' #    Set specific names and ages for the simulated species. The order
#' #    follows the order in object 'tree'.
#' #
#' #    Number of replicates: 2.
#'
#'
#'      names <- list( sp1  = "A", sp2  = "B", sp3  = "F", sp4  = "C",
#'                     sp5  = "H", sp6  = "D", sp7  = "G", sp8  = "E"
#'                    )
#'
#'      ages <- list( sp1 = 0,   sp2 = 0, sp3 = 0.1, sp4 = 0,
#'                    sp5 = 0.7, sp6 = 0, sp7 = 0.3, sp8 = 0
#'                   )
#'
#'      morpho::simulate_morpho( tree   = tree,   mtraits = 87,   reps = 2,
#'                               c      = 0.2,    psample = 5,    rho  = 0.5,
#'                               method = "chol", names   = names,
#'                               ages   = ages
#'                              )
#'
#' @export

simulate_morpho <- function( tree, mtraits, reps, c = NULL,
                             rho = NULL, R = NULL, psample = NULL,
                             method = c( "eigen", "chol" ),
                             out = NULL, names = NULL, ages = NULL,
                             ... ){

  #\\ Check objects tree, mtraits, and reps are provided

  if ( missing( tree ) ){
    stop( "Please provide an object of class \"phylo\" with a phylogenetic tree" )
  }

  if ( class(tree) != "phylo" ){
    stop( "Please use an object of class \"phylo\" with a phylogenetic tree" )
  }


  if ( missing( mtraits ) ){
    stop( "Please provide a numeric value with the amount of morphological
          continuous traits to be simulated" )
  }

  if ( missing( reps ) ){
    stop( "Please provide a numeric value with the amount of replicates
          you want to be simulated" )
  }

  #\\ Get amount of specimens

  s <- length(tree$tip.label)

  #\\ Get ... parameters in a list in case the user has provided
  #\\ them. Else, assign default values

  rtraitcont.pars <- list( ... )

  if ( "model" %in% names( rtraitcont.pars ) ){
    model <- rtraitcont.pars$model
  }
  else{
    model <- "BM"
  }

  if ( "sigma" %in% names( rtraitcont.pars ) ){
    sigma <- rtraitcont.pars$sigma
  }
  else{
    sigma <- 1
  }

  if ( "alpha" %in% names( rtraitcont.pars ) ){
    alpha <- rtraitcont.pars$alpha
  }
  else{
    alpha <- 1
  }

  if ( "theta" %in% names( rtraitcont.pars ) ){
    theta <- rtraitcont.pars$theta
  }
  else{
    theta <- 0
  }

  if ( "ancestor" %in% names( rtraitcont.pars ) ){
    ancestor <- rtraitcont.pars$ancestor
  }
  else{
    ancestor <- FALSE
  }

  if ( "root.value" %in% names( rtraitcont.pars ) ){
    root <- rtraitcont.pars$root.value
  }
  else{
    root <- 0
  }

  cat( "\n", "Simulating the data using the following parameters: ", "\n\n",
       "\tModel:    ", model,"\n",
       "\tSigma:    ", sigma, "\n",
       "\tAlpha:    ", alpha, "\n",
       "\tTheta:    ", theta, "\n",
       "\tRoot:     ", root, "\n",
       "\tAncestor: ", ancestor, "\n\n"
        )

  #\\ Simulate continous data with ape::rTraitCont
  #\\ Save everything in a 3D array with dimensions
  #\\ p (reps) x k (mtraits) x n (num. specimens)

  M             <- array( dim = c( s, mtraits, reps ) )
  dimnames( M ) <- list( tree$tip.label,
                         paste( "trait_", seq( 1:( mtraits ) ), sep="" ),
                         paste( "rep_", seq( 1:( reps ) ), sep="" )
                        )

  for ( i in seq( 1:reps ) ){

    # Add replicate i into array M

    M[ , , i] <- replicate( mtraits,
                            ape::rTraitCont( phy      = tree,     model      = model,
                                             sigma    = sigma,    root.value = root,
                                             ancestor = ancestor, alpha      = alpha,
                                             theta    = theta
                                            )
                           )
  }

  #\\ If population variance is provided ...

  if ( !is.null( c ) ){

    cat( "\n", "You have provided the population variance, c =", c, "\n" )

    # A) If the correlation matrix is provided ...

    if ( !is.null( R ) ){

      cat( "\n", "You have provided the correlation matrix", "\n" )

      if ( class( R ) != "matrix" ){
        stop( "Please use a correlation matrix of class \"matrix\"" )
      }

      # Sample num.species*traits samples from a multiv.norm.
      # distrib with mu = 0 and sd = c*R
      # ( Array with i noise matrices with s x traits dimensions )

      N             <- array( dim = c( s, mtraits, reps ) )
      dimnames( N ) <- list( paste( "Sp_", seq( 1:( s ) ), sep="" ),
                             paste( "trait_", seq( 1:( mtraits ) ), sep="" ),
                             paste( "rep_", seq( 1:( reps ) ), sep="" )
                            )

      for ( i in seq( 1:reps ) ){
        N[ , , i] <- t( mvtnorm::rmvnorm( s, mean  = rep( 0, nrow( R ) ),
                                          sigma = c*R )
                       )
      }

      # Generate population sample
      # ( 3D array with i pop matrices with psample x traits dimensions )

      P             <- array( dim = c( psample, mtraits, reps ) )
      dimnames( P ) <- list( paste( "Sp_", seq( 1:( psample ) ), sep="" ),
                                  paste( "trait_", seq( 1:( mtraits ) ), sep="" ),
                                  paste( "rep_", seq( 1:( reps ) ), sep="" )
                                 )

      for ( i in seq( 1:reps ) ){
        P[ , , i] <- t( mvtnorm::rmvnorm( psample, mean = rep( 0, nrow( R ) ),
                                          sigma = c*R )
                       )

      }

    }

    # B) If the correlation matrix is not provided and
    #    a matrix with the same rho parameter for all correlations
    #    is wanted...

    else if ( !is.null( rho ) ){

      cat( "\n", "You have provided the rho parameter", rho, "to generate the correlation matrix following the constant correlation model", "\n")

      R       <- matrix( rep( rho ), ncol = mtraits, nrow = mtraits )
      diag(R) <- 1

      # Sample num.species*traits samples from a multiv.norm.
      # distrib with mu = 0 and sd = c*R
      # ( Array with i noise matrices with s x traits dimensions )

      N             <- array( dim = c( s, mtraits, reps ) )
      dimnames( N ) <- list( paste( "Sp_", seq( 1:( s ) ), sep="" ),
                             paste( "trait_", seq( 1:( mtraits ) ), sep="" ),
                             paste( "rep_", seq( 1:( reps ) ), sep="" )
                            )

      for ( i in seq( 1:reps ) ){
        N[ , , i] <- t( mvtnorm::rmvnorm( s, mean  = rep( 0, nrow( R ) ),
                                          sigma = c*R )
                       )
      }

      # Generate population sample
      # ( 3D array with i pop matrices with psample x traits dimensions )

      P             <- array( dim = c( psample, mtraits, reps ) )
      dimnames( P ) <- list( paste( "Sp_", seq( 1:( psample ) ), sep="" ),
                             paste( "trait_", seq( 1:( mtraits ) ), sep="" ),
                             paste( "rep_", seq( 1:( reps ) ), sep="" )
                            )

      for ( i in seq( 1:reps ) ){
        P[ , , i] <- t( mvtnorm::rmvnorm( psample, mean = rep( 0, nrow( R ) ),
                                          sigma = c*R )
                       )

      }

    }


    # C) If the correlation matrix is not provided nor parameter rho ...

    else if ( is.null( rho ) & is.null( R ) ){

      cat( "\n", "You have not provided rho nor R. The data simulation will not account for correlation", "\n")

      # Sample num.species*traits samples from a multiv.norm.
      # distrib with mu = 0 and sd = c
      # ( Array with i noise matrices with s x traits dimensions )


      N             <- array( dim = c( s, mtraits, reps ) )
      dimnames( N ) <- list( paste( "Sp_", seq( 1:( s ) ), sep="" ),
                             paste( "trait_", seq( 1:( mtraits ) ), sep="" ),
                             paste( "rep_", seq( 1:( reps ) ), sep="" )
      )

      for ( i in seq( 1:reps ) ){
        N[ , , i] <- t( replicate( s, rnorm( mtraits, mean = 0, sd = c ) ) )
      }

      # Generate population sample
      # ( 3D array with i pop matrices with psample x traits dimensions )

      P             <- array( dim = c( psample, mtraits, reps ) )
      dimnames( P ) <- list( paste( "Sp_", seq( 1:( psample ) ), sep="" ),
                             paste( "trait_", seq( 1:( mtraits ) ), sep="" ),
                             paste( "rep_", seq( 1:( reps ) ), sep="" )
                            )

      for ( i in seq( 1:reps ) ){
        P[ , , i] <- t( replicate( psample, rnorm( mtraits, mean = 0, sd = c ) ) )
      }

    }


    # Generate noisy matrix (3D array with i reps)

    M.n <- M + N

    # Calculate var-covar of P (3D array with i reps)
    # and get variances (matrix with reps x mtraits dimensions)

    varcov.P             <- array( dim = c( mtraits, mtraits, reps ) )
    dimnames( varcov.P ) <- list( paste( "trait_", seq( 1:( mtraits ) ), sep="" ),
                                  paste( "trait_", seq( 1:( mtraits ) ), sep="" ),
                                  paste( "rep_", seq( 1:( reps ) ), sep="" )
                                 )

    var.P <- matrix( nrow = reps, ncol = mtraits )
    rownames( paste( "rep_", seq( 1:( reps ) ), sep="" ) )
    colnames( paste( "trait_", seq( 1:( mtraits ) ), sep="" ) )

    for ( i in seq( 1:reps ) ){
      varcov.P[ , , i] <- cov( P[ , , i] )
      var.P[i, ] <- diag( varcov.P[ , , i] )
    }

    # Scale matrix M.n (3D array with i reps)

    cat( "\n", "Scaling the data accounting for population variance ... ", "\n\n" )

    M.s             <- array( dim = c( s, mtraits, reps ) )
    dimnames( M.s ) <- list( paste( "Sp_", seq( 1:( s ) ), sep="" ),
                             paste( "trait_", seq( 1:( mtraits ) ), sep="" ),
                             paste( "rep_", seq( 1:( reps ) ), sep="" )
                            )

    for ( i in seq( 1:reps ) ){
      M.s[ , , i] <- M.n[ , , i] %*% diag( 1 / sqrt( var.P[i, ] ) )
    }

    # If a correlation matrix was provided, tranform M.s

    if ( !is.null( R ) | !is.null( rho ) ){

      # Calculate shrunk correlation matrix

      R.shrunk             <- array( dim = c( mtraits, mtraits, reps ) )
      dimnames( R.shrunk ) <- list( paste( "trait_", seq( 1:( mtraits ) ), sep="" ),
                                    paste( "trait_", seq( 1:( mtraits ) ), sep="" ),
                                    paste( "rep_", seq( 1:( reps ) ), sep="" )
                                  )

      for ( i in seq( 1:reps ) ){

        R.shrunk[ , , i] <- as.matrix( corpcor::cor.shrink( P[ , , i] ),
                                       verbose = F )

      }

      # Check method to decompose inverse of R.shrunk
      # is provided

      if ( missing(method) ){
        stop( "Please select a method to decompose the shrunk correlation matrix,
              either method = \"chol\" or method = \"eigen\" " )
      }

      # Create empty array for transformed data

      Z             <- array( dim = c( s, mtraits, reps ) )
      dimnames( Z ) <- list( paste( "Sp_", seq( 1:( s ) ), sep="" ),
                             paste( "trait_", seq( 1:( mtraits ) ), sep="" ),
                             paste( "rep_", seq( 1:( reps ) ), sep="" )
                           )

      # Match argument

      method <- match.arg(method)

      if ( method == "chol" ){

        cat( "\n", "Transforming the data using the Cholesky decomposition in order to account for correlation ... ", "\n" )

        for ( i in seq( 1:reps ) ){

          # chol returns upper triangular matrix

          U <- chol( R.shrunk[ , , i] )
          L <- t( U )

          # R.shrunk = L %*% U = L %*% t( L )
          # solve( R.shrunk ) = t( solve(L) ) %*% solve( L )
          # solve( R.shrunk ) = t( A )        %*% A
          #
          # all.equal( t( solve( L ) ) %*% solve( L ), solve( R.shrunk ) )

          Linv <- solve( L )

          # Transform scaled data,  M.s, so
          # Z = M.s %*% t( A ) = M.s %*% t( Linv )

          Z[ , , i] <- M.s[ , , i] %*% t( Linv )

        }

      }

      else if ( method == "eigen" ){

        cat( "\n", "Transforming the data using the eigen decomposition in order to account for correlation ... ", "\n" )

        for ( i in seq( 1:reps ) ){

          # eigen returns a list with the eigenvectors
          # and the eigenvalues
          #
          # solve( R.shrunk ) = R.shrunk.inv = t( A ) %*% A
          #
          # t( A ) = V %*% D
          # t( A ) = eigen( R.shrunk.inv )$vectors %*% diag( sqrt( eigen( R.shrunk.inv )$values ) )
          #
          # all.equal( tA %*% t( tA ), solve( R.shrunk ) )

          R.shrunk.inv <- solve( R.shrunk[ , , i] )
          tA <- eigen( R.shrunk.inv )$vectors %*% diag( sqrt( eigen( R.shrunk.inv )$values ) )

          # Transform scaled data, M.s, so
          # Z = M.s %*% t( A )

          Z[ , , i] <- M.s[ , , i] %*% tA

        }

      }

      # If the user wants a phylip alignment output ...

      if ( !is.null( out ) ){

        cat( "\n", "Writing output file/s ... ", "\n" )

        if ( !is.null( names ) ){
          names = names
        }
        else{
          names = NULL
        }

        if ( !is.null( ages ) ){
          ages    = ages
          max.age = 1
        }
        else{
          ages    = NULL
          max.age = NULL
        }

        for ( i in seq( 1: reps ) ){

            write_morpho( filename = paste( out, "_", i, sep = "" ),
                          proc     = Z[ , , i],
                          names    = names,
                          ages     = ages,
                          max.age  = max.age,
                          R.sh     = R.shrunk[ , , i],
                          scaled   = TRUE )
        }

      }


      # Return a list with the simulated data, noise matrix,
      # sample population matrix, noisy matrix, scaled matrix,
      # shrunk matrix, and transformed matrix

      cat( "\n" )
      return( list( M   = M,   N   = N,   P    = P,
                    M.n = M.n, M.s = M.s, R.sh = R.shrunk,
                    Z   = Z )
      )

    }

    # Otherwise, do not transform M.s and return the following ...

    else{

      # If the user wants a phylip alignment output ...

      if ( !is.null( out ) ){

        cat( "\n", "Writing output file/s ... ", "\n" )

        if ( !is.null( names ) ){
          names = names
        }
        else{
          names = NULL
        }

        if ( !is.null( ages ) ){
          ages    = ages
          max.age = 1
        }
        else{
          ages    = NULL
          max.age = NULL
        }

        for ( i in seq( 1: reps ) ){

          write_morpho( filename = paste( out, "_", i, sep = "" ),
                        proc     = M.s[ , , i],
                        names    = names,
                        ages     = ages,
                        max.age  = max.age,
                        scaled   = TRUE
                      )
        }

      }

      # Return a list with the simulated data, noise matrix,
      # sample population matrix, noisy matrix, and scaled matrix

      cat( "\n" )
      return( list( M = M, N   = N,
                    P = P, M.n = M.n, M.s = M.s )
             )
    }

  }

  #\\ If population variance is not provided ...

  else{

    # Do not scale nor transform the data (no Z, keep M as it is)

    cat( "\n", "This simulation does not account for population variance nor correlation", "\n" )

    # If the user wants a phylip alignment output ...

    if ( !is.null( out ) ){

      cat( "\n", "Writing output file/s ... ", "\n" )

      if ( !is.null( names ) ){
        names   = names
      }
      else{
        names = NULL
      }

      if ( !is.null( ages ) ){
        ages = ages
        max.age = 1
      }
      else{
        ages    = NULL
        max.age = NULL
      }

      for ( i in seq( 1: reps ) ){

        write_morpho( filename = paste( out, "_", i, sep = "" ),
                      proc     = M[ , , i],
                      names    = names,
                      ages     = ages,
                      max.age  = max.age
                     )
      }

    }

    # Return a list with the simulated data

    cat( "\n" )
    return( list( M = M) )

  }


}



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







